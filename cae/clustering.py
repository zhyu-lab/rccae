import numpy as np
import matplotlib.pyplot as plt
from sklearn.mixture import GaussianMixture
import operator


class GMM:
    def __init__(self, data, max_k):
        """
        self.X : data to cluster
        self.x,self.y,self.z：first three dimensions of self.X；
        self.max_num_cluster：the maximum number of clusters to consider
        """
        self.X = data
        self.x = self.X[:, 0]
        self.y = self.X[:, 1]
        self.z = self.X[:, 2]
        self.max_num_cluster = max_k

    def cluster(self):
        bic = []
        labels_all = []
        valid = []
        min_bic = []
        n_components = 1
        count = 0
        n_components_range = []
        # Gaussian mixture model is used for clustering, and BIC is used to select the optimal number of clusters
        while n_components <= self.max_num_cluster:
            n_components_range.append(n_components)
            gmm = GaussianMixture(n_components=n_components, covariance_type='full', n_init=50)
            gmm.fit(self.X)
            labels = gmm.predict(self.X)
            labels_u, counts = np.unique(labels, return_counts=True)
            if np.sum(counts < 3) / n_components > 0.2:
                valid.append(False)
            else:
                valid.append(True)
            labels_all.append([])
            labels_all[n_components-1].append(labels)
            s = gmm.bic(self.X)
            bic.append(s)
            if n_components == 1 or s < min_bic:
                min_bic = s
                count = 0
            elif s >= min_bic:
                count += 1
            if count >= 10:
                break
            n_components += 1

        n_components_range = np.array(n_components_range)
        n_components_range = n_components_range[valid]
        print(n_components_range)
        bic = np.array(bic)
        bic = bic[valid]
        labels_all = np.array(labels_all)
        labels_all = labels_all[valid, ]
        ind, val = min(enumerate(bic), key=operator.itemgetter(1))
        num_cluster = n_components_range[ind]
        labels = np.squeeze(labels_all[ind, ])

        # visualize clustering results
        fig = plt.figure("3D Scatter", facecolor="lightgray")
        ax3d = fig.add_subplot(projection="3d")
        plt.title('3D Scatter', fontsize=20)
        ax3d.set_xlabel('x', fontsize=14)
        ax3d.set_ylabel('y', fontsize=14)
        ax3d.set_zlabel('z', fontsize=14)
        plt.tick_params(labelsize=10)
        ax3d.scatter(self.x, self.y, self.z, s=20, c=[labels], cmap='viridis', marker="o")
        plt.show()
        return labels, num_cluster
