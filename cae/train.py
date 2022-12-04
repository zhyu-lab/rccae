import argparse
import os
import datetime
import torch
from torch.optim import lr_scheduler
import numpy as np
from matplotlib import pyplot as plt
from sklearn.metrics.cluster import v_measure_score, adjusted_rand_score
from sklearn.decomposition import PCA
import autoencoder
from autoencoder import AE
from genomedata import GenomeData
from clustering import GMM

os.environ["CUDA_VISIBLE_DEVICES"] = "0"

def main(args):
    start_t = datetime.datetime.now()

    gd = GenomeData(args.input)
    gd.load_data()
    gd.preprocess_data()

    data = gd.data_lrc_all
    data = np.transpose(data)
    data_bk = data.copy()
    print(data_bk.shape)
    data = data_bk.copy()

    def setup_seed(seed):
        torch.manual_seed(seed)
        torch.cuda.manual_seed_all(seed)
        np.random.seed(seed)

    setup_seed(args.seed)

    # define and create AE architecture
    ae_model = AE(data.shape[1], args.latent_dim, args.kernel_size).cuda()

    epochs = []
    train_loss = []
    test_acc = []

    optimizer = torch.optim.Adam(ae_model.parameters(), lr=args.lr)

    # Start training the model
    ae_model.train()
    for epoch in range(args.epochs):
        print("epoch", epoch)
        epochs.append(epoch)
        total_loss = 0
        for step, x in autoencoder.xs_gen(data, args.batch_size, 1):
            x = torch.from_numpy(x).cuda()
            x = x.unsqueeze(1)

            z, y = ae_model(x)
            x = x.squeeze()
            loss = autoencoder.mse_loss(y, x)
            total_loss += loss
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()

        print("epoch:", epoch, "| train loss: %.4f" % total_loss.data.cpu().numpy())
        train_loss.append(total_loss.data.cpu().numpy())

    output_dir = args.output
    if not os.path.isdir(output_dir):
        os.mkdir(args.output)
    ll_file = output_dir + '/loss.txt'
    if os.path.isfile(ll_file):
        os.remove(ll_file)
    file_o = open(ll_file, 'w')
    np.savetxt(file_o, np.c_[np.reshape(train_loss, (1, len(train_loss)))], fmt='%f', delimiter=',')
    file_o.close()

    # get latent representation of single cells after CAE training is completed
    z_hiddn = []
    x_cst = []
    data = data_bk.copy()

    ae_model.eval()
    for step, x in autoencoder.xs_gen(data, args.batch_size, 0):
        x = torch.from_numpy(x).cuda()
        x = x.unsqueeze(1)
        with torch.no_grad():
            z, y = ae_model(x)
            z = z.cpu().detach().numpy()
            y = y.cpu().detach().numpy()
            z_hiddn.append(z)
            x_cst.append(y)

    data_lrc = []
    data_hidden = []
    for i, z in enumerate(z_hiddn):
        if i == 0:
            features = z
        else:
            features = np.r_[features, z]
    for i, y in enumerate(x_cst):
        if i == 0:
            data_lrc = y
        else:
            data_lrc = np.r_[data_lrc, y]

    # use Gaussian mixture model to cluster the single cells
    print('clustering the cells...')
    if args.max_k <= 0:
        max_k = np.max([1, features.shape[0] // 5])
    else:
        max_k = np.min([args.max_k, features.shape[0] // 5])
    label_p, K = GMM(features, max_k).cluster()
    print('inferred number of clusters: {}'.format(K))

    # save results
    output_dir = args.output
    rc_file = output_dir + '/lrc.txt'
    if os.path.isfile(rc_file):
        os.remove(rc_file)
    file_o = open(rc_file, 'a')
    np.savetxt(file_o, np.c_[gd.bin_size], fmt='%d', delimiter=',')
    np.savetxt(file_o, np.c_[np.reshape(gd.data_chr_all, (1, len(gd.data_chr_all)))], fmt='%d', delimiter=',')
    np.savetxt(file_o, np.c_[np.reshape(gd.data_bin_all, (1, len(gd.data_bin_all)))], fmt='%d', delimiter=',')
    np.savetxt(file_o, np.c_[data_lrc], fmt='%.3f', delimiter=',')
    file_o.close()

    label_file = output_dir + '/label.txt'
    file_o = open(label_file, 'w')
    np.savetxt(file_o, np.c_[np.reshape(gd.barcodes, (1, len(gd.barcodes)))], fmt='%s', delimiter=',')
    np.savetxt(file_o, np.c_[np.reshape(label_p, (1, len(label_p)))], fmt='%d', delimiter=',')
    file_o.close()

    latent_file = output_dir + '/latent.txt'
    file_o = open(latent_file, 'w')
    np.savetxt(file_o, np.c_[features], fmt='%.3f', delimiter=',')
    file_o.close()

    end_t = datetime.datetime.now()
    print('elapsed time: ', (end_t-start_t).seconds)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="rcCAE")
    parser.add_argument('--epochs', type=int, default=100, help='number of epoches to train the CAE.')
    parser.add_argument('--batch_size', type=int, default=64, help='batch size.')
    parser.add_argument('--lr', type=float, default=0.0001, help='learning rate.')
    parser.add_argument('--max_k', type=int, default=0, help='the maximum number of clusters to consider.')
    parser.add_argument('--latent_dim', type=int, default=3, help='the latent dimension.')
    parser.add_argument('--kernel_size', type=int, default=7, help='convolutional kernel size.')
    parser.add_argument('--seed', type=int, default=0, help='random seed.')
    parser.add_argument('--input', type=str, default='', help='a file containing read counts, GC-content and mappability data.')
    parser.add_argument('--output', type=str, default='', help='a directory to save results.')
    args = parser.parse_args()
    main(args)
