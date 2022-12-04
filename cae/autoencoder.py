from matplotlib import pyplot as plt
import torch.nn as nn
import torch
import numpy as np
import math
from mpl_toolkits.axes_grid1 import host_subplot


def initialize_weights(net):
    for m in net.modules():
        if isinstance(m, nn.Conv2d):
            m.weight.data.normal_(0, 0.02)
            m.bias.data.zero_()
        elif isinstance(m, nn.ConvTranspose2d):
            m.weight.data.normal_(0, 0.02)
            m.bias.data.zero_()
        elif isinstance(m, nn.Linear):
            m.weight.data.normal_(0, 0.02)
            m.bias.data.zero_()


# Define loss function
def mse_loss(predict, target):
    loss = torch.nn.MSELoss()
    return loss(predict, target)


def xs_gen(cell_list, batch_size, random):
    if random == 1:
         np.random.shuffle(cell_list)
    steps = math.ceil(len(cell_list) / batch_size)
    for i in range(steps):
        batch_x = cell_list[i * batch_size: i * batch_size + batch_size]
        yield i, batch_x


class AE(torch.nn.Module):
    """
    This class implements a convolutional autoencoder
    """
    def __init__(self, in_dim, z_dim=3, k_size=7, device='cuda'):
        super(AE, self).__init__()
        self.device = device

        d = in_dim
        for i in range(3):
            d = np.floor((d - k_size) / 1 + 1)
        d = np.int32(d)

        self.encoder = nn.Sequential(
            nn.Conv1d(1, 128, k_size, stride=1),
            nn.LeakyReLU(),
            nn.Conv1d(128, 64, k_size, stride=1),
            nn.LeakyReLU(),
            nn.Conv1d(64, 32, k_size, stride=1),
            nn.LeakyReLU()
        )

        self.fc1 = nn.Linear(d * 32, z_dim)

        self.fc2 = nn.Linear(z_dim, d * 32)
        self.decoder = nn.Sequential(
            nn.LeakyReLU(),
            nn.ConvTranspose1d(32, 64, k_size, stride=1),
            nn.LeakyReLU(),
            nn.ConvTranspose1d(64, 128, k_size, stride=1),
            nn.LeakyReLU(),
            nn.ConvTranspose1d(128, 1, k_size, stride=1),
        )

        initialize_weights(self)

    def bottleneck(self, h):
        z = self.fc1(h)
        return z

    def encode(self, x):
        h = self.encoder(x)
        # print(h.size())
        h = h.view(h.size(0), -1)
        z = self.bottleneck(h)
        return z

    def decode(self, z):
        z = self.fc2(z)
        z = z.view(z.size(0), 32, int(z.size(1) / 32))
        # print(z.size())
        z = self.decoder(z)
        z = z.squeeze()
        # print(z.size())
        return z

    def forward(self, x):
        z = self.encode(x)
        y = self.decode(z)
        return z, y

