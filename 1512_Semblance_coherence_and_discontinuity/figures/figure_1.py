import numpy as np
import matplotlib.pyplot as plt

import basic_methods
from basic_methods import moving_window
from parent_directory import data, image_dir


def main():
    seismic = data.load_seismic()

    fig, axes = setup_figure()

    plot(axes[:,0], 'Seismic Data', '', seismic, 'A')
    plot(axes[:,1], 'Bahorich & Farmer (1995)', 'Cross-correlation',
         basic_methods.bahorich_coherence(seismic, 21), 'B')
    plot(axes[:,2], 'Marfurt, et al, (1998)', 'Semblance-based',
         moving_window(seismic, (3,3,9), basic_methods.marfurt_semblance), 'C')
    plot(axes[:,3], 'Gersztenkorn and Marfurt (1999)', 'Eigenstructure-based',
         moving_window(seismic, (3,3,9), basic_methods.eig), 'D')
    plot(axes[:,4], 'Dip Corrected', 'Eigenstructure-based',
         dip_corrected(seismic, (3, 3, 9), basic_methods.eig), 'E')

    fig.savefig(image_dir('figure_1.png'), dpi=200, bbox_inches='tight')
    plt.show()

def dip_corrected(seismic, window, func):
    surface = data.load_horizon()
    flat = basic_methods.flatten(seismic, surface, seismic.shape[-1])
    sembl = moving_window(flat, window, func)
    return basic_methods.unflatten(sembl, surface, seismic.shape)

def plot(axes, title, subtitle, data, letter=''):
    data = np.ma.masked_equal(data, 0)
    axes[0].set_title(title, size=12)
    axes[1].set_xlabel(subtitle, size=12)

    j0 = data.shape[1] // 2
    k0 = data.shape[2] // 2
    axes[0].imshow(data[:,:,k0].T, cmap='gray')
    axes[1].imshow(data[:,j0,:].T, cmap='gray')

    axes[0].annotate(letter, xy=(0.5, 1), xytext=(0, 20),
                     textcoords='offset points', xycoords='axes fraction',
                     ha='center', va='bottom', size=14)

def setup_figure():
    fig = plt.figure(figsize=(13,5))
    gs = plt.GridSpec(4, 5, left=0.05, right=0.95, hspace=0.05)

    axes = np.empty((2, 5), dtype=object)
    for col in range(5):
        axes[0,col] = fig.add_subplot(gs[:3, col])
        axes[0,col].set(anchor='S')
    for col in range(5):
        axes[1,col] = fig.add_subplot(gs[-1, col])
        axes[1,col].set(anchor='N')

    for ax in axes.flat:
        ax.set(xticks=[], yticks=[])

    return fig, axes

main()
