"""
Note that this is a very hackish script to put together this figure.
Forgive the sloppy approach.
"""
import numpy as np
import scipy.ndimage
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

import basic_methods
from parent_directory import image_dir

def main():
    fig, axes = setup_figure()

    functions = [in_phase_traces, out_of_phase_traces]
    labels = ['Different Amplitude', 'Shifted Traces']
    for i, (func, label) in enumerate(zip(functions, labels)):
        row = axes[i]
        wiggles(row[0], *func())
        cross_corr(row[1], func)

        for ax in row[2:]:
            setup_trace_xy(ax, *func())

        semblance(row[2], func)
        eigenstructure(row[3], func)

        row[0].set_ylabel(label, size=14)

    # Locatable axes with shared axis's don't observe the adjustable param??
    # Setting limits manually... Shouldn't be necessary...
    for ax in axes[0,2:]:
        ax.axis([-2.1, 2.1, -2.1, 2.1])
    for ax in axes[1,2:]:
        ax.axis([-1.1, 1.1, -1.1, 1.1])

    axes[0, 0].set_title('Input Traces')
    bottom_label(axes[1, 0], 'A')
    axes[0, 1].set_title('Bahorich & Farmer (1995)\nCross Correlation')
    bottom_label(axes[1, 1], 'B')
    axes[0, 2].set_title('Marfurt, et al (1998)\nSemblance')
    bottom_label(axes[1, 2], 'C')
    axes[0, 3].set_title('Gersztenkorn & Marfurt (1999)\nEigenstructure')
    bottom_label(axes[1, 3], 'D')

    fig.savefig(image_dir('figure_2.png'), dpi=200, bbox_inches='tight')

def generate_t(longer=False):
    t = np.linspace(0, 10*np.pi+2, 10)
    if longer:
        dt = t[1] - t[0]
        t = np.r_[np.arange(-4, 0) * dt + t[0], t, np.arange(1, 5) * dt + t[-1]]
    return t

def identical_traces(t=None):
    t = generate_t() if t is None else t
    y = np.sin(t + np.pi/8)
    return t, y, y

def in_phase_traces(t=None):
    t = generate_t() if t is None else t
    y = np.sin(t + np.pi/8)
    y2 = 2 * y
    return t, y, y2

def out_of_phase_traces(t=None):
    t = generate_t() if t is None else t
    y = np.sin(t + np.pi/8)
    y2 = 1 * np.sin(t + 2*np.pi/4)
    return t, y, y2

def setup_figure():
    fig = plt.figure(figsize=(12, 6))
    gs = plt.GridSpec(2, 14, left=0.05, right=0.95, hspace=0.05)

    rows = []
    for i in range(2):
        cols = [np.s_[:2], np.s_[2:5], np.s_[5:10], np.s_[10:]]
        row = [fig.add_subplot(gs[i, col]) for col in cols]
        rows.append(row)

    rows = np.array(rows, dtype=object)
    for ax in rows.flat:
        ax.set(xticks=[], yticks=[])

    return fig, rows

def wiggles(ax, t, y1, y2):
    shift = y1.max() - y2.min() #+ 0.1 * y2.ptp()
    y2 = y2 + shift

    ax.plot(y1, t, 'o', zorder=4, color='lightblue')
    ax.plot(y2, t, 'o', zorder=4, color='lightblue')

    t = scipy.ndimage.zoom(t, 30)
    y1 = scipy.ndimage.zoom(y1, 30)
    y2 = scipy.ndimage.zoom(y2, 30)

    ax.fill_betweenx(t, y1, where=y1 > 0, facecolor='black')
    ax.plot(y1, t, color='black')

    ax.fill_betweenx(t, y2, shift, where=y2 > shift, facecolor='black')
    ax.plot(y2, t, color='black')

    ax.margins(x=0.05)

def semblance(ax, func):
    _, y1, y2 = func()
    c = basic_methods.semblance2(np.vstack([y1, y2]))
    label_coherence(ax, 'C={:0.2f}'.format(c))
    x = np.linspace(-10, 10, 2)
    ax.plot(x, x, color='black', zorder=0, scalex=False, scaley=False)

def eigenstructure(ax, func):
    _, y1, y2 = func()
    c = basic_methods.eig(np.vstack([y1, y2]))
    label_coherence(ax, 'C={:0.2f}'.format(c))
    x = np.linspace(-10, 10, 2)
    a, _, _, _ = np.linalg.lstsq(y1[:,None], y2)
    ax.plot(x, a*x, color='black', zorder=0, scalex=False, scaley=False)


def cross_corr(ax, trace_func):
    tshort = generate_t()
    tlong = generate_t(longer=True)
    y1 = trace_func(tlong)[1]
    y2 = trace_func(tshort)[2]

    xcorr = np.correlate(y1, y2, mode='same')
    xcorr /= (y1.std() * y2.std() * y2.size)

    ax.plot(xcorr, color='black')
    ax.plot([xcorr.argmax()], [xcorr.max()], 'o', color='red')
    ax.axhline(0, color='gray', zorder=0)
    ax.axvline(xcorr.size // 2, color='gray', zorder=0)
    label_coherence(ax, 'C={:0.1f}'.format(xcorr.max()))

    ax.margins(y=0.2)

def setup_trace_xy(ax, t, y1, y2):
    ax.plot(y1, y2, 'o', color='lightblue')
    ax.set(aspect=1, adjustable='datalim')

    divider = make_axes_locatable(ax)

    hax = divider.append_axes('bottom', size='40%', pad=0, sharex=ax)
    vax = divider.append_axes('left', size='40%', pad=0, sharey=ax)

    for ax in [hax, vax]:
        ax.set(xticks=[], yticks=[], adjustable='datalim')
    ax.figure.add_axes(hax)
    ax.figure.add_axes(vax)

#    hax.set_xlabel('Trace 1 Amplitude')
#    vax.set_ylabel('Trace 2 Amplitude')

    vax.plot(t, y2, marker='o', color='black', mfc='lightblue')
    tt, y2 = scipy.ndimage.zoom(t, 30), scipy.ndimage.zoom(y2, 30)
    vax.fill_between(tt, y2, where=y2 > 0, facecolor='black')

    hax.plot(y1, t, marker='o', color='black', mfc='lightblue')
    tt, y1 = scipy.ndimage.zoom(t, 30), scipy.ndimage.zoom(y1, 30)
    hax.fill_betweenx(tt, y1, where=y1 > 0, facecolor='black')

    hax.margins(0.05)
    vax.margins(0.05)
    return hax, vax

def label_coherence(ax, label):
    ax.annotate(label, xy=(0,1), xytext=(5, -5),
                xycoords='axes fraction', textcoords='offset points',
                ha='left', va='top')

def bottom_label(ax, label):
    ax.annotate(label, xy=(0.5, 0), xytext=(0, 14),
                xycoords=('axes fraction', 'figure fraction'),
                textcoords='offset points',
                ha='center', va='bottom', size=14)

if __name__ == '__main__':
    main()
    plt.show()
