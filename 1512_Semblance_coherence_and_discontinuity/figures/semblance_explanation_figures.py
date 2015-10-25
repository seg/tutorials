import numpy as np
import matplotlib.pyplot as plt

import figure_2
from parent_directory import image_dir

def main():
    run(figure_2.semblance,
        'Semblance-Based\nMarfurt, et al (1998)',
        'semblance_{}.png')

    run(figure_2.eigenstructure,
        'Eigenstructure-Based\nGersztenkorn & Marfurt (1999)',
        'eigenstructure_{}.png')

def run(coherence_func, title, template):
    trace_funcs = [figure_2.identical_traces, figure_2.in_phase_traces,
                   figure_2.out_of_phase_traces]
    ylabels = ['Identical Traces', 'Different Amplitude', 'Shifted Traces']

    for trace_func, ylabel in zip(trace_funcs, ylabels):
        filename = template.format(ylabel.replace(' ', '_'))
        fig = compare(trace_func, figure_2.semblance, title, ylabel)
        fig.savefig(image_dir(filename), bbox_inches='tight', dpi=80)

def compare(trace_func, coherence_func, title, ylabel):
    fig, axes = setup_subplots()

    figure_2.wiggles(axes[0], *trace_func())
    figure_2.setup_trace_xy(axes[1], *trace_func())
    coherence_func(axes[1], trace_func)

    axes[0].margins(x=0.2, y=0)
    axes[0].set(xticks=[], yticks=[])
    axes[0].set_ylabel(ylabel, size='large')

    axes[1].set(title=title)
    fix_scaling(axes[1], trace_func)
    return fig

def setup_subplots():
    gs = plt.GridSpec(1, 12)
    fig = plt.figure(figsize=(9, 5))

    # Cheating quite a bit by making things overlap.
    # axes_grid + fixed aspect ratio + non-axes_grid leads to some hacks.
    ax1 = fig.add_subplot(gs[:2], anchor='E')
    ax2 = fig.add_subplot(gs[:], anchor='W')
    return fig, np.array([ax1, ax2], dtype='object')

def fix_scaling(ax, trace_func):
    # Locatable axes with shared axis's don't observe the adjustable param??
    # Hack to get around some limitations. Ideally use `ax.margins` instead.
    trace = np.hstack(trace_func()[1:])
    pad = 0.05 * trace.ptp()
    ax.axis([trace.min() - pad, trace.max() + pad,
             trace.min() - pad, trace.max() + pad])

if __name__ == '__main__':
    main()
