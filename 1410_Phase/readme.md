How to run through the tutorial
===============================
The code for this tutorial has been developed in GNU Octave a high level programming environment for linear algebra, signal processing and more. Octave is matlab compatible so this tutorial should also work in matlab should you have that available.

Expanded Tutorial Paper
-----------------------
A pdf of the expanded tutorial paper is [here](purves_tutorial_long.pdf) this includes plots of seismis attributes on sections are well as the trace plots included in the TLE tutorial.

Installation
------------
To install GNU octave, 
 - on Linux or MAC follow the instructions for your operating system on the [Octave website](http://www.gnu.org/software/octave/download.html) to get the appropriate packages. NOTE: some Linux distributions come with octave pre-installed.
 - on Windows download the latest binary installers distribution from [Octave Forge](http://sourceforge.net/projects/octave/files/Octave%20Windows%20binaries/)
 - on MacOS_X follow this [guide](http://wiki.octave.org/Octave_for_MacOS_X)

Octave Basics
-------------
If you are familar with Matlab, you already know (mostly) how to use Octave. Otherwise get started with this [cheatsheet](http://ais.informatik.uni-freiburg.de/teaching/ss13/robotics/exercises/cheatsheet.pdf) or dive in with [this one](http://folk.ntnu.no/joern/itgk/refcard-a4.pdf).

Generating nice plots
---------------------
I have recently updated the code to plot each graph / image in its own figure window. Octave will open a new figure at its default size, often not great for the graph in question. However, you can now resize the figure window to a better size and aspect ratio and then hit 'e' to replot the graph.

If you want to save the plots to image files, don't take screenshots instead use [print](http://wiki.octave.org/Octave_for_MacOS_X)

Running the Tutorial Code
-------------------------
Open octave from the command line or start menu shortcut then carry out the following steps on the octave command prompt.
 1. cd into the tutorials directory e.g. cd home\seg\tutorials\1410_Phase
 1. then run the following scripts inside Octave (tip: use 'tab' to autocomplete long names)
 1. configure_tutorial
 1. [plot_complex_attributes_on_a_trace](plot_complex_attributes_on_a_trace.m) - to produce figure 1 from the tutorial. This uses fftshifter.m to compute the hilbert transform, and subsequently envelope and phase.
 1. [plot_complex_attributes_on_a_slice](plot_complex_attributes_on_a_slice.m) - will repeat the trace based computation for each trace within a the section of data. SegyMAT is used to load the seismic section from the .sgy file included in the data folder, and we simply loop over the traces in the resulting 2d array.
 1. [plot_phase_at_envelope_peaks_on_a_trace](plot_phase_at_envelope_peaks_on_a_trace.m) - will produce figure 2 from the tutorial. This repeats the envelope and phase computation and then detects peaks in the envelope function. Envelope and phase are extracted at those locations. Phase values are then plotted against an computed "idealised phase".
 1. [plot_phase_at_envelope_peaks_on_a_section](plot_phase_at_envelope_peaks_on_a_slice.m) - applies the same algorithm from the last step to each trace in a section, creating new sparse envelope and phase sections. Note these could also be computed over all traces in a cube, however soem limiations of octave prevent us from doing so here.

Conclusions from working with Octave
------------------------------------
Octave provides an alternative environment to the commercial Matlab package. The degree of similarity and compatibility is good, certainly for ore functions. However, there are limitations in Octave that make working with it challenging, the most significant being the lack of precompiled 64 bit builds, which limits the amount of data that we can load and process to 2GB. In this tutorial, we've worked with traces and seismic sections that in this case easily fit within that limit. 

However, if you want to work with 3D data or generate maps over hilbert transformed data, then you'll need to put significant effort into managing memory and its likely your processing will be unnecessarily slow as a result, makign withing with 3D in Octave prohibitive.

Note on Plotting: Octave doesn't always autoscale the contents of its figure windows. Sometimes to get nice large/clear plots, you'll need to plot, resize the window and then plot again.

References
----------
 1. For convenience this repositoty includes the source code for [SegyMAT](http://segymat.sourceforge.net/)
 1. Data used in the tutorial has been taken from the Pensobscot survey available in the [OpenDTect Open Seismic Repository](https://opendtect.org/osr/pmwiki.php/Main/PENOBSCOT3DSABLEISLAND)
