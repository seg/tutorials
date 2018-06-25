# Contributing to the Geophysical Tutorials

The original proposal, by Tad Smith, editor (from TLE, December 2013):

> I’ve always been a fan of tutorial papers, as they provide the non-expert with a conceptual starting point from which to build an understanding of complex topics. Starting in 2014, TLE is going to start publishing one- to two-page tutorial columns, each of which will conceptually introduce and briefly discuss various aspects of geophysics (or other related disciplines). In this age of technical and electronic media overload, I’m hopeful the addition of a tutorial column will provide our readers with valuable new content that will help them continue to grow as geoscientists and perhaps inspire them to dig deeper into technical issues of interest. 

----

## Instructions to authors

### VISION

The idea is to be relentlessly practical and useful — open, shareable content that helps people.

The tutorials will have the following features:

- A short, easy-to-read introduction.
- One or two explanatory figures with captions.
- A snippet of code and links to other Python or MATLAB code for people to use.
- Everything under an open CC-BY or CC-BY-SA license, or equivalent.

The tutorials will appear in The Leading Edge. It’s possible we will receive others around the same time as yours, in which case we’ll make decisions about what order to run them in.

The tutorials will also appear in SEG Wiki, retaining the author attribution. Since they will be openly licensed, it’s possible the tutorials will appear elsewhere too, but only with author attribution. 


### SUBJECT

Anything you like. The idea is to be helpful to some reasonable number of people, especially those at the beginnings of their career or getting into a subject for the first time. Ask yourself questions like these:

- Why would anyone care about the topic? How is it applied?
- If there was just one key thing to know, what would it be?
- What’s an error you often see people make?
- What do you wish you’d known when you got into the subject?
- What helped you ‘get it’, or helps you conceptualize the problem?
- What are the next 3 or 4 most important things to know? (People like lists)
- What are the 3–5 most important things to read next?
- How can you illustrate a pertinent aspect of your ‘one thing to know’?


### FORMAT

First, read some of the existing tutorials! They have the following features:

- Title, author name, affiliation and contact details.
- NO more than 600 words for a one-pager, or 1000 for a two-pager. The absolute maximum we can accept is 2000 words.
- Include at least one figure, two is preferable, or three for a two- or three-pager.
- Include at least one code snippet or mini application, preferably in Python or MATLAB.
- Use only synthetic or open data, e.g. see http://wiki.seg.org/index.php/Open_data
- Your figures should ideally be at least 2000 pixels wide. Please save them in PDF or PNG format, rather than JPG or anything else. You can send them to us via Dropbox, Google Drive, or Hightail.
- You must be prepared to release some code to illustrate your method or lesson. We can help work out the terms, but a good starting point is the Apache 2.0 license. I have set up a GitHub account for SEG code, or we may use one of SEG's existing sites to share code.
- If you use Git and GitHub, feel free to clone the appropriate repo (look at the year) and add your files there.

Please write in something we can easily transform to MS Word (I know, sorry). For example:

- Open Office
- Google Docs
- Markdown
- LaTeX is OK, but please use Authorea orOverleaf so we can easily cast to Word.

Current best practice is to make the Jupyter Notebook into the actual manuscript. We have ways to hide code blocks and modify the input and output of blocks (using Jupyter Notebook cell 'tags'). This makes for a much smoother editing process. Please see the tutorials in 2018 to understand more about how we do this.


### DEADLINES

| Issue  | Deadline |
|---|---|
| October 2018 | 15 August 2018 |
| December 2018 | 15 October 2018 |
| February 2019 | 15 December 2018 | 
| April 2019 | 15 February 2019 | 
| June 2019 | 15 April 2019 |
| August 2019 | 15 June 2019 |

----

## Ideas for columns

### SEISMIC ACQUISITION

- Designing a land survey geometry
- Designing a marine survey geometry
- A toy ray-tracer, with interactive CUDA demo punchline

### SEISMIC PROCESSING

- 10 things to ask for from your processor
- What is a gather?
- What is noise? — S:N, coloured noise, geologic vs geophysical, coherent vs chaotic vs random

### SIGNAL PROCESSING

- What is frequency? or What is spectral decomposition? (Fourier, CWT, etc)
- What is time-frequency decomposition?
- Spatial resolution or just Seismic resolution — what is it really?

### SEISMIC INTERPRETATION

- Good displays
- Autotracking
- Seismic stratigraphy essentials
- How to wield the wedge model
- What is signal:noise?
- Choosing attributes
- Semblance, coherency, similarity, edge detect, curvature: what’s the diff?
- Spectral decomposition workflow

### SEISMIC ANALYSIS

- How to prove something worked (anecdotes are not enough)
- Starting an inversion project
- What is full waveform inversion
- What is anisotropy?
- Backus averaging for all

### SEISMIC PETROPHYSICS

- Asking for logs from the petrophysicist
- Velocity dispersion
- DT, DTS, and the sonic scanner / How does the sonic scanner work?
- Extracting more attributes from sonic scanner data
- What are checkshots?
- What to do with a VSP
