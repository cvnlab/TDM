# TDM (Temporal Decomposition through Manifold Fitting)

TDM is a MATLAB toolbox that implements a temporal-decomposition technique
for separating microvasculature- and macrovasculature-related responses
in fMRI time-series data. It is developed by Kendrick Kay (kendrick@post.harvard.edu).

The toolbox has the following dependencies. You must download (or clone)
these repositories and ensure that they are on your path.
- http://github.com/kendrickkay/knkutils
- http://github.com/kendrickkay/GLMdenoise

(In actuality, the GLMdenoise repository is not necessary for the core of the 
TDM technique, which is extracthrfmanifold.m. However, the example1.m
script makes use of the GLMdenoise toolbox to perform GLM fitting,
so you probably want to get the repository.)

The toolbox also relies on several external dependencies, a copy of which are
included in the folder "external". These external dependencies include:
- https://www.github.com/AntonSemechko/S2-Sampling-Toolbox
- https://research.ics.aalto.fi/ica/fastica/
(The fastica toolbox is necessary only if you want to use the ICA-based procedure.)

The total install time should be just a few minutes.

To use the toolbox, add it to your MATLAB path:
  addpath(genpath('TDM'));

To try the toolbox on an example dataset, change to the TDM directory and then type:
  example1;
This script will download the example dataset (if it has not already been
downloaded) and will go through the example analyses. The runtime of the 
analysis will likely be no more than about 5-10 minutes.

Outputs from the example script (via MATLAB's publish) are available here:
https://htmlpreview.github.io/?https://github.com/kendrickkay/TDM/blob/master/html/example1.html
https://htmlpreview.github.io/?https://github.com/kendrickkay/TDM/blob/master/html/example2.html

A video walk-through of the example1.m script is viewable here:
  https://www.youtube.com/watch?v=Sz13i-9EtmA
The video is quite comprehensive in what it covers and discusses.

Terms of use: This content is licensed under a BSD 3-Clause License.

If you use TDM in your research, please cite the following pre-print:
  Kay, K., Jamison, K.W., Zhang, R.-Y., Ugurbil, K. 
    A temporal decomposition method for identifying venous effects in task-based fMRI.
    Nature Methods (2020).

Additional information can be found at the OSF site https://osf.io/j2wsc/.

The code has been successfully run on MATLAB R2018a, but is likely compatible with
any recent version of MATLAB.

History of major code changes:
- 2020/03/30 - Version 1.1. See extracthrfmanifold.m for description of changes.
- 2019/12/08 - Initial version.

## CONTENTS

Contents:
- derivehrf.m - Use PCA to extract an HRF common to a set of timecourses
- example*.m - Example scripts
- external - External toolboxes included for your convenience
- extracthrfmanifold.m - The TDM technique
- html - Outputs from MATLAB's publish
- icadecomposehrf.m - Alternative ICA-based technique
- LICENSE - license information
- README.md - The file you are reading
- setup.m - A simple script that downloads the example data files and adds TDM to the MATLAB path
