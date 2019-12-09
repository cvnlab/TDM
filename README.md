# TDM (Temporal Decomposition through Manifold Fitting)

TDM is a MATLAB toolbox that implements a temporal-decomposition technique
for separating microvasculature- and macrovasculature-related responses
in fMRI time-series data. It is developed by Kendrick Kay (kendrick@post.harvard.edu).

The toolbox has the following dependencies. You must download (or clone)
these repositories and ensure that they are on your path.
- http://github.com/kendrickkay/knkutils
- http://github.com/kendrickkay/GLMdenoise

(In actuality, the GLMdenoise repository is not necessary for the core of the 
TDM technique, which is extracthrfmanifold.m. However, the example dataset
and example script make use of the GLMdenoise toolbox to perform GLM fitting,
so you probably want to get the repository.)

The toolbox also relies on several external dependencies, a copy of which are
included in the folder "external". These external dependencies include:
- https://www.mathworks.com/matlabcentral/fileexchange/66629-2-d-histogram-plot
- https://www.github.com/AntonSemechko/S2-Sampling-Toolbox
- https://research.ics.aalto.fi/ica/fastica/

(In actuality, the fastica toolbox is necessary only if you want to use
the ICA-based procedure. An example of this is done in the example1.m script.)

To use the toolbox, add it to your MATLAB path:
  addpath(genpath('TDM'));

To try the toolbox on an example dataset, change to the TDM directory and then type:
  example1;
This script will download the example dataset (if it has not already been
downloaded) and will go through the example analyses.

Outputs from the example script (via MATLAB's publish) are available here:
https://htmlpreview.github.io/?https://github.com/kendrickkay/TDM/blob/master/html/example1.html

A video walk-through of the example script is viewable here:
  https://www.youtube.com/watch?v=Sz13i-9EtmA

Terms of use: This content is licensed under a Creative Commons Attribution 3.0 
Unported License (http://creativecommons.org/licenses/by/3.0/us/). You are free 
to share and adapt the content as you please, under the condition that you cite 
the appropriate manuscript (see below).

If you use TDM in your research, please cite the following pre-print:
  Kay K, Jamison KW, Zhang R, and Ugurbil K (2019)
    TDM: a temporal decomposition method for removing venous effects from task-based fMRI.
    bioRxiv. 

History of major code changes:
- 2019/12/08 - Initial version.

## CONTENTS

Contents:
- derivehrf.m - Use PCA to extract an HRF common to a set of timecourses
- example*.m - Example scripts
- external - External toolboxes included for your convenience
- extracthrfmanifold.m - The TDM technique
- html - Outputs from MATLAB's publish
- icadecomposehrf.m - Alternative ICA-based technique
- README.md - The file you are reading
- setup.m - A simple script that downloads the example dataset and adds TDM to the MATLAB path

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Copyright (c) 2019, Kendrick Kay
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.

Redistributions in binary form must reproduce the above copyright notice, this
list of conditions and the following disclaimer in the documentation and/or
other materials provided with the distribution.

The names of its contributors may not be used to endorse or promote products 
derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
