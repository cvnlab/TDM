%% Example 2: Apply TDM on an example dataset (quick version)

% This is the same as Example 1 except that we skip straight to the
% part where extracthrfmanifold.m is called. This is achieved by loading
% exampledataset_inputsonly.mat, which contains the inputs to 
% extracthrfmanifold.m as computed from the Example 1 script.

%% Download dataset (if necessary) and add TDM to the MATLAB path

setup;

%% Load in the data

% Load in the data
load('exampledataset_inputsonly.mat','data','intensity','tr');

%% Run extracthrfmanifold.m

% Use extracthrfmanifold.m to derive estimates of the latent early
% and late timecourses. The 1 indicates to create figure windows.
results = extracthrfmanifold(data,intensity,tr,1);
%%
