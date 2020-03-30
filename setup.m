% This script downloads the example dataset (if it has not already
% been downloaded) and adds TDM to the MATLAB path.

% Download the dataset if it has not already been downloaded.
% (If you like, you may manually download the dataset from:
%    https://osf.io/j2wsc/)
files = {'exampledataset.mat' 'exampledataset_inputsonly.mat'};
osfids = {'d6z8s' '89mfu'};
for p=1:length(files)
  if ~exist(fullfile(pwd,files{p}),'file')
    fprintf('Downloading %s (please be patient).\n',files{p});
    urlwrite(sprintf('https://osf.io/%s/download',osfids{p}),files{p});
    fprintf('Downloading is done!\n');
  end
end
clear p files osfids;

% Add TDM to the MATLAB path (in case the user has not already done so).
path0 = strrep(which('setup.m'),[filesep 'setup.m'],'');
addpath(genpath(path0));
clear path0;
