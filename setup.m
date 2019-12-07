% This script downloads the example dataset (if it has not already
% been downloaded) and adds TDM to the MATLAB path.

% Download the dataset if it has not already been downloaded.
% (If you like, you may manually download the dataset from:
%    http://kendrickkay.net/TDM/exampledataset.mat)
files = {'exampledataset.mat'};
for p=1:length(files)
  if ~exist(fullfile(pwd,files{p}),'file')
    fprintf('Downloading %s (please be patient).\n',files{p});
    urlwrite(sprintf('http://kendrickkay.net/TDM/%s',files{p}),files{p});
    fprintf('Downloading is done!\n');
  end
end
clear p files;

% Add TDM to the MATLAB path (in case the user has not already done so).
path0 = strrep(which('setup.m'),[filesep 'setup.m'],'');
addpath(genpath(path0));
clear path0;
