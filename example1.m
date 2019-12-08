%% Example 1: Apply TDM on an example dataset

%% Download dataset (if necessary) and add TDM to the MATLAB path

setup;

%% Load in the data

% Load in the data
load('exampledataset.mat');

% Check the workspace
whos
%%

%% Inspect the data

% Check dimensionality of the data
data
%%

% Check dimensionality of the design matrices
design
%%

% Check stimulus duration
stimdur
%%

% Check TR
tr
%%

%% Fit Standard GLM

% Here, we quickly fit a vanilla ("Standard") GLM that uses a single canonical HRF.

% Define options (no bootstrapping, no denoising (in order to speed things up))
opt = struct('numboots',0,'numpcstotry',0);  % to use denoising, could set numpcstotry to 20

% Use GLMdenoisedata.m to fit a GLM that incorporates a single canonical HRF.
% We use 6 condition-splits (so we get 6 distinct betas for each condition).
resultsST = GLMdenoisedata(condition_split(design,6),data,stimdur,tr,'assume',[],opt, ...
                           'GLMdenoise_ST_figures');
%%

% Inspect the R2 (variance explained) map.
figure;
temp = resultsST.R2(imglookup);  % use indexing provided by imglookup to generate a surface map
temp(extrapmask) = NaN;          % some pixels are invalid; set these to NaN
imagesc(temp,[0 50]);
colormap(hot);
axis image tight ij;
%%

% This is a nearest-neighbor rendering of a flat patch of left hemisphere early visual cortex.
% Posterior to anterior is left to right, and inferior to superior is bottom to top.
% The map looks reasonably sane.

%% Fit FIR model

% The first step in performing TDM is to derive unconstrained estimates of response
% timecourses by fitting a finite impulse response (FIR) model to the data.

% Define options (no bootstrapping, no denoising (in order to speed things up))
opt2 = struct('numboots',0,'numpcstotry',0);  % to use denoising, could set numpcstotry to 20

% Use GLMdenoisedata.m to fit a finite impulse response (FIR) model that estimates
% a timecourse for every condition. We use 2 condition-splits (so we get 2 sets
% of FIR results). We estimate 31 time points (trial onset plus 30 time points);
% these time points come at a rate of 1 second since the data are prepared at 1 second.
resultsFIR = GLMdenoisedata(condition_split(design,2),data,stimdur,tr,'fir',30,opt2, ...
                            'GLMdenoise_FIR_figures');
%%

%% Inspect FIR timecourses for one vertex

% Before proceeding with TDM, let's inspect results from the FIR fitting
% to make sure things look reasonable.

% Make a plot
ix = 51344;
figure; hold on;
cmap0 = jet(6);
for p=1:2     % condition-split
  for q=1:6   % condition
    temp = squeeze(resultsFIR.modelmd(ix,:,:));  % (2 condition-splits * 6 conditions) x 31 time points
    plot(0:30,temp((q-1)*2+p,:),'-','Color',cmap0(q,:));
  end
end
straightline(0,'h','k-');
xlabel('Time after trial onset (s)');
ylabel('BOLD (%)');
%%

% In the figure, we plot two traces (reflecting the 2 condition-splits)
% for each of the 6 experimental conditions.

%% Determine an R2 threshold for the FIR results

% We want to collect timecourses from vertices that have some minimum
% level of signal-to-noise ratio for experimentally evoked BOLD responses.

% Use an automatic method for determining a reasonable threshold.
r2thresh = findtailthreshold(resultsFIR.R2(:));
%%

%% Extract FIR timecourses

% Get timecourses for vertices that pass the threshold
ix = resultsFIR.R2 > r2thresh;
timecourses = resultsFIR.modelmd(ix,:,:);  % N x (2 condition-splits * 6 conditions) x 31

% Get bias-corrected EPI intensities for these vertices
bcvalues = bc(ix);                         % N x 1

%% Perform TDM

% Now we will pass the timecourses to the core of the TDM technique;
% this core is implemented in extracthrfmanifold.m.

% Define options (we will use defaults for everything)
opt3 = struct();

% For this analysis, we have no reason to keep the 2 distinct condition-splits,
% so we will average timecourses across the splits (to improve signal-to-noise).
% We also change the dimensionality to N x time (31) x conditions (6).
temp = permute(blob(timecourses,2,2)/2,[1 3 2]);

% Use extracthrfmanifold.m to derive estimates of the latent early
% and late timecourses. The -1 indicates to also write .eps figures.
resultsTDM = extracthrfmanifold(temp,bcvalues,tr,{'TDM_figures' -1},opt3);

%% Inspect TDM figure outputs

%%
% Show manifold.png
%
% <<../TDM_figures/manifold.png>>
%%

% In the above figure, the black and gray dots indicate the
% Early and Late timecourses, respectively.
%%

%%
% Show modelfit.png
%
% <<../TDM_figures/modelfit.png>>
%%

% In the above figure, we see that 'Model Fit' is a reasonably
% good match to the data shown in 'Final'.
%%

%%
% Show timecourses.png
%
% <<../TDM_figures/timecourses.png>>
%%

% In the above figure, the Early (black) and Late (gray)
% timecourses look reasonable.

%% Fit TDM GLM

% Now that we have derived Early and Late timecourses, we can use them to
% re-fit the fMRI time-series data.

% First, we perform condition-splitting on the design matrix.
% We use 6 condition-splits (so that we get 6 betas for each condition).
design0 = condition_split(design,6);

% Next, we perform convolution with the Early and Late timecourses separately.
for p=1:length(design0)
  temp = [];
  for q=1:2
    temp = [temp conv2(full(design0{p}),resultsTDM.elhrf(q,:)')];
  end
  design0{p} = temp(1:size(design0{p},1),:);
end

% Let's take a quick look at the design matrix for the first run.
% The first set of columns correspond to the Early timecourse;
% the second set of columns correspond to the Late timecourse.
figure;
imagesc(design0{1});
colormap(gray);
%%

% Finally, we proceed to fit the GLM.

% Define options (no bootstrapping, no denoising (in order to speed things up))
opt4 = struct('numboots',0,'numpcstotry',0);  % to use denoising, could set numpcstotry to 20

% Use GLMdenoisedata.m to fit a GLM that incorporates the TDM-derived timecourses.
% Note that we already convolved in the HRFs, so we use an <hrfknobs> input of 1 which
% effectively convolves our design matrices with 1 (which doesn't change anything).
resultsTDMGLM = GLMdenoisedata(design0,data,stimdur,tr,'assume',1,opt4, ...
                              'GLMdenoise_TDMGLM_figures');

%% Inspect maps

% First, prepare a convenient betas matrix.
% The dimensionality is vertices x 3 depths x 6 condition-splits x 6 conditions x 2 loadings
betas0 = reshape(resultsTDMGLM.modelmd{2},[],3,6,6,2);

% Also prepare a convenient EPI intensity matrix.
% The dimensionality is vertices x 3 depths.
bc0 = reshape(bc,[],3);

% Look at Early and Late betas for each experimental condition.
% In the following figure, the rows indicate Early/Late,
% while the columns indicate different conditions.
depthix = 1;  % look at the first depth index (which is Depth 2)
figureprep([100 100 1200 600],1);
for q=1:2     % early/late
  for p=1:6   % condition
    subplot(2,6,(q-1)*6+p);
    temp = mean(betas0(:,depthix,:,p,q),3);
    temp = temp(imglookup);
    temp(extrapmask) = NaN;
    imagesc(temp,[-10 10]);
    colormap(cmapsign4);
    axis image tight ij off;
  end
end
%%

% Look at peak eccentricity tuning for each depth.
% In the following figure, the rows indicate Depth 2, 4, and 6,
% while the columns indicate Early/Late.
figureprep([100 100 1200 1000],1);
for q=1:2    % early/late
  for d=1:3  % depth index
    subplot(3,2,(d-1)*2+q); hold on;
    temp = centerofmass(posrect(mean(betas0(:,d,:,:,q),3)),4);
    temp = temp(imglookup);
    temp(extrapmask) = NaN;
    imagesc(temp,[1 6]);
    colormap(jet);
    axis image tight ij off;
  end
end
%%

% Look at bias-corrected EPI intensity
figureprep([100 100 1200 1000],1);
for d=1:3  % depth index
  subplot(3,1,d); hold on;
  temp = bc0(:,d);
  temp = temp(imglookup);
  temp(extrapmask) = NaN;
  imagesc(temp,[0 2]);
  colormap(gray);
  axis image tight ij off;
end
%%

%% Alternative ICA-based procedure

% We developed and tested an alternative ICA-based procedure for deriving latent timecourses.
% Although we believe the TDM method is preferable to the ICA-based procedure, we include here
% an example of performing that procedure for sake of completeness.

% Tweak the dimensionality of timecourses to N x time (31) x conditions (6) x splits (2)
tc0 = permute(reshape(timecourses,size(timecourses,1),2,size(timecourses,2)/2,[]),[1 4 3 2]);

% Define options (we will use defaults for everything)
opt5 = struct();

% Perform ICA-based procedure for deriving latent timecourses
[icahrf,icahrfmetrics,pcahrf,pcahrfmetrics,allics,chosenics] = ...
  icadecomposehrf(tc0,tr,[],{'ICA_figures' -1},opt5);
%%

% Inspect results
figure;
plot(0:30,icahrf');
straightline(0,'h','k-');
xlabel('Time (s)');
ylabel('Signal');
%%
