function [icahrf,icahrfmetrics,pcahrf,pcahrfmetrics,allics,chosenics] = icadecomposehrf(data,tr,fasticaopts,wantfigs,opt)

% function [icahrf,icahrfmetrics,pcahrf,pcahrfmetrics,allics,chosenics] = icadecomposehrf(data,tr,fasticaopts,wantfigs,opt)
%
% <data> is V x time x C x test/retest with the timecourses. The first time point
%   should correspond to 0 s. V should be >= 1 and should indicate different voxels.
%   C should be >= 1 and should indicate different experimental conditions.
%   test/retest has a dimensionality of 2 and should contain two independent
%   measurements of the timecourses. The units of <data> should be percent signal
%   change. In general, the timecourses provided in <data> should probably be filtered 
%   to omit noisy/irrelevant voxels.
% <tr> is the sampling rate of the timecourses in seconds (e.g. 1.0)
% <fasticaopts> (optional) is a cell vector of options to pass to fastica.m.
%   Default: {'stabilization' 'on'}.
% <wantfigs> (optional) is
%   0 means do not create any figures
%   1 means create figures in new figure windows
%   A means write figures to folder A
%   {A mode} is like A but with specification of the <mode> input to printnice.m
%   Default: 1.
% <opt> (optional) is a struct with:
%   <rng> (optional) is [A B] indicating the range (in seconds) over which to
%     ensure that the mean of HRF time points is positive. Default: [0 10].
%   <consistencythresh> (optional) is the R^2 percentage threshold such that if exceeded,
%     two timecourses are considered to be consistent. Default: 50.
%   <maxhrftouse> (optional) is the maximum number of timecourses to consider
%     for quantification of variance explained. The purpose of this is to reduce
%     execution time. Timecourses are randomly chosen. Default: 50000.
%   <fracconsider> (optional) is the fraction of the top-variance-explaining
%     ICs to consider. Default: 1/2.
%   <wantpause> (optional) is whether to issue a keyboard statement right before completion.
%     Default: 0.
%
% Basic outline of the procedure:
% - For each test/retest case, call fastica.m to obtain components (number is equal
%   to the number of time points in <data>).
% - Greedily order the components (within each case) such that they maximize the variance
%   explained in the data. (I.e., choose the one component that maximizes variance,
%   choose the next component such that when combined with the first component variance
%   is maximized, etc.)
% - For each component, ensure it is positive over the range given in <opt.rng> (by sign
%   flipping), and then normalize each component to unit length.
% - Calculate the amount of variance that each component in the test set can explain in
%   each component in the retest set.
% - Iterate through each component in the test set, greedily assigning the best match
%   in the retest set to that component. This creates a new order for the components
%   in the retest set.
% - From the components that exist for the test set, choose components that satisfy both:
%   (1) high consistency (<opt.consistencythresh>) between that component and its match
%       in the retest set
%   (2) within the top <opt.fracconsider> of components.
% - Of those components that remain, determine the pair of components that explain
%   as much variance in the data as possible.
% - Finally, order the components in the determined pair such that the earlier-peaking
%   component is first and the later-peaking component is second.
%
% History:
% - 2020/06/20 - minor optimizations
%
% Return:
% <icahrf> is 2 x time with the ICA-derived timecourses.
%   The first row is the early-peaking timecourse, and the second row is the
%   late-peaking timecourse. Each timecourse is normalized such that its
%   maximum value is 1.
% <icahrfmetrics> is 2 x peak/loc/risetime/falltime where
%   peak is the peak amplitude,
%   loc is the time to the peak in seconds, 
%   risetime is the time of rise (half of the peak) in seconds, and
%   falltime is the time of fall (half of the peak) in seconds.
%   Note that peak will be close to, but not necessarily equal to, 1.
% <pcahrf> is 1 x time with the first principal component timecourse.
%   The timecourse is normalized such that its maximum value is 1.
% <pcahrfmetrics> is 1 x peak/loc/risetime/falltime where
%   the various values have the same meaning as in <icahrfmetrics>.
% <allics> is time x time x 3. The first slice contains the
%   first-split ICs (time x orderedICs), and the second slice
%   contains the second-split ICs (time x orderedICs) which have
%   been matched to the ICs in the first split. Each IC timecourse
%   is sign-flipped (according to opt.rng) and is unit length.
%   The third slice contains the average of the first two slices
%   (and is what is used to generate <icahrf>).
% <chosenics> is 1 x 2 with the indices of the ICs deemed
%   to be the early HRF (first entry) and the late HRF (second entry).
%   These indices refer to the second dimension of <allics>.

%%%%% SETUP

% inputs
if ~exist('fasticaopts','var') || isempty(fasticaopts)
  fasticaopts = {'stabilization' 'on'};
end
if ~exist('wantfigs','var') || isempty(wantfigs)
  wantfigs = 1;
end
if ~exist('opt','var') || isempty(opt)
  opt = struct;
end

% opt inputs
if ~isfield(opt,'rng') || isempty(opt.rng)
  opt.rng = [0 10];
end
if ~isfield(opt,'consistencythresh') || isempty(opt.consistencythresh)
  opt.consistencythresh = 50;
end
if ~isfield(opt,'maxhrftouse') || isempty(opt.maxhrftouse)
  opt.maxhrftouse = 50000;
end
if ~isfield(opt,'fracconsider') || isempty(opt.fracconsider)
  opt.fracconsider = 1/2;
end
if ~isfield(opt,'wantpause') || isempty(opt.wantpause)
  opt.wantpause = 0;
end

% some more input stuff
if ischar(wantfigs)
  wantfigs = {wantfigs []};
end

% check that all the data are valid
assert(all(isfinite(data(:))));

% calc
numvox =  size(data,1);
numtime = size(data,2);
numcond = size(data,3);
numics = numtime;
hrfix0 = ceil(opt.rng(1)/tr)+1 : floor(opt.rng(2)/tr)+1;  % indices over which to ensure mean positive
time0 = linspacefixeddiff(0,tr,numtime);

% internal constants
numcondsplit = 2;  % number of "condition splits" is fixed at 2

%%%%% PERFORM PCA

% perform PCA
[pcahrf,peak,loc,risetime,falltime] = derivehrf(squish(permute(data,[4 3 1 2]),3),tr,opt.rng);

% deal with scale issues
mx = max(pcahrf);
pcahrf = pcahrf / mx;
pcahrfmetrics = [peak/mx loc risetime falltime];

%%%%% PERFORM ICA

% initialize
xx = zeros(numtime,numcond*numvox,numcondsplit,'single'); % time x C*V x test/retest. this is the data, now maximally independent. A*xx gives the original data.
A = zeros(numtime,numtime,numcondsplit); % time x time x test/retest. each column is a timecourse.
W = zeros(numtime,numtime,numcondsplit); % time x time x test/retest. each row is a vector on which the data are projected.
listsofar = zeros(numcondsplit,numtime);   % test/retest x time with the index of the IC chosen at each iteration.
finalvarexp = zeros(numcondsplit,numtime); % test/retest x time with variance explained thus far at each tieration

% process test and retest
for ccc=1:numcondsplit

  % prepare the data (note that we smash conditions together)
  thedata = squish(permute(data(:,:,:,ccc),[3 1 2]),2)';  % time x C*V

  % call fastica
  [xx(:,:,ccc),A(:,:,ccc),W(:,:,ccc)] = fastica(thedata,'numOfIC',numtime,fasticaopts{:});

  % compute variance explained!
  listsofar0 = [];                  % index of which one chosen at each iteration
  finalvarexp0 = zeros(1,numtime);  % list of variance explained at each iteration
  for q=1:numtime
    varx = NaN(1,numtime);
    for p=1:numtime
      if ismember(p,listsofar0)
        continue;
      end
      h = olsmatrix(A(:,[listsofar0 p],ccc))*thedata;  % 1 x C*V with the estimated weights
      varx(p) = calccod(flatten(A(:,[listsofar0 p],ccc)*h),flatten(thedata),[],[],0);  % compute variance explained
    end
    [mm,ii] = max(varx);
    finalvarexp0(q) = mm;
    listsofar0 = [listsofar0 ii];
  end
  
  % record
  listsofar(ccc,:) = listsofar0;
  finalvarexp(ccc,:) = finalvarexp0;

end

% clean up
clear thedata h;

%%%%% MASSAGE ICA RESULTS 

% get ICs (with greedy variance explained approach)
IC1 = A(:,listsofar(1,:),1);  % time x orderedICs
IC2 = A(:,listsofar(2,:),2);  % time x orderedICs

% flip ICs if necessary and unit-length-normalize
for p=1:numics
  sgn = 2*((mean(IC1(hrfix0,p))>0)-.5);  % 1 or -1
  IC1(:,p) = unitlength(sgn * IC1(:,p));
  sgn = 2*((mean(IC2(hrfix0,p))>0)-.5);  % 1 or -1
  IC2(:,p) = unitlength(sgn * IC2(:,p));
end

% compute cross correlation using 0-relative R2 regression (timecourses IC1 x timecourses IC2)
cmatrix = calcconfusionmatrix(IC2,IC1,6);

% iteratively determine the best match that links IC2 to IC1 [the order in IC1 matters! (the first several have privilege)]
matchesinsecond = zeros(1,size(cmatrix,1));
for p=1:size(cmatrix,1)
  remainder = setdiff(1:size(cmatrix,2),matchesinsecond);
  [~,ix] = max(cmatrix(p,remainder));
  matchesinsecond(p) = remainder(ix);
end

% reorder IC2 in the order that best matches IC1
IC2 = IC2(:,matchesinsecond);

% compute a version of the ICs that averages across the test/retest
ICs = (IC1+IC2)/2;

% let's re-compute cmatrix to deal with the ordering issue
cmatrix = calcconfusionmatrix(IC2,IC1,6);

% extract the 0-relative R2 metric that quantifies goodness of match between IC1 and IC2
consistencyR2 = diag(cmatrix);

% let's be heavy-handed and do some reduction.
% which ICs are both (1) consistent and (2) within the top fraction of variance explaining?
goodix = flatten(intersect(find(consistencyR2 >= opt.consistencythresh),1:round(numics*opt.fracconsider)));

% we want to find the best pair of ICs that explain as much of the data as possible.
% we choose from only those ICs that pass our minimal requirements.

% first let's prepare the data
thedata = single(squish(permute(data,[3 1 2 4]),2));  % C*V x time x test/retest
ix = picksubset(1:size(thedata,1),round(opt.maxhrftouse/2));
thedata = squish(permute(thedata(ix,:,:),[3 1 2]),2)';  % time x stuff

% now let's try all possible pairs
varvar = NaN(numics,numics);
for rr=goodix
  for cc=goodix
    if rr>=cc
      continue;
    end
    X = ICs(:,[rr cc]);
    thefit = (X*olsmatrix(X))*thedata;  % time x stuff
    varvar(rr,cc) = calccod(flatten(thefit),flatten(thedata),[],[],0);
  end
end
[~,ix] = max(varvar(:));
[mxROW,mxCOL] = ind2sub([numics numics],ix);
hrfix1 = mxROW;  % for now, let's guess that the early HRF is the row one
hrfix2 = mxCOL;

% compute some HRF metrics
[hrf1,peak1,loc1,risetime1,falltime1] = derivehrf(ICs(:,hrfix1)',tr,opt.rng,[],0);  % NO MEAN ADJUSTMENT
[hrf2,peak2,loc2,risetime2,falltime2] = derivehrf(ICs(:,hrfix2)',tr,opt.rng,[],0);  % NO MEAN ADJUSTMENT

% prepare outputs while dealing with scale issues
mx1 = max(hrf1);
mx2 = max(hrf2);
icahrf = [hrf1/mx1;
          hrf2/mx2];
icahrfmetrics = [peak1/mx1 loc1 risetime1 falltime1;
                 peak2/mx2 loc2 risetime2 falltime2];

% swap if necessary
if loc1 > loc2
  icahrf = icahrf([2 1],:);
  icahrfmetrics = icahrfmetrics([2 1],:);
  [hrfix1,hrfix2] = swap(hrfix1,hrfix2);
end

% clean up
clear thedata;

%%%%% PREPARE OUTPUTS

allics = cat(3,IC1,IC2,ICs);
chosenics = [hrfix1 hrfix2];

%%%%% HRF INSPECTIONS

if ~isequal(wantfigs,0)

  % prepare
  robustrng0 = robustrange(picksubset(flatten(data(:,:,:,1)),10000));
  yyy = squish(permute(data,[4 3 1 2]),3);  % stuff x time
  ispos0 = mean(yyy(:,hrfix0),2) > 0;
  
  % loop over positive and negative
  todo = {ispos0 ~ispos0};
  labels = {'POS' 'NEG'};
  for zz=1:length(todo)
    yyy0 = yyy(todo{zz},:);
    xxx = repmat(time0,[size(yyy0,1) 1]);
    figureprep([100 100 500 500],isequal(wantfigs,1));
    hold on;
    binx = time0;
    biny = linspace(robustrng0(1),robustrng0(2),100);
    [n,x,y] = hist2d(xxx,yyy0,binx,biny);
    imagesc(x(1,:),y(:,1),log(n));
    colormap(jet);
    axis tight;
    ax = axis;
    plot(time0,ax(4)/2*pcahrf,'c-','LineWidth',2);
    plot(time0,ax(4)/2*icahrf(1,:),'w-','LineWidth',2);
    plot(time0,ax(4)/2*icahrf(2,:),'g-','LineWidth',2);
    set(straightline(0,'h','y-'),'LineWidth',2);
    xlabel('Time (s)');
    ylabel('Response');
    if ~isequal(wantfigs,1)
      figurewrite(sprintf('hrfinspection%s',labels{zz}),[],wantfigs{2},wantfigs{1});
    end
  end

end

%%%%% VISUALIZE ICA RESULTS

if ~isequal(wantfigs,0)

  % pre-compute some error bar stuff
  mxs = max(allics(:,chosenics,3),[],1);  % 1 x 2
  tcs = allics(:,chosenics,1:2);          % time x 2 x 2
  tcs = bsxfun(@rdivide,tcs,mxs);         % divide by normalization factor
  se0 = std(tcs,[],3)/sqrt(2);            % time x 2

  % visualize
  figureprep([100 100 600 1100],isequal(wantfigs,1));

  subplot(4,2,1); hold on;
  bar(consistencyR2,1);
  axis([0 numics+1 0 100]);
  straightline(opt.consistencythresh,'h','r-');
  set(straightline(goodix,'v','b-'),'LineWidth',3);
  set(straightline(chosenics,'v','c-'),'LineWidth',2);
  title(sprintf('Consistency (found %d)',length(goodix)));
  ylabel('Variance (%)');

  subplot(4,2,2); hold on;
  bar(finalvarexp(1,:),1);
  axis([0 numics+1 0 100]);
  straightline(round(numics*opt.fracconsider)+0.5,'v','r-');
  set(straightline(goodix,'v','b-'),'LineWidth',3);
  set(straightline(chosenics,'v','c-'),'LineWidth',2);
  title('Greedy variance ordering from first split');
  ylabel('Variance (%)');

  subplot(4,2,3); hold on;
  imagesc(cmatrix,[0 100]); axis image tight;
  colormap(jet); colorbar;
  title('Matching ICs across splits');
  ylabel('ICs from split 1');
  xlabel('ICs from split 2');

  subplot(4,2,4); hold on;
  imagesc(varvar,[0 100]); axis image tight;
  colormap(jet); colorbar;
  scatter(mxCOL,mxROW,'ko');
  title('Pair-wise variance evaluation');
  ylabel('IC number');
  xlabel('IC number');

  subplot(4,1,3); hold on;
  cmap0 = jet(length(goodix));
  for p=1:length(goodix)
    plot(time0,IC1(:,goodix(p)),'r-','Color',cmap0(p,:),'LineWidth',2);
    plot(time0,IC2(:,goodix(p)),'r-','Color',cmap0(p,:),'LineWidth',2);
  end
  axis([0 time0(end) -.7 .7]);
  straightline(0,'h','k-');
  straightline(0:3:time0(end),'v','k:');
  xlabel('Time (s)');
  ylabel('Response');

  subplot(4,1,4); hold on;
  errorbar3(time0,icahrf(1,:),se0(:,1)','v',([0 0 0]+2*[1 1 1])/3);
  errorbar3(time0,icahrf(2,:),se0(:,2)','v',([.5 .5 .5]+2*[1 1 1])/3);
  plot(time0,icahrf(1,:),'r-','Color',[0 0 0],'LineWidth',3);
  plot(time0,icahrf(2,:),'r-','Color',[.5 .5 .5],'LineWidth',3);
  axis([0 time0(end) -.7 1.2]);
  straightline(0,'h','k-');
  straightline(0:3:time0(end),'v','k:');
  scatter(icahrfmetrics(1,2),icahrfmetrics(1,1),  'ro','filled','CData',[0 0 0]);
  scatter(icahrfmetrics(1,3),icahrfmetrics(1,1)/2,'ro','filled','CData',[0 0 0]);
  scatter(icahrfmetrics(1,4),icahrfmetrics(1,1)/2,'ro','filled','CData',[0 0 0]);
  scatter(icahrfmetrics(2,2),icahrfmetrics(2,1),  'ro','filled','CData',[.5 .5 .5]);
  scatter(icahrfmetrics(2,3),icahrfmetrics(2,1)/2,'ro','filled','CData',[.5 .5 .5]);
  scatter(icahrfmetrics(2,4),icahrfmetrics(2,1)/2,'ro','filled','CData',[.5 .5 .5]);
  xlabel('Time (s)');
  ylabel('Response');

  if ~isequal(wantfigs,1)
    figurewrite(sprintf('icainspection'),[],wantfigs{2},wantfigs{1});
  end

end

%%%%% PAUSE?

if opt.wantpause
  keyboard;
end
