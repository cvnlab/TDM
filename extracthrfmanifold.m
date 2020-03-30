function results = extracthrfmanifold(data,intensity,tr,wantfigs,opt)

% function results = extracthrfmanifold(data,intensity,tr,wantfigs,opt)
%
% <data> is V x time x C with the timecourses. The first time point
%   should correspond to 0 s. V should be >= 1 and should indicate different voxels.
%   C should be >= 1 and should indicate different experimental conditions.
%   The units of <data> should be percent signal change. In general, the
%   timecourses provided in <data> should probably be filtered to omit 
%   noisy/irrelevant voxels.
% <intensity> is V x 1 with the pixel intensity. The units should be around 1,
%   where, e.g., 0.9 means 10% darker than typical.
% <tr> is the sampling rate of the timecourses in seconds (e.g. 1.0)
% <wantfigs> (optional) is
%   0 means do not create any figures
%   1 means create figures in new figure windows
%   A means write figures to folder A
%   {A mode} is like A but with specification of the <mode> input to figurewrite.m
%   Default: 1.
% <opt> (optional) is a struct with:
%   <rng> (optional) is [A B] indicating the range (in seconds) over which to
%     ensure that the mean of HRF time points is positive. Default: [0 10].
%   <numspherehistbins> (optional) is the desired number of histogram bins for quantifying
%     counts on the unit sphere. Default: 100.
%   <bins> (optional) is the bin centers for quantifying the orthographic projection
%     of the unit sphere onto the 2D image plane. Must be evenly spaced. Default: -1.1:.02:1.1.
%   <sphereparticles> (optional) is the number of particles for one hemisphere of the unit
%     sphere. We subdivide after placing particles, and we use the resulting particles to
%     perform the DC adjustment. Default: 100.
%   <hrfbasis> (optional) is the 3-dimensional orthogonal basis within which to perform 
%     the analysis. If [] or not supplied, we perform PCA on <data> and use the top 3 PCs
%     as the basis.
%   <vlengthweight> (optional) is the non-negative weight to use for vector-length
%     before averaging with the density. Both density and vector-length are regularized such
%     that they range between 0 and 1 before the weighted-averaging process.  So, for
%     example, <vlengthweight> set to 2 means use a weighted average of density and
%     vector-length with weights of 1 and 2, respectively.  Default is 1 which means
%     simply average together density and vector-length.  If you set <vlengthweight> to
%     0, this in effect means only use the density.
%   <ignorenegative> (optional) is whether to ignore timecourses that have negative loadings
%     on the first PC. If 1, we ignore. If 0, we automatically mirror across the origin
%     to ensure that loadings on the first PC are positive. Default: 0.
%   <wantpause> (optional) is whether to issue a keyboard statement right before completion.
%     Default: 0.
%
% Here is the summary of the procedure:
% - (1) Perform PCA on the <data>. Flip the sign of the 1st PC if necessary to ensure that
%       the mean of the timecourse is positive over the range <opt.rng>. The purpose of
%       PCA is to (i) provide a simple comparison HRF in the output (results.pcahrf),
%       and to (ii) establish the dimensions of the 3-dimensional space that we will
%       operate in (see next step).
% - (2) We define a 3-dimensional space given by the top 3 PCs (the first PC points
%       straight towards your eye, the second PC points to the right, and the third PC 
%       points to the top). Alternatively, you can specify the dimensions to use
%       via <opt.hrfbasis>.
% - (3) Map all timecourses in <data> onto the unit sphere. This is accomplished by projecting
%       each timecourse onto the 3 dimensions of the space. To ensure we get only "positive"
%       timecourses, we mirror the coordinates across the origin if the loading on the first
%       dimension is negative. We perform unit-length-normalization on each vector 
%       (i.e. each set of 3 coordinate loadings), taking care to also record the original vector
%       length before normalization. After these steps, all of the timecourses in <data>
%       are now points that live on the forward-facing hemisphere of the unit sphere.
% - (4) Compute a 2D histogram of the density of the timecourses with bins given by <opt.bins>.
%       Note that this is an orthogonal projection that results in increased surface area for
%       bins on the outskirts of the projection.
% - (5) In a parallel fashion to step 4, compute 2D images for vector length and <intensity>.
%       Specifically, the former image contains the median vector length across the
%       timecourses that fall in each given histogram bin, and the latter image contains
%       the median intensity across the timecourses that fall in each given histogram bin.
% - (6) Create a "culled" version of the timecourses mentioned in step 3 by removing a DC bias. 
%       This is accomplished by counting the number of timecourses associated with each 
%       sphere particle (see below), computing a histogram using <opt.numspherehistbins> bins,
%       and then finding the bin with the most counts. We then stochastically subsample the
%       timecourses to achieve the desired reduction.
% - (7) Create a 2D density image for fitting purposes. We do this by rotating the culled
%       timecourses such that the first PC of the timecourses (as determined in step 1) is
%       aligned with the positive z-axis, truncating the points to select only those that
%       are facing up, computing a 2D histogram using <opt.bins>, and then normalizing (with
%       truncation) the resulting image such that 0 maps to 0 and the maximum maps to 1.
%       In summary, compared to the results of step 4, we obtain a DC-removed, rotated,
%       truncated, and normalized 2D density image.
% - (8) Prepare the 2D vector length image for fitting. We do this by rotating and truncating
%       (in the same way as done in step 7) and then computing median for each histogram bin.
%       We then determine the mode based on a histogram of the 2D image using bins given by
%       <opt.numspherehistbins>, and then normalize (with truncation) the 2D image such that 
%       the mode maps to 0 and the maximum maps to 1.
% - (9) We prepare the final 2D image to be fit by computing a weighted average of the 
%       2D density image (step 7) and the 2D vector length image (step 8). The weighting
%       is controlled by <opt.vlengthweight>, and the resulting image is guaranteed to have
%       values between 0 and 1.
% - (10) We fit a 2D oriented Gaussian to the result of step 9. The Gaussian is parameterized
%        by two parameters for the center, two parameters for the spread along the major and
%        minor axes, one parameter for the rotation, a gain parameter, and an offset parameter.
%        The fitting is performed in a way that interprets <m> as a probability distribution;
%        that is, elements of <m> that have higher density contribute more to the error metric.
% - (11) Based on the fitted Gaussian from step 10, extract the two points that correspond to
%        the mean +/- 1 standard deviation along the principal axis of the bivariate Gaussian.
% - (12) Reconstruct the timecourses that correspond to these two points, and compute some
%        simple timecourse metrics. The peak times of the two timecourses are used to decide 
%        which timecourse is the 'Early' timecourse and which is the 'Late' timecourse.
%
% Information on the sphere particles:
% - We distribute particles on the unit sphere (using S2 Sampling Toolbox), with 
%   <opt.sphereparticles> particles for one hemisphere of the sphere, and 
%   subdivide the resulting mesh. For the default value of 100 for 
%   <opt.sphereparticles>, this process produces about 1600 particles that
%   are located on one hemisphere of the unit sphere. We use these particles
%   as quasi-histogram bin centers.
%
% Return <results> as a struct. The fields of primary importance are:
% <elhrf> is 2 x time with the final timecourses.
%   The first row is the early-peaking timecourse, and the second row is the
%   late-peaking timecourse. Each timecourse is normalized such that its
%   maximum value is 1.
% <elhrfmetrics> is 2 x peak/loc/risetime/falltime where
%   peak is the peak amplitude,
%   loc is the time to the peak in seconds, 
%   risetime is the time of rise (half of the peak) in seconds, and
%   falltime is the time of fall (half of the peak) in seconds.
%   Note that peak will be close to, but not necessarily equal to, 1.
% <pcahrf> is 1 x time with the first principal component timecourse.
%   The timecourse is normalized such that its maximum value is 1.
% <pcahrfmetrics> is 1 x peak/loc/risetime/falltime where
%   the various values have the same meaning as in <elhrfmetrics>.
% <v> is time x time with the principal component timecourses given
%   by the columns of <v>.
% <ndensity> is bins x bins with a 2D image of timecourse density,
%   produced by computing a histogram (see step 4).
% <nvectorlength> is bins x bins with a 2D image of vector length (see step 5).
% <nintensity> is bins x bins with a 2D image of intensity (see step 5).
% <xx> is bins x bins with the x-coordinate for each bin
% <yy> is bins x bins with the y-coordinate for each bin
% <nprep> is the 2D density image prepared (regularized) for fitting (step 7)
% <nprep2> is the 2D vector length image prepared (regularized) for fitting (step 8)
% <nprep3> is the weighted average of <nprep> and <prep2> (step 9).
%   This is the final 2D image that the Gaussian is fit to.
% <params> is 1 x 7 with the estimated model parameters (fitgaussian2doriented.m)
% <R2> is the coefficient of determination (R^2) indicating the percentage of 
%   variance explained by the model fit.
% <gcoord> is 3 x 2 with x- and y-coordinates for 3 points: mean, early, late.
%   Note that the z-coordinate is implicit since the coordinates reside on
%   the unit sphere.
% <fullarc> is 181 x 3 with x-, y-, and z-coordinates for 181 points equally spaced 
%   along the semicircle defined by the early and late vectors. (The spacing in-between
%   points is exactly 1 degree.) This output might be useful for choosing alternative
%   early and late timecourses. The start of the semicircle is on the "early" side;
%   the end of the semicircle is on the "late" side.
% <specialangles> 1 x 2 with the exact angles corresponding to the early and late
%   timecourses with respect to the <fullarc>. The angles lie within the range 0-180.
%
% History:
% - 2020/03/30 - version 1.1. changes include:
%     (1) tweak documentation
%     (2) add opt.ignorenegative input option
%     (3) expose results.R2 output
%     (4) new outputs results.fullarc, results.specialangles
%     (5) new example usage (exampledataset_inputsonly.mat)
% - 2019/12/07 - official check-in (version 1.0).
% - 2019/04/27 - change nprep to save the normalized version of nprep.
%
% Example usage:
% - Please see example1.m for a full example.
% - A condensed version is provided in example2.m, and is as follows:
%     load('exampledataset_inputsonly.mat','data','intensity','tr');
%     results = extracthrfmanifold(data,intensity,tr,1);

%% %%%%% SETUP

% inputs
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
if ~isfield(opt,'numspherehistbins') || isempty(opt.numspherehistbins)
  opt.numspherehistbins = 100;
end
if ~isfield(opt,'bins') || isempty(opt.bins)
  opt.bins = -1.1:.02:1.1;
end
if ~isfield(opt,'sphereparticles') || isempty(opt.sphereparticles)
  opt.sphereparticles = 100;
end
if ~isfield(opt,'hrfbasis') || isempty(opt.hrfbasis)
  opt.hrfbasis = [];
end
if ~isfield(opt,'vlengthweight') || isempty(opt.vlengthweight)
  opt.vlengthweight = 1;
end
if ~isfield(opt,'ignorenegative') || isempty(opt.ignorenegative)
  opt.ignorenegative = 0;
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
hrfix0 = ceil(opt.rng(1)/tr)+1 : floor(opt.rng(2)/tr)+1;  % indices over which to ensure mean positive
time0 = linspacefixeddiff(0,tr,numtime);

% calc
dif0 = diff(opt.bins(1:2));
opt.binsEDGES = linspacefixeddiff(opt.bins(1)-dif0/2,dif0,length(opt.bins)+1);

%% %%%%% COMPUTE SPHERE PARTICLES

% uniformly distribute particles across the surface of a unit sphere,
% subdivide the resulting mesh, and then extract the forward-facing ones.
% SP will be particles x 3 with coordinates on the unit sphere.
fprintf('calculating particles...');
[SP,Tri,~,~] = ParticleSampleSphere('N',opt.sphereparticles*2);
fv = struct('faces',Tri,'vertices',SP); 
fv = SubdivideSphericalMesh(fv,2);   % subdivide to become denser
SP = fv.vertices;
SP = SP(SP(:,3)>0,:);  % particles x 3 (these are facing forward, which we define as the 3rd dimension, or z)
if sum(SP(:,3)==0) ~= 0
  warning('some particles are exactly at z=0!');
end
fprintf('done.\n');

%% %%%%% PERFORM PCA

% pca
[u,s,v] = svd(squish(permute(data,[3 1 2]),2),0);

% flip sign of first PC if necessary
sgn = 2*((mean(v(hrfix0,1)) > 0) - 0.5);
v(:,1) = sgn*v(:,1);
u(:,1) = sgn*u(:,1);

% compute some HRF metrics
[hrf,peak,loc,risetime,falltime] = derivehrf(v(:,1)',tr,opt.rng,[],0);  % NO MEAN ADJUSTMENT

% prepare outputs while dealing with scale issues
mx = max(hrf);
pcahrf = [hrf/mx];
pcahrfmetrics = [peak/mx loc risetime falltime];

%% %%%%% COMPUTE HRF BASIS

if isempty(opt.hrfbasis)
  opt.hrfbasis = v(:,1:3);
end

%% %%%%% MAP TO UNIT SPHERE

% project onto basis
m0 = squish(permute(data,[3 1 2]),2) * opt.hrfbasis;  % N x 3 with coordinates in the HRF basis space

% if opt.ignorenegative is 0, flip each timecourse such that loading on first component is positive.
% otherwise, just discard the timecourses that load negative on the first component.
if opt.ignorenegative == 0
  m0 = bsxfun(@times,m0,sign(m0(:,1)));
else
  negix = m0(:,1) < 0;
  m0 = m0(~negix,:);
end

% map to the unit sphere
m0n = unitlength(m0,2);

% calc vector length
m0l = vectorlength(m0,2);

% explicitly define coordinates for timecourses (x+ is PC2, y+ is PC3, z+ (towards you) is PC1)
XXX0 = m0n(:,2);
YYY0 = m0n(:,3);
ZZZ0 = m0n(:,1);

%% %%%%% PREPARE INTENSITY

% prepare the intensity (i.e. the darkness).
% note that we repeat the intensity for the different conditions.
d0 = squish(permute(repmat(intensity,[1 1 numcond]),[3 1 2]),2);  % N x 1

% we have to make sure to perform the subselection if we are ignoring negative timecourses.
if opt.ignorenegative
  d0 = d0(~negix);
end

%% %%%%% COMPUTE ORTHOGRAPHICALLY MAPPED DATA

% calculate density of timecourses on unit sphere in orthographic 2D space.
% due to spherical curvature, bins on the outskirts have much more surface area.
% thus, counts in these bins will be abnormally high.
[ndensity,xx,yy] = hist2d(XXX0,YYY0,opt.bins,opt.bins);

% calculate median vector length (i.e. how wiggly the timecourses seem to be in
% terms of percent signal change) and median intensity. we use exactly the
% same orthographic bins as used for the density. note that we prepare the
% matrix such that the columns and rows proceed from negative to positive values.
nvectorlength = [];
nintensity = [];
for cc=1:length(opt.binsEDGES)-1
  ok = find(XXX0 > opt.binsEDGES(cc) & XXX0 <= opt.binsEDGES(cc+1));
  for rr=1:length(opt.binsEDGES)-1
    ok2 = find(YYY0(ok) > opt.binsEDGES(rr) & YYY0(ok) <= opt.binsEDGES(rr+1));
    nvectorlength(rr,cc) = median(m0l(ok(ok2)));
    nintensity(rr,cc) =    median(d0(ok(ok2)));
  end
end

%% %%%%% CREATE ALTERNATIVE VERSION WITH DC BIAS REMOVED

% with pure noise timecourses, this will tend to place density all over the sphere.
% we need to remove the DC bias to help the fitting of the density which will
% occur in a later step.

% perform a histogram on the unit sphere (calculate number of
% timecourses assigned to each sphere particle)
[~,iix] = max([XXX0 YYY0 ZZZ0] * SP',[],2);  % what sphere particle index achieves the max dot product?
histonsph = zeros(size(SP,1),1);  % P x 1 with counts
for p=1:size(SP,1)
  histonsph(p) = sum(iix==p);
end

% figure out the mode (i.e. the number of counts that is most common across
% the various particle positions). the assumption is that this is approximately
% the "noise" level.
[hn,hx] = hist(histonsph,opt.numspherehistbins);  
[~,iii] = max(hn);
themode = hx(iii);

% stochastically subsample the timecourses such that we remove the DC "noise" level
tokeep = logical(ones(size(XXX0,1),1));
for p=1:size(SP,1)
  originalnum = histonsph(p);  % original number of particles that lie in this bin
  desirednum = round(max(0,originalnum-themode));  % the desired number of particles for this bin
  tokeep(iix==p) = permutedim(logical([ones(desirednum,1); zeros(originalnum-desirednum,1)]));
end

% create an alternative version of the coordinates (they are culled!)
XXX0c = XXX0(tokeep);
YYY0c = YYY0(tokeep);
ZZZ0c = ZZZ0(tokeep);

%% %%%%% FIT THE DISTRIBUTION

% determine rotation to z+ (based on the first PC from the original PCA).
% note that because we will be operating on the XXX coordinates,
% we have to project the first PC onto our basis and reorder 2 3 1.
rotmatrix = xyzrotatetoz(unitlength(v(:,1)' * opt.hrfbasis(:,[2 3 1])));

% rotate the culled timecourse coordinates
newcoord = rotmatrix*[XXX0c YYY0c ZZZ0c ones(size(XXX0c))]';  % 4 x L

% truncate to select only those facing us
ok = newcoord(3,:) > 0;

% compute density of timecourses (DC-removed, rotated, truncated)
[nprep,xx,yy] = hist2d(newcoord(1,ok),newcoord(2,ok),opt.bins,opt.bins);

% we normalize such that 0 maps to 0 and the max value maps to 1.
% this does truncate values outside that range.
nprep = normalizerange(nprep,0,1,0);

% calculate median vector length in the rotated and truncated space
temp = rotmatrix*[XXX0 YYY0 ZZZ0 ones(size(XXX0))]';
ok = temp(3,:) > 0;
XXX0t = temp(1,ok);
YYY0t = temp(2,ok);
nprep2 = [];
for cc=1:length(opt.binsEDGES)-1
  ok = find(XXX0t > opt.binsEDGES(cc) & XXX0t <= opt.binsEDGES(cc+1));
  for rr=1:length(opt.binsEDGES)-1
    ok2 = find(YYY0t(ok) > opt.binsEDGES(rr) & YYY0t(ok) <= opt.binsEDGES(rr+1));
    nprep2(rr,cc) = median(m0l(ok(ok2)));
  end
end

% find the mode of the resulting vector lengths
[hn,hx] = hist(nprep2(:),opt.numspherehistbins);
[~,iii] = max(hn);
themode = hx(iii);

% prepare nprep2 by mapping the mode to 0 and the max value to 1.
% (be careful, this is tricky. is the max value robust?)
% intuitively, this will take the most common "color" and make that the zero point.
% also, at the very end, we set the NaNs to 0.
nprep2 = nanreplace(normalizerange(nprep2,0,1,themode));  

% finally, the data we will fit is the weighted-average of nprep and nprep2.
nprep3 = (1*nprep + opt.vlengthweight*nprep2) / (1+opt.vlengthweight);

% fit the bivariate Gaussian to this density image,
% interpreting the image as a probability density to weight the error metric
[params,R2] = fitgaussian2doriented(nprep3,[],[],1);

% compute the model fit
[xxALT,yyALT] = meshgrid(1:size(nprep3,2),1:size(nprep3,1));
modelfit = evalgaussian2doriented(params,xxALT,yyALT);

%% %%%%% GATHER SOME INFORMATION ON THE FIT

% gcoord will be 3 x 2. the three rows indicate the mean
% the mean plus principal axis vector, and the mean minus
% principal axis vector.

% the center of the Gaussian
gcoord = [params(1) params(2)];  % [X Y] in matrix/image space

% 1-SD along the principal axis, both positive and negative
ang = -params(7)/180*pi;
rot0 = [cos(ang) -sin(ang);
        sin(ang)  cos(ang)];
if params(3) > params(4)
  vec = [params(3) 0] * rot0;  % sx is bigger
else
  vec = [0 params(4)] * rot0;  % sy is bigger
end
gcoord(2,:) = gcoord(1,:) + vec;
gcoord(3,:) = gcoord(1,:) - vec;

% we also want contours of the 1-sd and 2-sd ellipsoids
allangs = [linspacecircular(0,2*pi,100) 2*pi];
ellipsoid = [];  % 101 x 2 (XY) x 2 (ellipsoids)
for p=1:length(allangs)
  ellipsoid(p,:,1) = gcoord(1,:) + [cos(allangs(p))*params(3)   sin(allangs(p))*params(4)  ] * rot0;
  ellipsoid(p,:,2) = gcoord(1,:) + [cos(allangs(p))*params(3)*2 sin(allangs(p))*params(4)*2] * rot0;
end

%% %%%%% BACKTRANSFORM THE FIT-RELATED QUANTITIES

% define a function such that it subtracts center of the matrix space,
% scales such that field of view of matrix space is matched to field of view
% of the original binning scheme
bfun = @(coord) (coord - (1+size(nprep3,1))/2) / size(nprep3,1) * (opt.binsEDGES(end)-opt.binsEDGES(1));

% translate the results of the Gaussian fitting to our original binning space
gcoord =    bfun(gcoord);
ellipsoid = bfun(ellipsoid);

% check whether endpoints are on the sphere
assert(all(sum(gcoord(2:3,:).^2,2) < 1),'principal axis extends beyond unit circle');

% define function to undo the rotation-to-z+
undofun = @(x) inv(rotmatrix) * [x sqrt(1-sum(x.^2,2)) ones(size(x,1),1)]';  % result is 4 x N

% keep a copy
gcoordFIT = gcoord;
ellipsoidFIT = ellipsoid;

% go from binning space to original unrotated space
gcoord =    undofun(gcoord);
ellipsoid = undofun(squish(permute(ellipsoid,[3 1 2]),2));

% finalize
gcoord = gcoord(1:2,:)';                                           % 3 [mean,plus,minus] x 2 (XY)
ellipsoid = permute(reshape(ellipsoid(1:2,:),[2 2 101]),[3 1 2]);  % 101 x 2 (XY) x 2 (ellipsoids)

%% %%%%% RECONSTRUCT TIMECOURSES AND DECIDE EARLY/LATE

% guess
hrfix1 = 2;  % for now, let's guess that the early HRF is the 2nd
hrfix2 = 3;

% reconstruct
hrf1 = [gcoord(hrfix1,:) sqrt(1-sum(gcoord(hrfix1,:).^2))] * opt.hrfbasis(:,[2 3 1])';  % 1 x time
hrf2 = [gcoord(hrfix2,:) sqrt(1-sum(gcoord(hrfix2,:).^2))] * opt.hrfbasis(:,[2 3 1])';  % 1 x time

% compute some HRF metrics
[hrf1,peak1,loc1,risetime1,falltime1] = derivehrf(hrf1,tr,opt.rng,[],0);  % NO MEAN ADJUSTMENT
[hrf2,peak2,loc2,risetime2,falltime2] = derivehrf(hrf2,tr,opt.rng,[],0);  % NO MEAN ADJUSTMENT

% prepare outputs while dealing with scale issues
mx1 = max(hrf1);
mx2 = max(hrf2);
elhrf = [hrf1/mx1;
         hrf2/mx2];
elhrfmetrics = [peak1/mx1 loc1 risetime1 falltime1;
                peak2/mx2 loc2 risetime2 falltime2];

% swap if necessary
if loc1 > loc2
  elhrf = elhrf([2 1],:);
  elhrfmetrics = elhrfmetrics([2 1],:);
  gcoord([2 3],:) = gcoord([3 2],:);  % NOTE: gcoord is now in the order of [mean,early,late]
end

%% %%%%% PREPARE FULL SEMICIRCLE FOR OUTPUT

% figure out equator points
eltemp = [gcoord(2:3,:) sqrt(1-sum(gcoord(2:3,:).^2,2))];         % 2 x 3 with coordinates of early and late
[~,~,orthbasis] = svd(eltemp,0);                                  % orthbasis(:,1:2) is 3x2 with 2 PCs
wt2 = [1 -1] * sqrt(1./(orthbasis(3,2)^2/orthbasis(3,1)^2+1));    % weight on second PC (1 x 2)
wt1 = -wt2*orthbasis(3,2)/orthbasis(3,1);                         % weight on first PC (1 x 2)
equatorpoints = orthbasis(:,1)*wt1 + orthbasis(:,2)*wt2;          % 3 x 2 with coordinates of the two equator points

% figure out the pole that points towards z+
if dot(equatorpoints(:,1),orthbasis(:,1))==1
  pole = unitlength(projectionmatrix(equatorpoints(:,1))*orthbasis(:,2));
else
  pole = unitlength(projectionmatrix(equatorpoints(:,1))*orthbasis(:,1));
end
if pole(3) < 0
  pole = -pole;  % 3 x 1 with the pole
end

% let's suppose the second column in equatorpoints is the positive x+ direction.
% then, if the x-coordinate of the early is less than the x-coordinate of the late,
% we don't have to do anything.
if eltemp(1,:)*equatorpoints(:,2) < eltemp(2,:)*equatorpoints(:,2)
else
  equatorpoints = equatorpoints(:,[2 1]);  % otherwise, we have to swap
end

% now, the first and second columns of equatorpoints correspond
% to "start" and "finish", respectively.

% construct the full arc
fullarc = constructpathonsphere([equatorpoints(:,1) pole equatorpoints(:,2)]',1);
if size(fullarc,1) < 181
  fullarc(end+1,:) = equatorpoints(:,2);  % deal with precision error
end

% compute the angles of the early and the late
specialangles = (acos(eltemp*equatorpoints(:,1))/pi*180)';  % 1 x 2

% visualization check
if 0
  figure; hold on;
  scatter3(eltemp(1,1),eltemp(1,2),eltemp(1,3),'ro');
  scatter3(eltemp(2,1),eltemp(2,2),eltemp(2,3),'bo');
  scatter3(equatorpoints(1,1),equatorpoints(2,1),equatorpoints(3,1),'cx');
  scatter3(equatorpoints(1,2),equatorpoints(2,2),equatorpoints(3,2),'cd');
  scatter3(pole(1),pole(2),pole(3),'ko');
  for zz=1:size(fullarc,1)
    scatter3(fullarc(zz,1),fullarc(zz,2),fullarc(zz,3),'k.');
  end
  axis equal;
  view(3);
end

%% %%%%% FINALLY, LET'S MAKE SOME FIGURES

% Make a figure that shows inspections of timecourses
if ~isequal(wantfigs,0)

  vlen = vectorlength(elhrf,2);
  figureprep([100 100 500 250],isequal(wantfigs,1));
  h = [];
  h(1) = plot(time0,opt.hrfbasis(:,1),'r-');
  h(2) = plot(time0,opt.hrfbasis(:,2),'g-');
  h(3) = plot(time0,opt.hrfbasis(:,3),'b-');
  h(4) = plot(time0,v(:,1),           'm-','LineWidth',2);
  h(5) = plot(time0,elhrf(1,:)/vlen(1),'k-','LineWidth',2);
  h(6) = plot(time0,elhrf(2,:)/vlen(2),'k-','LineWidth',2,'Color',[.5 .5 .5]);
  scatter(elhrfmetrics(1,2),elhrfmetrics(1,1)/vlen(1),  'ro','filled','CData',[0 0 0]);
  scatter(elhrfmetrics(1,3),elhrfmetrics(1,1)/2/vlen(1),'ro','filled','CData',[0 0 0]);
  scatter(elhrfmetrics(1,4),elhrfmetrics(1,1)/2/vlen(1),'ro','filled','CData',[0 0 0]);
  scatter(elhrfmetrics(2,2),elhrfmetrics(2,1)/vlen(2),  'ro','filled','CData',[.5 .5 .5]);
  scatter(elhrfmetrics(2,3),elhrfmetrics(2,1)/2/vlen(2),'ro','filled','CData',[.5 .5 .5]);
  scatter(elhrfmetrics(2,4),elhrfmetrics(2,1)/2/vlen(2),'ro','filled','CData',[.5 .5 .5]);
  xlabel('Time (s)');
  ylabel('Signal');
  axis([0 time0(end) -.7 .7]);
  straightline(0,'h','k-');
  h2 = straightline(0:3:time0(end),'v','k-');
  set(h2,'Color',[.7 .7 .7]);
  uistack(h2,'bottom');
  legend(h,{'Basis 1' 'Basis 2' 'Basis 3' 'PC1' 'Early' 'Late'},'Location','EastOutside');
  if ~isequal(wantfigs,1)
    figurewrite(sprintf('timecourses'),[],wantfigs{2},wantfigs{1});
  end

end

% Make a figure that inspects the model fit
if ~isequal(wantfigs,0)

  figureprep([100 100 800*2 500],isequal(wantfigs,1));
  allthings = cat(3,ndensity,nprep,nprep2,nprep3,modelfit);
  titles = {'Original' 'Regularized' 'RegularizedVL' 'Final' sprintf('Model Fit (R2=%.2f)',R2)};
  for sb=1:5
    subplot(1,5,sb); hold on;
    imagesc(xx(1,:),yy(:,1),allthings(:,:,sb));
    colormap(hot);
    cax = caxis;
    caxis([0 cax(2)]);
    axis equal tight; 
    set(straightline(0,'h','w:'),'Color',[.95 .95 .95]);
    set(straightline(0,'v','w:'),'Color',[.95 .95 .95]);
    xlabel('Loading on Basis 2');
    ylabel('Loading on Basis 3');
    drawellipse(0,0,0,1,1,[],[],'w-');
    title(titles{sb});
    if sb==1
      plot(ellipsoid(:,1,1),ellipsoid(:,2,1),'b-');
      plot(ellipsoid(:,1,2),ellipsoid(:,2,2),'c-');
    else
      plot(ellipsoidFIT(:,1,1),ellipsoidFIT(:,2,1),'b-');
      plot(ellipsoidFIT(:,1,2),ellipsoidFIT(:,2,2),'c-');
    end
  end
  if ~isequal(wantfigs,1)
    figurewrite(sprintf('modelfit'),[],wantfigs{2},wantfigs{1});
  end

end

% Make a figure that shows the final results
if ~isequal(wantfigs,0)

  % compute mixture
  ang = linspace(0,90,21);
  coordscontinuum = [];  % points x 3
  for mm=1:length(ang)
    mix0 = cos(ang(mm)/180*pi) * [gcoord(2,:) sqrt(1-sum(gcoord(2,:).^2))] + ...
           sin(ang(mm)/180*pi) * [gcoord(3,:) sqrt(1-sum(gcoord(3,:).^2))];  % 3 x 1, unit vector
    coordscontinuum(mm,:) = unitlength(mix0);
  end

  % make figure
  figureprep([100 100 800*2 300*2],isequal(wantfigs,1));
  allthings = cat(3,ndensity,nvectorlength,nintensity);
  titles = {'Density' 'Vector Length' 'EPI intensity'};
  cmaps = {hot(256) hot(256) gray(256)};
  for sb=1:3
    subplot(1,3,sb); hold on;
    switch sb
    case {1 2}
      mn0 = 0;
      mx0 = max(flatten(allthings(:,:,sb)));
    case {3}
      tmp = prctile(flatten(allthings(:,:,sb)),[25 75]);
      mn0 = 1-1.5*max(abs(1-tmp));  % for 25th and 75th, max deviation from 1, 1.5x that, and then 1 +/- 
      mx0 = 1+1.5*max(abs(1-tmp));
    end
    im0 = cmaplookup(allthings(:,:,sb),mn0,mx0,[],cmaps{sb});
    image(xx(1,:),yy(:,1),im0);
    axis equal tight; 
    set(straightline(0,'h','w:'),'Color',[.95 .95 .95]);
    set(straightline(0,'v','w:'),'Color',[.95 .95 .95]);
    xlabel('Loading on Basis 2');
    ylabel('Loading on Basis 3');
    drawellipse(0,0,0,1,1,[],[],'w-');
    title(sprintf('%s (%.2f to %.2f)',titles{sb},mn0,mx0));
      % ellipsoids
    plot(ellipsoid(:,1,1),ellipsoid(:,2,1),'b-');
    plot(ellipsoid(:,1,2),ellipsoid(:,2,2),'c-');
      % continuum
    plot(coordscontinuum(:,1),coordscontinuum(:,2),'k-');
      % first PC
    pc0 = v(:,1)' * opt.hrfbasis(:,[2 3]);
    scatter(pc0(1),pc0(2),'mo','filled');
      % early and late
    scatter(gcoord(2,1),gcoord(2,2),'ko','filled');
    set(scatter(gcoord(3,1),gcoord(3,2),'ko','filled'),'CData',[.5 .5 .5]);
    %  % plot every 0.5 sd from the full arc
    %halfsd = (specialangles(2)-mean(specialangles))/2;
    %angs = [fliplr(mean(specialangles):-halfsd:0) mean(specialangles)+halfsd:halfsd:180];
    %scatter(fullarc(round(angs)+1,1),fullarc(round(angs)+1,2),'r.');

  end
  if ~isequal(wantfigs,1)
    figurewrite(sprintf('manifold'),[],wantfigs{2},wantfigs{1});
  end

end

%% %%%%% PREPARE OUTPUTS

allvars = {'elhrf' 'elhrfmetrics' 'pcahrf' 'pcahrfmetrics' 'gcoord' 'gcoordFIT' 'ellipsoid' 'ellipsoidFIT' ...
 'ndensity' 'nvectorlength' 'nintensity' 'nprep' 'nprep2' 'nprep3' 'modelfit' 'params' 'R2' ...
 'v' 'xx' 'yy' 'xxALT' 'yyALT' 'fullarc' 'specialangles' 'opt'};
clear results;
for p=1:length(allvars)
  results.(allvars{p}) = eval(allvars{p});
end

%% %%%%% PAUSE?

if opt.wantpause
  keyboard;
end
