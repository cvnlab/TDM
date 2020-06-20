function [hrf,peak,loc,risetime,falltime] = derivehrf(data,tr,rng,wantmetrics,wantmeanadj)

% function [hrf,peak,loc,risetime,falltime] = derivehrf(data,tr,rng,wantmetrics,wantmeanadj)
%
% <data> is V x time with one or more timecourses.
%   The first time point should correspond to 0 s.
% <tr> is the sampling rate of the timecourses in seconds (e.g. 1.0)
% <rng> (optional) is [A B] indicating the range (in seconds) over which to
%   ensure that the mean of HRF time points is positive. Default: [0 10].
% <wantmetrics> (optional) is
%   0 means just calculate the HRF via PCA. (The HRF is returned as unit length.)
%   1 means to also scale the HRF to the mean absolute beta
%   2 means to also calculate the HRF peak and loc, as well as rise and fall times
%   Default: 2.
% <wantmeanadj> (optional) is whether to use zero-mean strategy. Default: 1.
%
% Use PCA to extract an HRF common to the timecourses in <data>.
% We ensure the returned HRF is mean-positive over <rng>.
%
% If <wantmeanadj>, we use a zero-mean strategy (set the mean of each timecourse
%   to zero before performing PCA) to deal with weird DC effects, and we 
%   adjust the derived HRF to start from 0 (by subtracting the first time point).
%
% We use sinc interpolation to find the peak, and we use simple linear
% interpolation to calculate rise and fall times.
%
% The PCA strategy nicely handles diversity of timecourses in <data>,
% handles negative timecourses well, and is sensitive to units (and thus
% will tend to pay more attention to timecourses with large responses).
%
% Outputs are returned as NaN if not requested to be computed.
% 
% Return:
% <hrf> is 1 x time
% <peak> is the peak amplitude
% <loc> is the time to the peak in seconds
% <risetime> is the time of rise (half of the peak) in seconds
% <falltime> is the time of fall (half of the peak) in seconds
%
% Example:
% hrfs = repmat(getcanonicalhrf(4,1),[100 1]);
% hrfs = hrfs + randn(size(hrfs))*0.3;
% [hrf,peak,loc,risetime,falltime] = derivehrf(hrfs,1);
% figure; hold on;
% plot(0:size(hrfs,2)-1,hrfs,'r-');
% plot(0:size(hrfs,2)-1,hrf,'k-','LineWidth',2);
% scatter(loc,peak,'go','filled');
% scatter(risetime,peak/2,'bo','filled');
% scatter(falltime,peak/2,'bo','filled');
% straightline(0,'h','b-');
% xlabel('Time (s)');

% internal constants
fctr = 100;  % factor by which to upsample to find peak

% inputs
if ~exist('rng','var') || isempty(rng)
  rng = [0 10];
end
if ~exist('wantmetrics','var') || isempty(wantmetrics)
  wantmetrics = 2;
end
if ~exist('wantmeanadj','var') || isempty(wantmeanadj)
  wantmeanadj = 1;
end

% fill in dummy values (because in some cases we do not compute them)
peak = NaN;
loc = NaN;
risetime = NaN;
falltime = NaN;

% perform PCA
if wantmeanadj

  % do PCA on the zero-mean version
  fun = @(x) x'*x;
  [~,~,v] = svd(fun(zeromean(data,2)),0);

  % get the first PC and adjust it to be relative to the value at time 0
  hrf = v(:,1)' - v(1,1);         % 1 x time

else

  % do PCA
  [~,~,v] = svd(data'*data,0);

  % get the first PC
  hrf = v(:,1)';                  % 1 x time

end

% flip the HRF if necessary
ix = ceil(rng(1)/tr)+1 : floor(rng(2)/tr)+1;
sgn = 2*((mean(hrf(ix))>0)-.5);     % 1 or -1
hrf = sgn * hrf;                    % apply flip to HRF

% do the scaling?
if wantmetrics >= 1

  % fit the HRF to the original data
  h = olsmatrix(hrf')*data';       % 1 x N [this is the beta weight on the unit-length HRF]

  % scaling is just the mean abs beta
  sc = mean(abs(h));              % determine the scaling that we should apply
  hrf = sc*hrf;                   % go ahead and scale

  % do various HRF metrics?
  if wantmetrics >= 2

    % calculate the peak
    [peak,loc] = calcpeak(hrf,2,'lanczos3',fctr,1);

    % calc the rise time
    ix = 1:floor(loc);
    temp = find(hrf(ix)>peak/2);
    if isempty(temp) || firstel(temp)==1
      risetime = NaN;
    else
      ix2 = [firstel(temp)-1 firstel(temp)];
      risetime = interp1(hrf(ix2),ix2,peak/2,'linear');
    end

    % calc the fall time
    ix = ceil(loc):length(hrf);
    temp = find(hrf(ix)<peak/2);
    if isempty(temp) || firstel(temp)==1
      falltime = NaN;
    else
      ix2 = (ceil(loc)-1) + [firstel(temp)-1 firstel(temp)];
      falltime = interp1(hrf(ix2),ix2,peak/2,'linear');
    end

    % convert to seconds
    loc = (loc-1)*tr;
    risetime = (risetime-1)*tr;
    falltime = (falltime-1)*tr;
  
  end

end
