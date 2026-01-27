function [newdata, outlierindices, outlierthresh, outperiods, numoutlierperiods,...
    outliers] = interpoveroutliers_041819(data,nstd, peakval, varargin)
% replaces outliers nstd standard deviations above mean with linear
% interpolation over data before and after the outliers.. excludes data
% above peakval for outliers estimation
%
% INPUTS
%   data -- vector (Nx1 or 1xN) with data to remove outliers
%   nstd -- number of standard deviations to exclude outliers, 1x1 or 1x2,
%       1x2 means different standard deviation for data below mean (first #) and data
%       above mean (2nd #).  since spikes go down, may have more permissive criteria
%       (higher nstd) for data below mean
%   peakval -- exclude above peakval and below -peakval top computer mean
%       and std, suggest 2mV
%
% OPTIONS
%   'samplesaroundoutlier' -- samples to cut out around outlier in number of
%       data points (eg 2 samples to cut out 1ms if sampling rate is 2000)
%
% OUTPUTS
%   newdata -- data with outliers removed
%   outlierindices -- indices that were replaced
%   outlierthresh --threshold for outliers, eg mean +/- nstd*stdev
%
% ASinger 8/21/15
% Updated SP 6/19/18 to include more outlier indices
% Updated ALP 4/18/19 to include outlier periods in output and saving

%set options
samplesaroundoutlier = 0;
for option = 1:2:length(varargin)-1
    switch varargin{option}
        case 'samplesaroundoutlier'
            samplesaroundoutlier = varargin{option+1};
    end
end

tempout = data;

%get outlier
mn = mean( data(data<peakval & data>-peakval) );
stdev = std(data (data<peakval & data>-peakval));
if length(nstd) == 1
outliers = data > mn + nstd*stdev | data < mn - nstd * stdev; %outliers = 1 when data is above or below thresh
elseif length(nstd) == 2
outliers = data > mn + nstd(2)*stdev | data < mn - nstd(1) * stdev; %outliers = 1 when data is above or below thresh
end    

%replace outliers
if any(outliers)
    runs = contiguous(outliers, 1); %get runs of 1s
    outperiods = runs{1,2}; %first column is start of outlier period, second column is end
    if samplesaroundoutlier~=0
        tempoutperiods = [outperiods(:,1)-samplesaroundoutlier outperiods(:,2)+samplesaroundoutlier];
        tempoutperiods(outperiods<1) = 1;
        tempoutperiods(outperiods>size(data,2)) = size(data,2);
        outperiods = combineExcludePeriods(tempoutperiods,tempoutperiods);
        %outperiods = tempoutperiods;
    end
    
    for p = 1:size(outperiods,1) %changed from 1:size(outperiods,1)
        if outperiods(p,1) == 1 %added ALP 2/14/2020 for cases when the first outlier starts at 1
            lindata = linspace(data(outperiods(p,1)), data(outperiods(p,2)+1), outperiods(p,2)-outperiods(p,1)+1);
        else
            lindata = linspace(data(outperiods(p,1)-1), data(outperiods(p,2)+1), outperiods(p,2)-outperiods(p,1)+1 );
        end
        tempout(outperiods(p,1):outperiods(p,2)) = lindata;
    end
    numoutlierperiods = size(outperiods,1);
else
    numoutlierperiods = 0;
    outperiods = [];
end

%make outputs
newdata = tempout;
% outlierindices = find(outliers); %changed SP 6.19.18 to below so that
% samples around outliers would also be included in ripple exclusion code
outlierindices = [];
for i = 1:numoutlierperiods
    outlierindices = [outlierindices [outperiods(i,1): outperiods(i,2)]];
end

if length(nstd) == 1
    outlierthresh = [mn - nstd*stdev mn + nstd*stdev ];
elseif length(nstd) == 2
    outlierthresh = [mn - nstd(1)*stdev mn + nstd(2)*stdev ];
end


