function [outlierthresh, numoutlierperiods] = makefiltereeg_Intan(anprocesseddatadir, index, freq, eegsamprate, outliernstd)
%
% INPUTS
% freq -- frequencies for bandpass filtering, usually [1 300] for typical
%   eeg
% eegsamprate -- final sampling rate of eeg, usually 2000
% Updated ALP 4/18/19 to include outlier periods in saving

%load file
load([anprocesseddatadir, 'raweeg', num2str(index(3))])
samprate = raweeg{index(1)}{index(2)}{index(3)}.samprate;

%replace outliers in LFP 
[neweeg, outlierindices, outlierthresh, outlierperiods, numoutlierperiods, ~] = interpoveroutliers_041819(raweeg{index(1)}{index(2)}{index(3)}.data, outliernstd, 2000); %changed from 2 (for mV) to 2000 (for uV)

%filter eeg
eeg = raweeg;
eeg{index(1)}{index(2)}{index(3)}.outlierindices = outlierindices;
eeg{index(1)}{index(2)}{index(3)}.outlierthresh = outlierthresh;
eeg{index(1)}{index(2)}{index(3)}.outlierperiods = outlierperiods;
eeg{index(1)}{index(2)}{index(3)}.numoutlierperiods = numoutlierperiods;


eeg = filtereeg(eeg, neweeg, index,  anprocesseddatadir, [freq], eegsamprate);


