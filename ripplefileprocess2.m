function ripplefileprocess2(directoryname,  eegfilename, filterdir, varargin)
%RIPPLEDAYPROCESS(directoryname,fileprefix,days, options)
%
%Applies a ripple filter to all epochs for each day and saves the data in
%in the EEG subdirectory of the directoryname folder. note data will be
%downsampled to make sampling rate 2000 before filtering to match filter sampling
% For files with 3 Number index instead of 4.  
%
%directoryname - example '/data99/user/animaldatafolder', a folder
%                containing processed matlab data for the animal
%
%index - eegfile to be filtered, eg 'eeg31.mat.' if empty,
%           filter all the *eeg* files in directoryname
%
%options -
%
%		'f', matfilename
%			specifies the name of the mat file containing the
%			ripple filter to use
%			(default /usr/local/filtering/ripplefilter.mat).
%			Note that the filter must be called 'ripplefilter'.
%

%
% AS updated 2/14/13

f = '';
defaultfilter = [filterdir, 'ripplefilter.mat'];


%system = 2;

%set variable options
for option = 1:2:length(varargin)-1
    switch varargin{option}
        case 'f'
            f = varargin{option+1};
    end
end

minint = -32768;

% if the filter was not specified, load the default
if isempty(f)
    load(defaultfilter); 
else
    eval(['load ', [filterdir, f] ]);
end

% create the list of files for this day that we should filter
if ~isempty(eegfilename)
    flist{1} = sprintf('%s%s', directoryname, eegfilename);
else
    tmpflist = dir(sprintf('%s*eeg*.mat', directoryname));
    flist = cell(size(tmpflist));
    for i = 1:length(tmpflist)
        flist{i} = sprintf('%s%s', directoryname, tmpflist(i).name);
    end
end

% go through each file in flist and filter it
for fnum = 1:length(flist)
    
    %load the eeg file
    load(flist{fnum});
    index = eeg{end}{end}{end}.index;
    
    
    %make samprate 2000
    if eeg{index(1)}{index(2)}{index(3)}.samprate ~= 2000
        down = eeg{index(1)}{index(2)}{index(3)}.samprate/2000;
        if rem(down, 2)==0
            eeg{index(1)}{index(2)}{index(3)}.data = eeg{index(1)}{index(2)}{index(3)}.data(1:down:end, :);
        else
            error('sampling rate is not divisible by 2000')
        end
    end
    
    if ~isempty(eeg{index(1)}{index(2)}{index(3)}.data)
        % filter it and save the result
        ripple{index(1)}{index(2)}{index(3)} = filtereeg2_Intan(eeg{index(1)}{index(2)}{index(3)}, ...
            ripplefilter, 'int16', 0);
        clear eegrec
        %add index including channel
        ripple{index(1)}{index(2)}{index(3)}.index = index;
        
        % save the resulting file
%         ripplefile = sprintf('%sEEG/ripple%02d-%04d-%02d.mat', ...
%             directoryname,  index(1), index(2), index(3));
        ripplefile = [directoryname, '/EEG/', 'ripple', num2str(index(3))]; %ALP 6/5/18
        save(ripplefile, 'ripple');
        clear ripple
    else
        warning(['could not filter index ', num2str(index), ' cuz no data'])
    end
end

