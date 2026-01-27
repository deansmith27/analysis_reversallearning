function ripple = ripplefileprocess_DC(eeg, filter)
%RIPPLEDAYPROCESS(directoryname,fileprefix,days, options)
%based on ripplefileprocess2, updated for compatibility with full pipeline
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

defaultfilter = filter;

% go through each file in flist and filter it
for fnum = 1:size(eeg,1)%each lfp data
    
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
        ripple.index = index;
        
    else
        warning(['could not filter index ', num2str(index), ' cuz no data'])
    end
end

