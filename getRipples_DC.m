function [rawDataBySessionNeural] = getRipples_DC(dirs, params, saveNeuralPath, plotRawTraceForRipples, plotRipples)
%adapted from filtereeg2_Intan.m, extractripples3.m,
% ripplepostfileprocess2, findpowerratioripplevsabove2
%last checked JLK 1/8/26
%DC changing to fully incorporate ripplefileprocess, extractripples3, 
%ripplepostfileprocess2, getBestRippleChan_simple, plotLFPperiods

%% extract ripples across session %%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% load session data %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load([saveNeuralPath '\rawDataBySessionNeural.mat'])
load([saveNeuralPath '\sessionPyrLayerInfo.mat'])
%tmp code
if isfield(rawDataBySessionNeural,'ripplesBad')
    rawDataBySessionNeural = rmfield(rawDataBySessionNeural,'ripplesBad');
end
if isfield(rawDataBySessionNeural,'ripplesGood')
    rawDataBySessionNeural = rmfield(rawDataBySessionNeural,'ripplesGood');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% initialize params and load filters %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
theta = [];
beta = [];
delta = [];
tdbratio = [];
ripple = [];
ripples = [];
ripplesBad = [];
chCtr = 0;
layerChans = sessionPyrLayerInfo.pyrLayerCA1;%CA1 only
%layerChans = [sessionPyrLayerInfo.pyrLayerCA3 sessionPyrLayerInfo.pyrLayerCA1];%CA3+CA1
load([dirs.code 'ripplefilter.mat'])
load([dirs.code 'thetafilter.mat'])
load([dirs.code 'betafilter.mat'])
load([dirs.code 'deltafilter.mat'])

%%%%% add filter description and parameters %%%%%
ripple.descript = ripplefilter.descript;
ripple.kernel = ripplefilter.kernel;
ripple.samprate = ripplefilter.samprate;
smoothing_width_rip = 0.004; % 4 ms
smoothing_kernel_rip = gaussian(smoothing_width_rip*ripple.samprate, ceil(8*smoothing_width_rip*ripple.samprate));
smoothing_width_tdb = 1; % 1 s
smoothing_kernel_tdb = gaussian(smoothing_width_tdb*ripple.samprate, ceil(8*smoothing_width_tdb*ripple.samprate));
minRipDur = round(params.ripple.minRipDur * ripple.samprate);
tmpDat = rawDataBySessionNeural.lfpData;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% clean square-wave, lick related noise from raw data %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% get channel from original data to find noise %%%%%
noiseChan = layerChans(end)+110;%110 channels above top CA1 pyr layer chan
tmpDat = rawDataBySessionNeural.lfpData;
tmpDatBandFilt = bandpass(tmpDat(noiseChan,:),[100 300],params.lfp_samprate_down);
%rip filter
tmpRipDat = filtfilt(ripplefilter.kernel, 1 , tmpDat(noiseChan,:));

%%%%% find noise in the channel %%%%%
%peaks based on bandpass filtered data
[~,tmpPksInds] = findpeaks(tmpDatBandFilt,"MinPeakHeight",0.005);
tmpPksInds = tmpPksInds(diff(tmpPksInds)>20);
%troughs based on bandpass filtered data
tmpMnDat = movmean(tmpDat(noiseChan,:),100);
tmpMidInds = [];
[~,tmpTroughsInds]  = findpeaks(-tmpMnDat,'MinPeakProminence',0.01);
for thisI = 1:length(tmpTroughsInds)-1
    tmpMidInds(thisI) = round(tmpTroughsInds(thisI) + ((tmpTroughsInds(thisI+1)-tmpTroughsInds(thisI))/2));
end

%%%%% smooth noise with a line in all channels %%%%%
tmpDat2 = tmpDat;
for chan = 1:size(tmpDat2,1)
    for mInd = 1:length(tmpMidInds)-1
        tmpInds = tmpMidInds(mInd):tmpMidInds(mInd+1);
        if any(ismember(tmpTroughsInds,tmpInds)) && sum(ismember(tmpPksInds,tmpInds)) == 2%1 trough, 2 peaks
           pk2trough = abs(tmpTroughsInds(find(ismember(tmpTroughsInds,tmpInds))) - tmpPksInds(find(ismember(tmpPksInds,tmpInds),2,'first')));
           if pk2trough < 150%both peaks close to trough             
                tmpLinInds = tmpPksInds(find(ismember(tmpPksInds,tmpInds)));
                tmpDat2(chan,tmpLinInds(1):tmpLinInds(2)) = linspace(tmpDat2(chan,tmpLinInds(1)),tmpDat2(chan,tmpLinInds(2)), length(tmpDat2(chan,tmpLinInds(1):tmpLinInds(2))));
           else%one peak far from trough
               tmpPksIndsByTrough = tmpPksInds(find(ismember(tmpPksInds,tmpInds),2,'first'));
               tmpPksIndsByTrough = tmpPksIndsByTrough(pk2trough<150);
               if tmpTroughsInds(find(ismember(tmpTroughsInds,tmpInds))) > tmpPksIndsByTrough%trough before close peak
                   tmpLinInds = (tmpTroughsInds(find(ismember(tmpTroughsInds,tmpInds))) - pk2trough(pk2trough<150)) : (tmpTroughsInds(find(ismember(tmpTroughsInds,tmpInds))) + pk2trough(pk2trough<150) + 40);
                   tmpDat2(chan,tmpLinInds) = linspace(tmpDat2(chan,tmpLinInds(1)),tmpDat2(chan,tmpLinInds(end)), length(tmpDat2(chan,tmpLinInds)));
               elseif tmpTroughsInds(find(ismember(tmpTroughsInds,tmpInds))) < tmpPksIndsByTrough%trough after close peak
                   tmpLinInds = (tmpTroughsInds(find(ismember(tmpTroughsInds,tmpInds))) - pk2trough(pk2trough<150) -40) : (tmpTroughsInds(find(ismember(tmpTroughsInds,tmpInds))) + pk2trough(pk2trough<150));
                   tmpDat2(chan,tmpLinInds) = linspace(tmpDat2(chan,tmpLinInds(1)),tmpDat2(chan,tmpLinInds(end)), length(tmpDat2(chan,tmpLinInds)));
               end%peak before/after trough
           end
        elseif any(ismember(tmpTroughsInds,tmpInds)) && sum(ismember(tmpPksInds,tmpInds)) == 4%1 trough, 4 peaks
            tmpLinInds = tmpPksInds(find(ismember(tmpPksInds,tmpInds)));
            tmpDat2(chan,tmpLinInds(1):tmpLinInds(2)) = linspace(tmpDat2(chan,tmpLinInds(1)),tmpDat2(chan,tmpLinInds(2)), length(tmpDat2(chan,tmpLinInds(1):tmpLinInds(2))));
            tmpDat2(chan,tmpLinInds(3):tmpLinInds(4)) = linspace(tmpDat2(chan,tmpLinInds(3)),tmpDat2(chan,tmpLinInds(4)), length(tmpDat2(chan,tmpLinInds(3):tmpLinInds(4))));
        elseif any(ismember(tmpTroughsInds,tmpInds)) && sum(ismember(tmpPksInds,tmpInds)) == 1%1 trough, 1 peak
            pk2trough = abs(tmpTroughsInds(find(ismember(tmpTroughsInds,tmpInds))) - tmpPksInds(find(ismember(tmpPksInds,tmpInds))));
            if tmpTroughsInds(find(ismember(tmpTroughsInds,tmpInds))) > tmpPksInds(find(ismember(tmpPksInds,tmpInds))) && pk2trough < 150%trough before 1 peak
                tmpLinInds = (tmpTroughsInds(find(ismember(tmpTroughsInds,tmpInds))) - pk2trough) : (tmpTroughsInds(find(ismember(tmpTroughsInds,tmpInds))) + pk2trough + 40);
                tmpDat2(chan,tmpLinInds) = linspace(tmpDat2(chan,tmpLinInds(1)),tmpDat2(chan,tmpLinInds(end)), length(tmpDat2(chan,tmpLinInds)));
            elseif tmpTroughsInds(find(ismember(tmpTroughsInds,tmpInds))) < tmpPksInds(find(ismember(tmpPksInds,tmpInds))) && pk2trough < 150%trough after 1 peak
                tmpLinInds = (tmpTroughsInds(find(ismember(tmpTroughsInds,tmpInds))) - pk2trough -40) : (tmpTroughsInds(find(ismember(tmpTroughsInds,tmpInds))) + pk2trough);
                tmpDat2(chan,tmpLinInds) = linspace(tmpDat2(chan,tmpLinInds(1)),tmpDat2(chan,tmpLinInds(end)), length(tmpDat2(chan,tmpLinInds)));
            end%peak before/after trough
        elseif any(ismember(tmpTroughsInds,tmpInds)) && sum(ismember(tmpPksInds,tmpInds)) == 3%1 trough, 3 peak
           pk2trough = abs(tmpTroughsInds(find(ismember(tmpTroughsInds,tmpInds))) - tmpPksInds(find(ismember(tmpPksInds,tmpInds),3,'first')));
           tmpPksIndsByTrough = tmpPksInds(find(ismember(tmpPksInds,tmpInds),3,'first'));
           tmpLinInds = tmpPksIndsByTrough(find(abs(diff(pk2trough))==min(abs(diff(pk2trough))))):tmpPksIndsByTrough(find(abs(diff(pk2trough))==min(abs(diff(pk2trough))))+1);
           tmpDat2(chan,tmpLinInds) = linspace(tmpDat2(chan,tmpLinInds(1)),tmpDat2(chan,tmpLinInds(end)), length(tmpDat2(chan,tmpLinInds)));
        end%number of peaks
    end%middle index
end%chan
%rip filter

%tmpRipDat2 = filtfilt(ripplefilter.kernel, 1 , tmpDat2(noiseChan,:));

if plotRawTraceForRipples%DC check if this plot is needed/useful
    %%%% plot %%%%%
    %noise chan
    figure; hold on
    % plot(tmpDat(noiseChan,:))
    % plot(tmpDat2(noiseChan,:))
    % plot(tmpPksInds, tmpDat2(noiseChan,tmpPksInds), '*c')
    % plot(tmpTroughsInds, tmpDat2(noiseChan,tmpTroughsInds), '*b')
    % plot(tmpMidInds, tmpDat2(noiseChan,tmpMidInds), '*r')
    % plot(tmpRipDat)
    % plot(tmpRipDat2)
    %random HIP chan
    tmpRipDatHIP = filtfilt(ripplefilter.kernel, 1 , tmpDat(layerChans(5),:));
    tmpRipDatHIP2 = filtfilt(ripplefilter.kernel, 1 , tmpDat(layerChans(10),:));
    plot(tmpDat(layerChans(5),:)+0.5)
    plot(tmpDat(layerChans(10),:)+0.5)
    % plot(tmpPksInds, tmpDat(layerChans(10),tmpPksInds)+0.5, '*c')
    % plot(tmpTroughsInds, tmpDat(layerChans(10),tmpTroughsInds)+0.5, '*b')
    % plot(tmpMidInds, tmpDat(layerChans(10),tmpMidInds)+0.5, '*r')
    plot(tmpRipDatHIP+0.5)
    plot(tmpRipDatHIP2+0.5)
end%if plotRawTraceForRipples

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% find ripples per channel %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ch = layerChans
    chCtr = chCtr+1;
    %ripple filter (150-250 Hz)
    ripple.filtdata(chCtr,:) = filtfilt(ripplefilter.kernel, 1 , tmpDat(ch,:));

    %Hilbert transform ripple data
    hdata = hilbert(ripple.filtdata(chCtr,:));
    ripple.phase(chCtr,:) = angle(hdata);
    ripple.env(chCtr,:) = abs(hdata);
    ripple.env(chCtr,:) = smoothvect(ripple.env(chCtr,:), smoothing_kernel_rip);%smooth env
    clear hdata

    for t = 1:length(params.ripple.nstdEnv)%loop through detection thresholds

        %find baseline, std, and thresh for this channel's ripple env
        envBaseline = mean(ripple.env(chCtr,:));
        envStd = std(ripple.env(chCtr,:));
        envThresh = envBaseline + params.ripple.nstdEnv(t) * envStd;

        %extract the events if this is a valid trace
        if (envThresh > 0) && any(find(ripple.env(chCtr,:)<envBaseline))
            %find ripple events
            tmprip = [];
            tmprip = extractevents(ripple.env(chCtr,:), envThresh, envBaseline, 0, minRipDur, 0)';
            %start, middle (of energy), and end indices and times
            ripples(chCtr,t).startind = tmprip(:,1);
            ripples(chCtr,t).midind = tmprip(:,8);
            ripples(chCtr,t).endind = tmprip(:,2);
            ripples(chCtr,t).starttime = rawDataBySessionNeural.lfpTime(1) + ripples(chCtr).startind / ripple.samprate;
            ripples(chCtr,t).midtime = rawDataBySessionNeural.lfpTime(1) + ripples(chCtr).midind / ripple.samprate;
            ripples(chCtr,t).endtime = rawDataBySessionNeural.lfpTime(1) + ripples(chCtr).endind / ripple.samprate;
            %ripple characteristics
            ripples(chCtr,t).peak = tmprip(:,3);
            ripples(chCtr,t).energy = tmprip(:,7);
            ripples(chCtr,t).maxthresh = (tmprip(:,9) - envBaseline) / envStd;
        else
            ripples(chCtr,t).startind = [];
            ripples(chCtr,t).midind = [];
            ripples(chCtr,t).endind = [];
            ripples(chCtr,t).starttime = [];
            ripples(chCtr,t).midtime = [];
            ripples(chCtr,t).endtime = [];
            ripples(chCtr,t).peak = [];
            ripples(chCtr,t).energy = [];
            ripples(chCtr,t).maxthresh = [];
        end%> thresh + baseline

        %add other info to ripples struct
        ripples(chCtr,t).timeRange = [0 length(ripple.env(chCtr,:))/ripple.samprate] + rawDataBySessionNeural.lfpTime(1);
        ripples(chCtr,t).samprate = ripple.samprate;
        ripples(chCtr,t).envThreshold = envThresh;
        ripples(chCtr,t).envBaseline = envBaseline;
        ripples(chCtr,t).envStd = envStd;
        ripples(chCtr,t).minimum_duration = params.ripple.minRipDur;
        ripples(chCtr,t).chan = ch;

    end%thresholds

end%ch

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% collect ripples across channels %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% concatenate ripples across channels %%%%%
allRipMidInd1stSD = []; allRipMidInd2ndSD = []; allRipCh1stSD = []; allRipCh2ndSD = [];
for ch = 1:length(layerChans)
    %3 SD Threshold
    allRipMidInd1stSD = [allRipMidInd1stSD; ripples(ch,1).midind];
    allRipCh1stSD = [allRipCh1stSD; transpose(ones(1,length(ripples(ch,1).midind)))*ch];
    %5 SD Threshold
    allRipMidInd2ndSD = [allRipMidInd2ndSD; ripples(ch,2).midind];
    allRipCh2ndSD = [allRipCh2ndSD; transpose(ones(1,length(ripples(ch,2).midind)))*ch];
end%ch

%%%%% sort the ripples %%%%%
[allRipMidInd1stSD, I1stSD] = sort(allRipMidInd1stSD);
allRipCh1stSD = allRipCh1stSD(I1stSD);
[allRipMidInd2ndSD, I2ndSD] = sort(allRipMidInd2ndSD);
allRipCh2ndSD = allRipCh2ndSD(I2ndSD);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% find good ripples across channels %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% reduce to unique ripples detected at 2nd (potentially higher) SD %%%%%
allUniqRipMidInd2ndSD = allRipMidInd2ndSD(diff(allRipMidInd2ndSD)>(params.ripple.minRipDur * params.lfp_samprate_down * 4));
allUniqRipCh2ndSD = allRipCh2ndSD(diff(allRipMidInd2ndSD)>(params.ripple.minRipDur * params.lfp_samprate_down * 4));

%%%%% find ripples detected at 2nd SD close to ripples detected at 1st SD on 2+ other channels %%%%%
allUniqRipMidInd2ndSDGood = []; allUniqRipCh2ndSDGood = []; goodRipCtr = 0;
for r = 1:length(allUniqRipMidInd2ndSD)
    if sum(allRipMidInd1stSD<allUniqRipMidInd2ndSD(r)+(params.ripple.minRipDur * params.lfp_samprate_down * 4) &...
            allRipMidInd1stSD>allUniqRipMidInd2ndSD(r)-(params.ripple.minRipDur * params.lfp_samprate_down * 4)) >= 3%3 because this ripple + 2 more
        goodRipCtr = goodRipCtr + 1;
        allUniqRipMidInd2ndSDGood(goodRipCtr) = allUniqRipMidInd2ndSD(r);
        allUniqRipCh2ndSDGood(goodRipCtr) = allUniqRipCh2ndSD(r);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% collect all info for good ripples %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ripplesGood = [];
goodRipCtrPerCh = zeros(1,length(layerChans));
for r = 1:length(allUniqRipMidInd2ndSDGood)
    goodRipCtrPerCh(allUniqRipCh2ndSDGood(r)) = goodRipCtrPerCh(allUniqRipCh2ndSDGood(r))+1;
    goodRip = find(ripples(allUniqRipCh2ndSDGood(r),2).midind==allUniqRipMidInd2ndSDGood(r));
    ripplesGood(allUniqRipCh2ndSDGood(r)).startind(goodRipCtrPerCh(allUniqRipCh2ndSDGood(r))) = ripples(allUniqRipCh2ndSDGood(r),2).startind(goodRip);
    ripplesGood(allUniqRipCh2ndSDGood(r)).midind(goodRipCtrPerCh(allUniqRipCh2ndSDGood(r))) = ripples(allUniqRipCh2ndSDGood(r),2).midind(goodRip);
    ripplesGood(allUniqRipCh2ndSDGood(r)).endind(goodRipCtrPerCh(allUniqRipCh2ndSDGood(r))) = ripples(allUniqRipCh2ndSDGood(r),2).endind(goodRip);
    ripplesGood(allUniqRipCh2ndSDGood(r)).starttime(goodRipCtrPerCh(allUniqRipCh2ndSDGood(r))) = ripples(allUniqRipCh2ndSDGood(r),2).starttime(goodRip);
    ripplesGood(allUniqRipCh2ndSDGood(r)).midtime(goodRipCtrPerCh(allUniqRipCh2ndSDGood(r))) = ripples(allUniqRipCh2ndSDGood(r),2).midtime(goodRip);
    ripplesGood(allUniqRipCh2ndSDGood(r)).endtime(goodRipCtrPerCh(allUniqRipCh2ndSDGood(r))) = ripples(allUniqRipCh2ndSDGood(r),2).endtime(goodRip);
    ripplesGood(allUniqRipCh2ndSDGood(r)).peak(goodRipCtrPerCh(allUniqRipCh2ndSDGood(r))) = ripples(allUniqRipCh2ndSDGood(r),2).peak(goodRip);
    ripplesGood(allUniqRipCh2ndSDGood(r)).energy(goodRipCtrPerCh(allUniqRipCh2ndSDGood(r))) = ripples(allUniqRipCh2ndSDGood(r),2).energy(goodRip);
    ripplesGood(allUniqRipCh2ndSDGood(r)).maxthresh(goodRipCtrPerCh(allUniqRipCh2ndSDGood(r))) = ripples(allUniqRipCh2ndSDGood(r),2).maxthresh(goodRip);
    ripplesGood(allUniqRipCh2ndSDGood(r)).timeRange = ripples(allUniqRipCh2ndSDGood(r),2).timeRange;
    ripplesGood(allUniqRipCh2ndSDGood(r)).samprate = ripples(allUniqRipCh2ndSDGood(r),2).samprate;
    ripplesGood(allUniqRipCh2ndSDGood(r)).envThreshold = ripples(allUniqRipCh2ndSDGood(r),2).envThreshold;
    ripplesGood(allUniqRipCh2ndSDGood(r)).envBaseline = ripples(allUniqRipCh2ndSDGood(r),2).envBaseline;
    ripplesGood(allUniqRipCh2ndSDGood(r)).envStd = ripples(allUniqRipCh2ndSDGood(r),2).envStd;
    ripplesGood(allUniqRipCh2ndSDGood(r)).minimum_duration = ripples(allUniqRipCh2ndSDGood(r),2).minimum_duration;
    ripplesGood(allUniqRipCh2ndSDGood(r)).chan = ripples(allUniqRipCh2ndSDGood(r),2).chan;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% appply ripple criteria (optional) %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ripplesBad = [];
if params.ripple.applyCriteria
    chCtr = 0;
    for ch = 1:length(ripplesGood)
        chCtr = chCtr+1;
        if any(ripplesGood(chCtr).midind)%if channel has ripples

            %find time around ripple middle (2 s total)
            ripperiods = [ripplesGood(chCtr).midind'-params.ripple.timeAroundRip*ripple.samprate, ripplesGood(chCtr).midind'+params.ripple.timeAroundRip*ripple.samprate];
            ripperiods(ripperiods<=0) = 1;

            %initilize variables
            highFreqVals = [];
            highFreqRatio = [];
            meantdb = [];
            excludeForHighFreq = zeros(1,length(ripplesGood(chCtr).midind));
            excludeForMeanTDB = zeros(1,length(ripplesGood(chCtr).midind));
            excludeForSpeed = zeros(1,length(ripplesGood(chCtr).midind));
            excludeForRawDataNoise = zeros(1,length(ripplesGood(chCtr).midind));
            excludeForNoSpikesPresent = zeros(1,length(ripplesGood(chCtr).midind));
            ripSpikes = nan(100,length(rawDataBySessionNeural.apData),50);

            %find theta/delta+beta ratio for this channel
            if params.ripple.applyTDBRatio
                % theta filter (6-10 Hz)
                theta.filtdata(chCtr,:) = filtfilt(thetafilter.tf.num, 1 , rawDataBySessionNeural.lfpData(ch,:));
                % beta filter (12-30 Hz)
                beta.filtdata(chCtr,:) = filtfilt(betafilter.tf.num, 1 , rawDataBySessionNeural.lfpData(ch,:));
                % delta filter (1-4 Hz)
                delta.filtdata(chCtr,:) = filtfilt(deltafilter.tf.num, 1 , rawDataBySessionNeural.lfpData(ch,:));
                tdbratio.filtdata(chCtr,:) = theta.filtdata(chCtr,:) ./ (delta.filtdata(chCtr,:)+beta.filtdata(chCtr,:));
                tdbratio.filtdata(chCtr,:) = smoothvect(tdbratio.filtdata(chCtr,:), smoothing_kernel_tdb);%smooth data
                tdbratio.baseline(chCtr) = mean(tdbratio.filtdata(chCtr,:));
            end%TDB ratio

            %find all outlier indicies for this channel
            thisChOutInd = [];
            for outInd = 1:size(rawDataBySessionNeural.lfpOutlierInd,1)
                thisChOutInd = [thisChOutInd rawDataBySessionNeural.lfpOutlierInd(outInd,1):rawDataBySessionNeural.lfpOutlierInd(outInd,2)];
            end

            %loop through ripples
            for r = 1:length(ripplesGood(chCtr).midind)

                %threshold for ratio of high frequencies (150-250/250-350 Hz)
                if params.ripple.applyHighFreqRatio
                    % compute psd
                    [Pxx,F] = pwelch(detrend(rawDataBySessionNeural.lfpData(ch,ripperiods(r,1):ripperiods(r,2))),[],[],[],ripple.samprate);
                    % compute max & mean power in 150-250 and 250-350
                    subfreqs{1} = find(F>=params.ripple.freqNumerator(1) & F<=params.ripple.freqNumerator(2)); %ripple band
                    subfreqs{2} = find(F>params.ripple.freqDenominator(1) & F<=params.ripple.freqDenominator(2)); %above ripple, prob movement
                    for f = 1:2 %for 150-250 or 250-400
                        highFreqVals(1,f) = mean(Pxx(subfreqs{f}));
                    end%f
                    highFreqRatio(r) = highFreqVals(1)./highFreqVals(2);
                    excludeForHighFreq(r) = highFreqRatio(r)<params.ripple.ratioThresh;
                end%high freq ratio

                %average theta/(delta+beta) ratio
                if params.ripple.applyTDBRatio
                    meantdb(r) = mean(tdbratio.filtdata(chCtr,ripperiods(r,1):ripperiods(r,2)));%average tdbratio during time around center of rip
                    excludeForMeanTDB(r) = meantdb(r)>tdbratio.baseline(chCtr);
                end%TDB ratio

                %speed
                if params.ripple.applySpeed
                    meanRipSpeed(r) = mean(rawDataBySessionNeural.speed(find(rawDataBySessionNeural.lfpTime>=ripperiods(r,1) & rawDataBySessionNeural.lfpTime<=ripperiods(r,2))));
                    excludeForSpeed(r) = meanRipSpeed(r)>params.speedTh;
                end%speed

                %noise in raw data
                % remove ripple if any index lies within an outlier period
                if any(ismember(ripperiods(r,1):ripperiods(r,2),thisChOutInd))
                    excludeForRawDataNoise(r) = 1;
                end

                % noise deflection in raw data; find differences in raw data values above a threshold
                if any(abs(diff(rawDataBySessionNeural.lfpData(ch,ripperiods(r,1):ripperiods(r,2)))) >= ...
                        std(diff(rawDataBySessionNeural.lfpData(ch,ripperiods(r,1):ripperiods(r,2))))*params.ripple.nstdNoise)
                    excludeForRawDataNoise(r) = 1;
                end%noise deflection
                % % %plot criterion
                % % figure
                % % plot(rawDataBySessionNeural.lfpData(ch,ripperiods(r,1):ripperiods(r,2)))
                % % figure; hold on
                % % plot(diff(rawDataBySessionNeural.lfpData(ch,ripperiods(r,1):ripperiods(r,2))))
                % % plot(ones(1,length(rawDataBySessionNeural.lfpData(ch,ripperiods(r,1):ripperiods(r,2)))) * mean(diff(rawDataBySessionNeural.lfpData(ch,ripperiods(r,1):ripperiods(r,2)))))
                % % plot(ones(1,length(rawDataBySessionNeural.lfpData(ch,ripperiods(r,1):ripperiods(r,2)))) * std(diff(rawDataBySessionNeural.lfpData(ch,ripperiods(r,1):ripperiods(r,2))))*params.ripple.nstdNoise)
                % % plot(ones(1,length(rawDataBySessionNeural.lfpData(ch,ripperiods(r,1):ripperiods(r,2)))) * -std(diff(rawDataBySessionNeural.lfpData(ch,ripperiods(r,1):ripperiods(r,2))))*params.ripple.nstdNoise)

                %spikes present
                % note: looks for spikes from clusters across all channels, not just clusters from the current channel (ch) this loop
                if params.ripple.applyMUA
                    for clu = 1:length(rawDataBySessionNeural.apData)
                        if any(rawDataBySessionNeural.apData(clu).spikeInds/params.samprate >= ripperiods(r,1)/params.lfp_samprate_down &...
                                rawDataBySessionNeural.apData(clu).spikeInds/params.samprate <= ripperiods(r,2)/params.lfp_samprate_down)
                            numSpikes = [];
                            numSpikes = length(find(rawDataBySessionNeural.apData(clu).spikeInds/params.samprate >= ripperiods(r,1)/params.lfp_samprate_down &...
                                rawDataBySessionNeural.apData(clu).spikeInds/params.samprate <= ripperiods(r,2)/params.lfp_samprate_down));
                            %save the spikes for each ripple
                            ripSpikes(r,clu,1:numSpikes) = rawDataBySessionNeural.apData(clu).spikeInds(find(rawDataBySessionNeural.apData(clu).spikeInds/params.samprate >= ripperiods(r,1)/params.lfp_samprate_down &...
                                rawDataBySessionNeural.apData(clu).spikeInds/params.samprate <= ripperiods(r,2)/params.lfp_samprate_down));
                        else
                            ripSpikes(r,clu,:) = nan;
                        end%find spikes
                    end%clu
                    if sum(~isnan(ripSpikes(r,:,:)),'all') > 0%spikes
                        excludeForNoSpikesPresent(r) = 0;
                    else%no spikes
                        excludeForNoSpikesPresent(r) = 1;
                    end%if sum(~isnan(ripSpikes(r,:,:)),'all') > 0
                else
                    excludeForNoSpikesPresent(r) = 0;
                end%spikes present

            end%r

            %determine ripples to exclude
            ripplesToExclue = excludeForHighFreq | excludeForMeanTDB | excludeForSpeed | excludeForRawDataNoise | excludeForNoSpikesPresent;
            %save start indices of bad ripples to separate struct
            ripplesBad(chCtr).midind = ripplesGood(chCtr).midind(ripplesToExclue);
            ripplesBad(chCtr).exclusionReason = [excludeForHighFreq(ripplesToExclue); excludeForMeanTDB(ripplesToExclue); ...
                excludeForSpeed(ripplesToExclue);excludeForRawDataNoise(ripplesToExclue); excludeForNoSpikesPresent(ripplesToExclue)];
            %remove ripples that meet exclusion criteria
            ripplesGood(chCtr).startind(ripplesToExclue) = [];
            ripplesGood(chCtr).midind(ripplesToExclue) = [];
            ripplesGood(chCtr).endind(ripplesToExclue) = [];
            ripplesGood(chCtr).starttime(ripplesToExclue) = [];
            ripplesGood(chCtr).midtime(ripplesToExclue) = [];
            ripplesGood(chCtr).endtime(ripplesToExclue) = [];
            ripplesGood(chCtr).peak(ripplesToExclue) = [];
            ripplesGood(chCtr).energy(ripplesToExclue) = [];
            ripplesGood(chCtr).maxthresh(ripplesToExclue) = [];
            %save spikes for this channel
            ripplesGood(chCtr).ripSpikes = ripSpikes(~ripplesToExclue,:,:);

        else

            excludeForHighFreq = [];
            excludeForMeanTDB = [];
            excludeForSpeed = [];
            excludeForRawDataNoise = [];
            excludeForNoSpikesPresent = [];
        end%if any(ripples(ch).midind)

    end%ch

end%apply ripple criteria

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% add to sturct and save data %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rawDataBySessionNeural.ripplesGood = ripplesGood;
rawDataBySessionNeural.ripplesBad = ripplesBad;
filename = [saveNeuralPath '\' 'rawDataBySessionNeural.mat'];
save(filename, 'rawDataBySessionNeural', '-v7.3')

%% extract ripples across laps %%
%Note: Rest sessions do not have lap structs

if isfile([saveNeuralPath '\rawDataByLapNeural.mat'])

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% load session data %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    load([saveNeuralPath '\rawDataByLapNeural.mat'])
    %tmpcode
    if isfield(rawDataByLapNeural,'ripplesGood')
        rawDataByLapNeural = rmfield(rawDataByLapNeural,'ripplesGood');
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% collect ripple info per lap %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for ii = 1:length(rawDataByLapNeural)
        if ~isempty(rawDataByLapNeural(ii).lfpTime)
            for ch = 1:length(rawDataBySessionNeural.ripplesGood)
                lpRipCtr = zeros(3,1);
                for r = 1:length(rawDataBySessionNeural.ripplesGood(ch).startind)
                    if ismember(rawDataBySessionNeural.ripplesGood(ch).startind(r), rawDataByLapNeural(ii).lfpTime(1):rawDataByLapNeural(ii).lfpTime(end))
                        lpRipCtr(1) = lpRipCtr(1)+1;
                        rawDataByLapNeural(ii).ripplesGood(ch).startind(lpRipCtr(1)) = rawDataBySessionNeural.ripplesGood(ch).startind(r);
                        rawDataByLapNeural(ii).ripplesGood(ch).midind(lpRipCtr(1)) = rawDataBySessionNeural.ripplesGood(ch).midind(r);
                        rawDataByLapNeural(ii).ripplesGood(ch).endind(lpRipCtr(1)) = rawDataBySessionNeural.ripplesGood(ch).endind(r);
                        rawDataByLapNeural(ii).ripplesGood(ch).starttime(lpRipCtr(1)) = rawDataBySessionNeural.ripplesGood(ch).starttime(r);
                        rawDataByLapNeural(ii).ripplesGood(ch).midtime(lpRipCtr(1)) = rawDataBySessionNeural.ripplesGood(ch).midtime(r);
                        rawDataByLapNeural(ii).ripplesGood(ch).endtime(lpRipCtr(1)) = rawDataBySessionNeural.ripplesGood(ch).endtime(r);
                        rawDataByLapNeural(ii).ripplesGood(ch).peak(lpRipCtr(1)) = rawDataBySessionNeural.ripplesGood(ch).peak(r);
                        rawDataByLapNeural(ii).ripplesGood(ch).energy(lpRipCtr(1)) = rawDataBySessionNeural.ripplesGood(ch).energy(r);
                        rawDataByLapNeural(ii).ripplesGood(ch).maxthresh(lpRipCtr(1)) = rawDataBySessionNeural.ripplesGood(ch).maxthresh(r);
                        rawDataByLapNeural(ii).ripplesGood(ch).ripSpikes(lpRipCtr(1),:,:) = rawDataBySessionNeural.ripplesGood(ch).ripSpikes(r,:,:);
                    end%if ismember
                end%r

                if lpRipCtr(1) > 0
                    rawDataByLapNeural(ii).ripplesGood(ch).timeRange = rawDataBySessionNeural.ripplesGood(ch).timeRange;
                    rawDataByLapNeural(ii).ripplesGood(ch).samprate = rawDataBySessionNeural.ripplesGood(ch).samprate;
                    rawDataByLapNeural(ii).ripplesGood(ch).envThreshold = rawDataBySessionNeural.ripplesGood(ch).envThreshold;
                    rawDataByLapNeural(ii).ripplesGood(ch).envBaseline = rawDataBySessionNeural.ripplesGood(ch).envBaseline;
                    rawDataByLapNeural(ii).ripplesGood(ch).envStd =  rawDataBySessionNeural.ripplesGood(ch).envStd;
                    rawDataByLapNeural(ii).ripplesGood(ch).minimum_duration = rawDataBySessionNeural.ripplesGood(ch).minimum_duration;
                    rawDataByLapNeural(ii).ripplesGood(ch).chan = rawDataBySessionNeural.ripplesGood(ch).chan;
                else
                    rawDataByLapNeural(ii).ripplesGood(ch).startind = [];
                    rawDataByLapNeural(ii).ripplesGood(ch).midind = [];
                    rawDataByLapNeural(ii).ripplesGood(ch).endind = [];
                    rawDataByLapNeural(ii).ripplesGood(ch).starttime = [];
                    rawDataByLapNeural(ii).ripplesGood(ch).midtime = [];
                    rawDataByLapNeural(ii).ripplesGood(ch).endtime = [];
                    rawDataByLapNeural(ii).ripplesGood(ch).peak = [];
                    rawDataByLapNeural(ii).ripplesGood(ch).energy = [];
                    rawDataByLapNeural(ii).ripplesGood(ch).maxthresh = [];
                    rawDataByLapNeural(ii).ripplesGood(ch).ripSpikes = [];
                    rawDataByLapNeural(ii).ripplesGood(ch).timeRange = [];
                    rawDataByLapNeural(ii).ripplesGood(ch).samprate = [];
                    rawDataByLapNeural(ii).ripplesGood(ch).envThreshold = [];
                    rawDataByLapNeural(ii).ripplesGood(ch).envBaseline = [];
                    rawDataByLapNeural(ii).ripplesGood(ch).envStd = [];
                    rawDataByLapNeural(ii).ripplesGood(ch).minimum_duration = [];
                    rawDataByLapNeural(ii).ripplesGood(ch).chan = [];
                end%~isempty

            end%ch
        end%if ~isempty(rawDataByLapNeural(ii).lfpTime)
    end%lap

    %%%%%%%%%%%%%%%%%%%%%
    %%%%% save data %%%%%
    %%%%%%%%%%%%%%%%%%%%
    filename = [saveNeuralPath '\' 'rawDataByLapNeural.mat'];
    save(filename, 'rawDataByLapNeural', '-v7.3')

end%if isfile([saveNeuralPath '\rawDataByLapNeural.mat'])

%% extract ripples across trials %%
%Note: Rest sessions and active sessions with <= 1 trial do not have trial structs

if isfile([saveNeuralPath '\rawDataByTrialNeural.mat'])

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% load session data %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    load([saveNeuralPath '\rawDataByTrialNeural.mat'])
    %tmpcode
    for znType = 1:size(rawDataByTrialNeural,1) %1=reward, 2=nonreward, 3 = alt nonreward
        for znNum = 1:size(rawDataByTrialNeural,2)
            if isfield(rawDataByTrialNeural{znType,znNum},'ripplesGood')
                rawDataByTrialNeural{znType,znNum} = rmfield(rawDataByTrialNeural{znType,znNum},'ripplesGood');
            end
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% collect ripple info per trial %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    lpRipCtr = zeros(3,100,20);
    for znType = 1:size(rawDataByTrialNeural,1) %1=reward, 2=nonreward, 3 = alt nonreward
            for lp = 1:size(rawDataByTrialNeural{znType,end},2)

                if isempty(rawDataByTrialNeural{znType,end})
                    continue
                end
            
                if ~isempty(rawDataByTrialNeural{znType,end}(lp).lfpTime)
                    for ch = 1:length(rawDataBySessionNeural.ripplesGood)
                        
                        for r = 1:length(rawDataBySessionNeural.ripplesGood(ch).startind)
                            if ismember(rawDataBySessionNeural.ripplesGood(ch).startind(r), rawDataByTrialNeural{znType,end}(lp).lfpTime(1):rawDataByTrialNeural{znType,end}(lp).lfpTime(end))
                                lpRipCtr(znType,lp,ch) = lpRipCtr(znType,lp,ch)+1;
                                rawDataByTrialNeural{znType,end}(lp).ripplesGood(ch).startind(lpRipCtr(znType,lp,ch)) = rawDataBySessionNeural.ripplesGood(ch).startind(r);
                                rawDataByTrialNeural{znType,end}(lp).ripplesGood(ch).midind(lpRipCtr(znType,lp,ch)) = rawDataBySessionNeural.ripplesGood(ch).midind(r);
                                rawDataByTrialNeural{znType,end}(lp).ripplesGood(ch).endind(lpRipCtr(znType,lp,ch)) = rawDataBySessionNeural.ripplesGood(ch).endind(r);
                                rawDataByTrialNeural{znType,end}(lp).ripplesGood(ch).starttime(lpRipCtr(znType,lp,ch)) = rawDataBySessionNeural.ripplesGood(ch).starttime(r);
                                rawDataByTrialNeural{znType,end}(lp).ripplesGood(ch).midtime(lpRipCtr(znType,lp,ch)) = rawDataBySessionNeural.ripplesGood(ch).midtime(r);
                                rawDataByTrialNeural{znType,end}(lp).ripplesGood(ch).endtime(lpRipCtr(znType,lp,ch)) = rawDataBySessionNeural.ripplesGood(ch).endtime(r);
                                rawDataByTrialNeural{znType,end}(lp).ripplesGood(ch).peak(lpRipCtr(znType,lp,ch)) = rawDataBySessionNeural.ripplesGood(ch).peak(r);
                                rawDataByTrialNeural{znType,end}(lp).ripplesGood(ch).energy(lpRipCtr(znType,lp,ch)) = rawDataBySessionNeural.ripplesGood(ch).energy(r);
                                rawDataByTrialNeural{znType,end}(lp).ripplesGood(ch).maxthresh(lpRipCtr(znType,lp,ch)) = rawDataBySessionNeural.ripplesGood(ch).maxthresh(r);
                                rawDataByTrialNeural{znType,end}(lp).ripplesGood(ch).ripSpikes(lpRipCtr(znType,lp,ch),:,:) = rawDataBySessionNeural.ripplesGood(ch).ripSpikes(r,:,:);
                            end%if ismember
                        end%r
    
                        if lpRipCtr(znType,lp,ch) > 0
                            rawDataByTrialNeural{znType,end}(lp).ripplesGood(ch).timeRange = rawDataBySessionNeural.ripplesGood(ch).timeRange;
                            rawDataByTrialNeural{znType,end}(lp).ripplesGood(ch).samprate = rawDataBySessionNeural.ripplesGood(ch).samprate;
                            rawDataByTrialNeural{znType,end}(lp).ripplesGood(ch).envThreshold = rawDataBySessionNeural.ripplesGood(ch).envThreshold;
                            rawDataByTrialNeural{znType,end}(lp).ripplesGood(ch).envBaseline = rawDataBySessionNeural.ripplesGood(ch).envBaseline;
                            rawDataByTrialNeural{znType,end}(lp).ripplesGood(ch).envStd =  rawDataBySessionNeural.ripplesGood(ch).envStd;
                            rawDataByTrialNeural{znType,end}(lp).ripplesGood(ch).minimum_duration = rawDataBySessionNeural.ripplesGood(ch).minimum_duration;
                            rawDataByTrialNeural{znType,end}(lp).ripplesGood(ch).chan = rawDataBySessionNeural.ripplesGood(ch).chan;
                        else
                            rawDataByTrialNeural{znType,end}(lp).ripplesGood(ch).startind = [];
                            rawDataByTrialNeural{znType,end}(lp).ripplesGood(ch).midind = [];
                            rawDataByTrialNeural{znType,end}(lp).ripplesGood(ch).endind = [];
                            rawDataByTrialNeural{znType,end}(lp).ripplesGood(ch).starttime = [];
                            rawDataByTrialNeural{znType,end}(lp).ripplesGood(ch).midtime = [];
                            rawDataByTrialNeural{znType,end}(lp).ripplesGood(ch).endtime = [];
                            rawDataByTrialNeural{znType,end}(lp).ripplesGood(ch).peak = [];
                            rawDataByTrialNeural{znType,end}(lp).ripplesGood(ch).energy = [];
                            rawDataByTrialNeural{znType,end}(lp).ripplesGood(ch).maxthresh = [];
                            rawDataByTrialNeural{znType,end}(lp).ripplesGood(ch).ripSpikes = [];
                            rawDataByTrialNeural{znType,end}(lp).ripplesGood(ch).timeRange = [];
                            rawDataByTrialNeural{znType,end}(lp).ripplesGood(ch).samprate = [];
                            rawDataByTrialNeural{znType,end}(lp).ripplesGood(ch).envThreshold = [];
                            rawDataByTrialNeural{znType,end}(lp).ripplesGood(ch).envBaseline = [];
                            rawDataByTrialNeural{znType,end}(lp).ripplesGood(ch).envStd = [];
                            rawDataByTrialNeural{znType,end}(lp).ripplesGood(ch).minimum_duration = [];
                            rawDataByTrialNeural{znType,end}(lp).ripplesGood(ch).chan = [];
                        end%~isempty
                    end%ch
                end%~isempty(rawDataByTrialNeural{znType,end}(tr).lfpTime)

            end%tr
    end%znType

    %%%%%%%%%%%%%%%%%%%%%
    %%%%% save data %%%%%
    %%%%%%%%%%%%%%%%%%%%%
    filename = [saveNeuralPath '\' 'rawDataByTrialNeural.mat'];
    save(filename, 'rawDataByTrialNeural', '-v7.3')

end%if isfile([saveNeuralPath '\rawDataByTrialNeural.mat'])

%% plot %%
if plotRipples

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%% find channel with max number of ripples %%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    numGoodRipPerCh = []; numBadRipPerCh = [];
    for ch = 1:length(layerChans)
        numGoodRipPerCh(ch) = length(ripplesGood(ch).midind);
        numBadRipPerCh(ch) = length(ripplesBad(ch).midind);
    end
    [~, maxNumRipChInd] = max(numGoodRipPerCh);

    %%%%%%%%%%%%%%%%
    %%%%% plot %%%%%
    %%%%%%%%%%%%%%%%
    %plot whole raw data and ripple-filtered traces with stars for start of each detected ripple
    figure; hold on
    chSpace = 500;
    p1 = plot(rawDataBySessionNeural.lfpData(layerChans(maxNumRipChInd),1:length(ripple.filtdata(maxNumRipChInd,:)))*500+1*chSpace);%raw trace
    plot(tmpDat2(layerChans(maxNumRipChInd),1:length(ripple.filtdata(maxNumRipChInd,:)))*500+1*chSpace);%noise corrected raw trace
    p2 = plot(ripple.filtdata(maxNumRipChInd,:)*500+1*chSpace);%ripple filtered trace
    p3 = plot(ripple.env(maxNumRipChInd,:)*500+1*chSpace);%ripple env
    p4 = plot(1:length(ripple.filtdata(maxNumRipChInd,:)), ones(1,length(ripple.filtdata(maxNumRipChInd,:)))* ripplesGood(maxNumRipChInd).envBaseline*500+1*chSpace);%env baseline
    p5 = plot(1:length(ripple.filtdata(maxNumRipChInd,:)), ones(1,length(ripple.filtdata(maxNumRipChInd,:)))* ripplesGood(maxNumRipChInd).envThreshold*500+1*chSpace);%env threshold
    for ch = 1:length(layerChans)
        p6 = plot(ripplesGood(ch).midind,ripple.filtdata(maxNumRipChInd,ripplesGood(ch).midind)*500+1*chSpace,'*g');%detected good ripples
        if ~isempty(ripplesBad)
            if ~isempty(ripplesBad(ch).midind)
                p7 = plot(ripplesBad(ch).midind(logical(ripplesBad(ch).exclusionReason(3,:))),...
                    ripple.filtdata(maxNumRipChInd,ripplesBad(ch).midind(logical(ripplesBad(ch).exclusionReason(3,:))))*500+1*chSpace,'*r');%detected speed bad ripples
                p8 = plot(ripplesBad(ch).midind(logical(ripplesBad(ch).exclusionReason(4,:))),...
                    ripple.filtdata(maxNumRipChInd,ripplesBad(ch).midind(logical(ripplesBad(ch).exclusionReason(4,:))))*500+1*chSpace,'*c');%detected raw data bad ripples
                p9 = plot(ripplesBad(ch).midind(logical(ripplesBad(ch).exclusionReason(5,:))),...
                    ripple.filtdata(maxNumRipChInd,ripplesBad(ch).midind(logical(ripplesBad(ch).exclusionReason(5,:))))*500+1*chSpace+10,'*b');%detected MUA bad ripples
            end%if ~isempty(ripplesBad(ch).midind)
        end%if ~isempty(ripplesBad)
    end%ch
    %legend
    legNames = {'RawTrace', 'RipTrace', 'RipEnv', 'EnvBase', 'EnvThresh'};
    if ~isempty(p6)
        legNames{end+1} = 'GoodRip';
    end
    if ~isempty(p7)
        legNames{end+1} = 'SpeedExc';
    end
    if ~isempty(p8)
        legNames{end+1} = 'NoiseExc';
    end
    if ~isempty(p9)
        legNames{end+1} = 'MUAExc';
    end
    leg = legend([p1 p2 p3 p4 p5 p6 p7 p8 p9], legNames,'Box','off',Location='southeast');
    legTitle = get(leg,'title');
    set(legTitle, 'String','Max Channel:', 'FontSize',8);

    % % %plot each layer channel and its respective good ripples 
    % % figure; hold on
    % % for chan = maxNumRipChInd%1:length(layerChans)
    % %     chSpace = 500;
    % %     plot(rawDataBySessionNeural.lfpData(layerChans(chan),1:length(ripple.filtdata(chan,:)))*500+chan*chSpace);%raw trace
    % %     plot(ripple.filtdata(chan,:)*500+chan*chSpace);%ripple filtered trace
    % %     plot(ripple.env(chan,:)*500+chan*chSpace);%ripple env
    % %     if ~isempty(ripplesGood(chan).envBaseline)%if ripples detected in this chan
    % %         plot(1:length(ripple.filtdata(chan,:)), ones(1,length(ripple.filtdata(chan,:)))* ripplesGood(chan).envBaseline*500+chan*chSpace);%env baseline
    % %         plot(1:length(ripple.filtdata(chan,:)), ones(1,length(ripple.filtdata(chan,:)))* ripplesGood(chan).envThreshold*500+chan*chSpace);%env threshold
    % %         plot(ripplesGood(chan).midind,ripple.filtdata(chan,ripplesGood(chan).midind)*500+chan*chSpace,'*g');%detected good ripples
    % %     end%if ~isempty(ripplesGood(chan).envBaseline)
    % % end
    % % 
    % % %plot spiking on top of filtered ripple trace
    % % for ch = maxNumRipChInd%1:length(layerChans)
    % %     for r = 1:length(ripplesGood(ch).startind)
    % %         %gather spikes for all clusters
    % %         spikeTrainLFPTime = round(squeeze(ripplesGood(ch).ripSpikes(r,:,:)) / params.samprate * params.lfp_samprate_down);
    % % 
    % %         %initialize
    % %         thisRipPossSp = repmat(1:ripplesGood(ch).endind(r),size(spikeTrainLFPTime,1),1);%matrix for all possible spike times
    % %         thisRipRealSp = zeros(size(thisRipPossSp,1), size(thisRipPossSp,2));%matrix of zeroes for all real spike times
    % % 
    % %         %find spike times for each cluster
    % %         for clu = 1:size(spikeTrainLFPTime,1)
    % %             thisCluSpInd = ~isnan(spikeTrainLFPTime(clu,:));
    % %             thisRipRealSp(clu,spikeTrainLFPTime(clu,thisCluSpInd)) = ripple.filtdata(ch,spikeTrainLFPTime(clu,thisCluSpInd));
    % %         end
    % %         thisRipRealSp(thisRipRealSp==0) = nan;%turn zeros into NaNs
    % % 
    % %         %(optional) reduce to only putative pyramidal cells roughly in hippocampus
    % %         hipChansOnly = 1;
    % %         if hipChansOnly
    % %             for clu = 1:size(spikeTrainLFPTime,1)
    % %                 if ismember(rawDataBySessionNeural.apData(clu).maxChan,[1:150]) &&...
    % %                         contains(rawDataBySessionNeural.apData(clu).putativeCellType, 'Pyr')
    % %                 else
    % %                     thisRipRealSp(clu,:) = nan;
    % %                 end%if
    % %             end%chan
    % %         end%pyrLayerChansOnly
    % % 
    % %         %plot
    % %         colors = cbrewer('div','RdBu', size(spikeTrainLFPTime,1)+3); colors = flip(colors);
    % %         figure; hold on
    % %         plot(thisRipPossSp(ch,ripplesGood(ch).startind(r):ripplesGood(ch).endind(r)), rawDataBySessionNeural.lfpData(ch,ripplesGood(ch).startind(r):ripplesGood(ch).endind(r)));%raw trace
    % %         plot(thisRipPossSp(ch,ripplesGood(ch).startind(r):ripplesGood(ch).endind(r)), ripple.filtdata(ch,ripplesGood(ch).startind(r):ripplesGood(ch).endind(r)));%filtered ripple trace
    % %         img = arrayfun( @(x) plot(thisRipPossSp(x,ripplesGood(ch).startind(r):ripplesGood(ch).endind(r)),thisRipRealSp(x,ripplesGood(ch).startind(r):ripplesGood(ch).endind(r)),'|',...
    % %             'Color', abs(colors(x,:)), 'DisplayName', num2str(x), 'LineWidth', 5), 1:size(spikeTrainLFPTime,1));%each cluster individual color
    % %     end%ripples
    % % end%ch

end%if plotRipples

end%function
