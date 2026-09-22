function choiceOutput = choiceAnalyzer_DS(allindex, dirs, params, choiceTrialType)
% Created from Stephanie Prince's Python ChoiceAnalyzer.
%
% Main question:
% Can pre-zone behavior predict whether the mouse will lick in the upcoming
% zone on that trial?
%
% Inputs to the LSTM:
%   1. Smoothed lick rate BEFORE zone entry
%   2. Smoothed speed BEFORE zone entry
%
% Binary target:
%   0 = no lick in the focal zone on that trial
%   1 = at least one lick in the focal zone on that trial
%
% IMPORTANT:
% For reward trials (choiceTrialType = 1), the focal zone is the reward
% zone, so the target is exactly "lick in reward zone" versus "no lick in
% reward zone".
%
% If choiceTrialType = 0, reward, control, and alternate-control trials are
% all included. In that mode, licksInZone refers to the focal zone for that
% trial type, because getTrialBehaviorStats_linearDC calculates licksInZone
% separately for reward and control zones.
%
% The model is environment-agnostic. Track/environment identity is saved as
% metadata for later comparisons, but it is NOT given to the LSTM.
%
% choiceTrialType:
%   0 = all available zone trial types
%   1 = reward trials only
%   2 = nonreward/control trials only
%   3 = alternate nonreward trials only
%
% Required toolbox:
%   Deep Learning Toolbox (trainnet, sequenceInputLayer, lstmLayer)

%% Check inputs %%
if nargin < 4 || isempty(choiceTrialType)
    choiceTrialType = 0;
end

if ~ismember(choiceTrialType, [0 1 2 3])
    error('choiceTrialType must be 0, 1, 2, or 3.');
end

if exist('trainnet', 'file') ~= 2
    error(['trainnet was not found. GPT_choiceAnalyzer_DS requires a MATLAB ' ...
        'version with Deep Learning Toolbox and trainnet support.']);
end

%% Choice Analysis parameters %%
% These follow the original ChoiceAnalyzer defaults where possible.
choiceParams.miniBatchSize = 32;
choiceParams.maxEpochs = 20;
choiceParams.initialLearnRate = 0.1;
choiceParams.l2Regularization = 0;
choiceParams.numFolds = 6;
choiceParams.numRepeats = 3;
choiceParams.randomSeed = 21;
choiceParams.numLSTMUnits = 10;

% Trial type labels are metadata only. They are not model inputs.
trialTypeNames = {'reward', 'control', 'altControl'};

if choiceTrialType == 0
    trialTypesToUse = [1 2 3];
    analysisName = 'all';
    fprintf(['Choice Analysis: using all available zone trial types. ' ...
        'For control trials, the target means licking in the control zone.\n']);
elseif choiceTrialType == 1
    trialTypesToUse = 1;
    analysisName = 'reward';
elseif choiceTrialType == 2
    trialTypesToUse = 2;
    analysisName = 'control';
else
    trialTypesToUse = 3;
    analysisName = 'altControl';
end

%% Gather trial data across active sessions %%
% Each inputData cell contains one trial.
% Rows = spatial/time-order bins before the zone.
% Columns = [lickRateSmooth, velocCountsSmooth].
inputData = {};
targetData = [];
trialInfo = struct([]);

fprintf('\nGathering Choice Analysis trials...\n');

for i = 1:size(allindex,1)

    % Only use active VR sessions.
    if size(allindex,2) >= 6 && allindex(i,6) ~= 2
        continue
    end

    %%%%% session info %%%%%
    subj = [params.iden num2str(allindex(i,1))];
    sessDate = num2str(allindex(i,2));
    sessNum = num2str(allindex(i,3));
    trackInfo = allindex(i,4);

    saveBehaviorPath = fullfile(dirs.saveoutputstructs, ...
        'Data', 'Behavior', 'sessionData', subj, ...
        [sessDate '_' sessNum '_' num2str(trackInfo)]);

    statsFile = fullfile(saveBehaviorPath, 'statsByRewardTrial.mat');

    if ~isfile(statsFile)
        fprintf('Skipping %s_%s_%s: statsByRewardTrial.mat not found.\n', ...
            subj, sessDate, sessNum);
        continue
    end

    loadedStats = load(statsFile, 'statsByRewardTrial');

    if ~isfield(loadedStats, 'statsByRewardTrial') || isempty(loadedStats.statsByRewardTrial)
        fprintf('Skipping %s_%s_%s: statsByRewardTrial is empty.\n', ...
            subj, sessDate, sessNum);
        continue
    end

    statsByRewardTrial = loadedStats.statsByRewardTrial;

    % Loop through the requested zone trial types.
    for znType = trialTypesToUse

        % Some environments do not have alternate-control trials.
        if znType > size(statsByRewardTrial,1)
            continue
        end

        % getTrialBehaviorStats_linearDC stores the chronological list of
        % individual trials in the last non-empty cell of each zone-type row.
        nonEmptyCells = find(~cellfun(@isempty, statsByRewardTrial(znType,:)));

        if isempty(nonEmptyCells)
            continue
        end

        combinedCell = nonEmptyCells(end);
        currTrials = statsByRewardTrial{znType, combinedCell};

        if isempty(currTrials) || ~isstruct(currTrials)
            continue
        end

        % The cells before the combined cell correspond to individual zone
        % numbers. This lets us reconstruct the zone number for metadata.
        numZones = combinedCell - 1;

        for tr = 1:numel(currTrials)

            currTrial = currTrials(tr);

            % Required fields for this analysis.
            if ~isfield(currTrial, 'lickRateSmooth') || ...
                    ~isfield(currTrial, 'velocCountsSmooth') || ...
                    ~isfield(currTrial, 'binEdges') || ...
                    ~isfield(currTrial, 'lickBehavior') || ...
                    ~isfield(currTrial.lickBehavior, 'licksInZone')
                continue
            end

            lickRateSmooth = double(currTrial.lickRateSmooth(:));
            speedSmooth = double(currTrial.velocCountsSmooth(:));
            binEdges = double(currTrial.binEdges(:));

            % Use the number of bins shared by both behavioral signals.
            numBins = min(length(lickRateSmooth), length(speedSmooth));

            if numBins == 0 || length(binEdges) < numBins
                continue
            end

            lickRateSmooth = lickRateSmooth(1:numBins);
            speedSmooth = speedSmooth(1:numBins);

            % binEdges are defined relative to the focal zone.
            % Negative bins occur before zone entry; 0 is zone onset.
            preZoneIndex = binEdges(1:numBins) < 0;

            if ~any(preZoneIndex)
                continue
            end

            % Only pre-zone behavior is given to the model. This prevents
            % the LSTM from seeing the lick inside the zone that defines the
            % binary target.
            currInput = [lickRateSmooth(preZoneIndex), speedSmooth(preZoneIndex)];

            % Clean occasional missing values without changing trial length.
            [currInput, usableTrial] = cleanSequence(currInput);

            if ~usableTrial
                continue
            end

            % Binary target: did at least one lick happen in the focal zone?
            currTarget = double(currTrial.lickBehavior.licksInZone > 0);

            inputData{end+1,1} = currInput;
            targetData(end+1,1) = currTarget;

            % Save metadata for later comparisons. None of these fields are
            % supplied to the LSTM.
            thisTrial = length(targetData);
            trialInfo(thisTrial).animal = allindex(i,1);
            trialInfo(thisTrial).date = allindex(i,2);
            trialInfo(thisTrial).sessionNum = allindex(i,3);
            trialInfo(thisTrial).track = trackInfo;
            trialInfo(thisTrial).trialType = znType;
            trialInfo(thisTrial).trialTypeName = trialTypeNames{znType};

            if numZones > 0
                trialInfo(thisTrial).zoneNumber = mod(tr-1, numZones) + 1;
                trialInfo(thisTrial).trialNumberWithinZone = floor((tr-1) / numZones) + 1;
            else
                trialInfo(thisTrial).zoneNumber = NaN;
                trialInfo(thisTrial).trialNumberWithinZone = tr;
            end

        end
    end
end

%% Check gathered data %%
numTrials = length(targetData);

if numTrials == 0
    error('No usable trials were found for Choice Analysis.');
end

numNoLick = sum(targetData == 0);
numLick = sum(targetData == 1);

fprintf('Choice Analysis trials found: %d\n', numTrials);
fprintf('No-lick trials: %d\n', numNoLick);
fprintf('Lick trials: %d\n', numLick);

if numNoLick == 0 || numLick == 0
    error(['Choice Analysis requires both target classes. The selected data ' ...
        'contains only one class.']);
end

% Repeated stratified K-fold requires enough observations in each class.
smallestClass = min(numNoLick, numLick);
numFolds = min(choiceParams.numFolds, smallestClass);

if numFolds < 2
    error('Not enough trials in the smaller target class for cross-validation.');
end

if numFolds < choiceParams.numFolds
    fprintf('Using %d folds instead of %d because the smaller class has %d trials.\n', ...
        numFolds, choiceParams.numFolds, smallestClass);
end

%% Build the LSTM %%
% The target is one outcome per trial, so OutputMode="last" is used.
% This returns one predicted lick probability for the whole pre-zone
% behavioral sequence.
numFeatures = 2;

layers = [
    sequenceInputLayer(numFeatures)
    lstmLayer(choiceParams.numLSTMUnits, OutputMode="last")
    fullyConnectedLayer(1)
    sigmoidLayer
    ];

%% Repeated stratified cross-validation %%
% Each trial is held out once per repetition. Predictions are averaged
% across repetitions to give one cross-validated probability per trial.
repeatPredictions = nan(numTrials, choiceParams.numRepeats);
foldAccuracy = nan(choiceParams.numRepeats, numFolds);
foldLoss = nan(choiceParams.numRepeats, numFolds);
trainingHistory = cell(choiceParams.numRepeats, numFolds);
foldTestIndex = cell(choiceParams.numRepeats, numFolds);

for r = 1:choiceParams.numRepeats

    rng(choiceParams.randomSeed + r - 1);

    % Categorical labels make cvpartition stratify the binary classes.
    cv = cvpartition(categorical(targetData), 'KFold', numFolds);

    fprintf('\nChoice Analysis repetition %d of %d\n', r, choiceParams.numRepeats);

    for k = 1:numFolds

        trainIndex = training(cv, k);
        testIndex = test(cv, k);

        inputTrain = inputData(trainIndex);
        inputTest = inputData(testIndex);
        targetTrain = targetData(trainIndex);
        targetTest = targetData(testIndex);

        % Normalize lick rate and speed using TRAINING data only.
        [inputTrain, inputTest] = normalizeSequenceData(inputTrain, inputTest);

        currBatchSize = min(choiceParams.miniBatchSize, sum(trainIndex));

        options = trainingOptions("adam", ...
            MaxEpochs=choiceParams.maxEpochs, ...
            MiniBatchSize=currBatchSize, ...
            InitialLearnRate=choiceParams.initialLearnRate, ...
            L2Regularization=choiceParams.l2Regularization, ...
            Shuffle="every-epoch", ...
            Verbose=false, ...
            Plots="none");

        [trainedNet, info] = trainnet(inputTrain, targetTrain, layers, ...
            "binary-crossentropy", options);

        prediction = minibatchpredict(trainedNet, inputTest, ...
            MiniBatchSize=currBatchSize);
        prediction = predictionToDouble(prediction);

        predictedClass = double(prediction >= 0.5);
        accuracy = mean(predictedClass == targetTest);
        loss = binaryLoss(targetTest, prediction);

        repeatPredictions(testIndex, r) = prediction;
        foldAccuracy(r,k) = accuracy;
        foldLoss(r,k) = loss;
        trainingHistory{r,k} = info;
        foldTestIndex{r,k} = find(testIndex);

        fprintf('  Fold %d of %d: accuracy = %.3f, loss = %.3f\n', ...
            k, numFolds, accuracy, loss);

    end
end

%% Average repeated hold-out predictions %%
meanPrediction = mean(repeatPredictions, 2, 'omitnan');
predictedClass = double(meanPrediction >= 0.5);
correctPrediction = predictedClass == targetData;
overallAccuracy = mean(correctPrediction);
overallLoss = binaryLoss(targetData, meanPrediction);

fprintf('\nCross-validated Choice Analysis accuracy: %.3f\n', overallAccuracy);
fprintf('Cross-validated Choice Analysis loss: %.3f\n', overallLoss);

%% Train one final model on all selected data %%
% Cross-validated predictions above should be used for performance.
% This final model is saved only so the same trained analyzer can later be
% applied to additional trials if desired.
[allInputNormalized, ~, finalNormalization] = normalizeSequenceData(inputData, {});

finalBatchSize = min(choiceParams.miniBatchSize, numTrials);
finalOptions = trainingOptions("adam", ...
    MaxEpochs=choiceParams.maxEpochs, ...
    MiniBatchSize=finalBatchSize, ...
    InitialLearnRate=choiceParams.initialLearnRate, ...
    L2Regularization=choiceParams.l2Regularization, ...
    Shuffle="every-epoch", ...
    Verbose=false, ...
    Plots="none");

rng(choiceParams.randomSeed);
[finalModel, finalTrainingHistory] = trainnet(allInputNormalized, targetData, ...
    layers, "binary-crossentropy", finalOptions);

%% Put trial information and predictions together %%
trialTable = struct2table(trialInfo);
trialTable.target = targetData;
trialTable.predictedProbability = meanPrediction;
trialTable.predictedClass = predictedClass;
trialTable.correct = correctPrediction;

%% Save output %%
choiceOutput = struct();
choiceOutput.analysisName = analysisName;
choiceOutput.choiceTrialType = choiceTrialType;
choiceOutput.featureNames = {'lickRateSmooth', 'velocCountsSmooth'};
choiceOutput.preZoneOnly = true;
choiceOutput.environmentUsedAsInput = false;
choiceOutput.trialTypeUsedAsInput = false;
choiceOutput.target = targetData;
choiceOutput.inputData = inputData;
choiceOutput.trialTable = trialTable;
choiceOutput.repeatPredictions = repeatPredictions;
choiceOutput.prediction = meanPrediction;
choiceOutput.predictedClass = predictedClass;
choiceOutput.overallAccuracy = overallAccuracy;
choiceOutput.overallLoss = overallLoss;
choiceOutput.foldAccuracy = foldAccuracy;
choiceOutput.foldLoss = foldLoss;
choiceOutput.foldTestIndex = foldTestIndex;
choiceOutput.trainingHistory = trainingHistory;
choiceOutput.finalModel = finalModel;
choiceOutput.finalNormalization = finalNormalization;
choiceOutput.finalTrainingHistory = finalTrainingHistory;
choiceOutput.params = choiceParams;

if choiceTrialType == 1
    choiceOutput.targetDescription = ...
        '1 = at least one lick in reward zone; 0 = no lick in reward zone';
else
    choiceOutput.targetDescription = ...
        ['1 = at least one lick in focal zone; 0 = no lick in focal zone. ' ...
         'Use choiceTrialType = 1 for a strict reward-zone target.'];
end

choiceSavePath = fullfile(dirs.saveoutputstructs, 'Data', 'Behavior', 'ChoiceAnalysis');

if ~isfolder(choiceSavePath)
    mkdir(choiceSavePath)
end

choiceSaveFile = fullfile(choiceSavePath, ['choiceOutput_' analysisName '.mat']);
save(choiceSaveFile, 'choiceOutput', '-v7.3');

fprintf('Choice Analysis saved to: %s\n', choiceSaveFile);

end


%% Helper: clean one behavioral sequence %%
function [sequence, usableTrial] = cleanSequence(sequence)
% Replace isolated missing values by interpolation. If an entire feature is
% missing, the trial cannot be used.

usableTrial = true;

for feature = 1:size(sequence,2)

    currFeature = sequence(:,feature);
    currFeature(~isfinite(currFeature)) = NaN;

    if all(isnan(currFeature))
        usableTrial = false;
        return
    end

    currFeature = fillmissing(currFeature, 'linear', 'EndValues', 'nearest');
    sequence(:,feature) = currFeature;

end

end


%% Helper: normalize using training data only %%
function [trainData, testData, normalization] = normalizeSequenceData(trainData, testData)
% Data format for trainnet sequence input is time-by-features.

allTrainData = vertcat(trainData{:});

featureMean = mean(allTrainData, 1, 'omitnan');
featureStd = std(allTrainData, 0, 1, 'omitnan');

% Prevent division by zero if a feature has no variation in this fold.
featureStd(featureStd == 0 | ~isfinite(featureStd)) = 1;
featureMean(~isfinite(featureMean)) = 0;

for i = 1:length(trainData)
    trainData{i} = (trainData{i} - featureMean) ./ featureStd;
end

for i = 1:length(testData)
    testData{i} = (testData{i} - featureMean) ./ featureStd;
end

normalization.mean = featureMean;
normalization.std = featureStd;

end


%% Helper: convert network prediction to normal MATLAB double column %%
function prediction = predictionToDouble(prediction)

if isa(prediction, 'dlarray')
    prediction = extractdata(prediction);
end

if isa(prediction, 'gpuArray')
    prediction = gather(prediction);
end

prediction = double(prediction(:));

end


%% Helper: binary cross-entropy on held-out predictions %%
function loss = binaryLoss(target, prediction)

epsVal = 1e-15;
prediction = min(max(prediction, epsVal), 1 - epsVal);
target = double(target(:));
prediction = double(prediction(:));

loss = -mean(target .* log(prediction) + ...
    (1 - target) .* log(1 - prediction));

end
