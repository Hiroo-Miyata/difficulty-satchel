close all; clear all; addpath(genpath("./utils/plotting/"));

[dates, rootFolder, axisOutputFileNames, processedFileNames, analysisBin, params] = pp_getParam("Dataset", "BL", "timeperiod", "GC", "neuralDataType", "non_stitched");
ndates = length(dates);

outputFolder = rootFolder + "/projection-2D-eachday/"; makeDir(outputFolder); 
outputFolderProjection = rootFolder + "/rewardAxis-projection-eachday/"; makeDir(outputFolderProjection);
outputFolderParams= rootFolder + "/rewardAxis-params-eachday/"; makeDir(outputFolderParams);
outputFolderStretching = rootFolder + "/stretching/"; makeDir(outputFolderStretching);
outputFolderEachNeuron = rootFolder + "/projection-time-eachday/"; makeDir(outputFolderEachNeuron); 

Yall = cell(ndates, 3, 2, 8);
for d = 1:ndates
date = dates(d);
load(processedFileNames(d));
axisOutputFileName = axisOutputFileNames(d);

% remove calibration part 
trialNums = [trialData.trial];
trialData = trialData(trialNums > 24);
% Get labels
directionLabels = [trialData.directionLabel]; directions = unique(directionLabels); ndirections = length(directions);
rewardLabels = [trialData.rewardLabel]; rewards = unique(rewardLabels); nrewards = length(rewards);
difficultyLabels =[trialData.targetSizeLabel]; difficulties = unique(difficultyLabels); ndifficulties = length(difficulties);
focusDifficultyLabels = [trialData.focusDifficultyLabel]; focusDifficulties = unique(focusDifficultyLabels); nfocusDifficulties = length(focusDifficulties);
focusBlockLabels = [trialData.block]; focusBlocks = unique(focusBlockLabels); nfocusBlocks = length(focusBlocks);
rewColors = [1 0 0; 1 0.6470 0; 0 0 1];
diffColors  = [0 0.447 0.741; 0.466 0.674 0.188]; %tiny huge: blue and green
direColors = {[1 .5 .5],[.75 .75 .5],[.5 1 .5],[.25 .75 .5],[0 .5 .5],[0.25 0.25 .5],[0.5 0 .5],[0.75 0.25 .5]};
trialNums =  [trialData.trial];
delayTimes = [trialData.delayLength]; delayTimes = round(delayTimes / 50) * 50;

% Get data and avg within dir x rew
neuralActivity = cat(3, trialData(:).firingRates);
neuralData = mean(neuralActivity(:, analysisBin, :), 2);
neuralData = squeeze(neuralData)';
[ntrials,nneurons] = size(neuralData);
neuralData_byParameters_mean = nan(ndirections,nrewards,ndifficulties,nneurons);
for l = 1:ndirections
    for r = 1:nrewards
        for j = 1:ndifficulties
            curInds = directionLabels==directions(l) & rewardLabels==rewards(r) & difficultyLabels==difficulties(j);
            neuralData_byParameters_mean(l,r,j,:) = mean(neuralData(curInds,:));
        end; clear j
    end; clear r
end; clear l

% reward axis of both
meanNeuralData = squeeze(nanmean(neuralData_byParameters_mean,[1 3]));
% meanNeuralData = squeeze(meanNeuralData(:, 2, :));
[wR,zR,eigVls_R] = pca(meanNeuralData-mean(meanNeuralData, 1),'numcomponents',nrewards-1);
if zR(nrewards, 1) < zR(1, 1)
    wR(:, 1) = -wR(:, 1);
end

% reward axis of Tiny
meanNeuralData = squeeze(nanmean(neuralData_byParameters_mean,[1]));
meanNeuralData = squeeze(meanNeuralData(:, 1, :));
[wRTiny,zRTiny,eigVls_RTiny] = pca(meanNeuralData-mean(meanNeuralData, 1),'numcomponents',nrewards-1);
if zRTiny(nrewards, 1) < zRTiny(1, 1)
    wRTiny(:, 1) = -wRTiny(:, 1);
end

% reward axis of Huge
meanNeuralData = squeeze(nanmean(neuralData_byParameters_mean,[1]));
meanNeuralData = squeeze(meanNeuralData(:, 2, :));
[wRHuge,zRHuge,eigVls_RHuge] = pca(meanNeuralData-mean(meanNeuralData, 1),'numcomponents',nrewards-1);
if zRHuge(nrewards, 1) < zRHuge(1, 1)
    wRHuge(:, 1) = -wRHuge(:, 1);
end

% difficulty axis
meanNeuralData = squeeze(nanmean(neuralData_byParameters_mean,[1 2]));
[wDif,zDif,eigVls_Dif]  = pca(meanNeuralData-mean(meanNeuralData, 1),'numcomponents',ndifficulties-1);
if zDif(ndifficulties, 1) < zDif(1, 1)
    wDif(:, 1) = -wDif(:, 1);
end

% R-Diff axis
meanNeuralData = squeeze(nanmean(neuralData_byParameters_mean,[1]));
meanNeuralData = reshape(meanNeuralData, [nrewards*ndifficulties, nneurons]);
[wRDiff,zRDiff,eigVls_RDiff] = pca(meanNeuralData-mean(meanNeuralData, 1),'numcomponents',nrewards*ndifficulties-1);
if zRDiff(ndifficulties, 1) < zRDiff(1, 1)
    wRDiff(:, 1) = -wRDiff(:, 1);
end
if zRDiff(nrewards+1, 2) < zRDiff(1, 2)
    wRDiff(:, 2) = -wRDiff(:, 2);
end

meanFiringRates = mean(neuralData, 1);

save(axisOutputFileName, "wR", "wRTiny", "wRHuge", "wDif", "wRDiff", "eigVls_R", "eigVls_RTiny", "eigVls_RHuge", "eigVls_Dif", "eigVls_RDiff", "meanFiringRates");

%% angle of wDif and wR
dot_product = dot(wR(:,1), wDif(:,1));
angle = acos(dot_product) * 180 / pi;
if angle > 90
    angle = 180 - angle;
end
disp(angle);

%% plot reward PC1 and difficulty PC1 by 2D scatter plot
neuralData_onDPC = (neuralData - mean(neuralData, 1)) * wRDiff(:,[1 2]);
figure;
hold on;
for i=1:ndifficulties
    for j=1:nrewards
        curInds = rewardLabels==rewards(j) & difficultyLabels==difficulties(i);
        X(i,j) = mean(neuralData_onDPC(curInds, 1));
        Xerr(i,j) = std(neuralData_onDPC(curInds, 1))/sqrt(sum(curInds));
        Y(i,j) = mean(neuralData_onDPC(curInds, 2));
        Yerr(i,j) = std(neuralData_onDPC(curInds, 2))/sqrt(sum(curInds));
    end
    h(i) = errorbar(X(i,:), Y(i,:), Yerr(i,:), Yerr(i,:), Xerr(i,:), Xerr(i,:), "Color", diffColors(i,:), 'LineWidth', 2.5); hold on;
    scatter(X(i,:), Y(i,:), 100, rewColors, "filled"); hold on;
end
set(gca, 'fontsize', 20, 'fontname', 'arial', 'tickdir', 'out', 'fontweight', 'bold');
ylim([-30 20]); xlim([-30 20]);
xlabel("RDiff Axis: PC1 (Hz): "+round(eigVls_RDiff(1)/sum(eigVls_RDiff)*100, 2)+"%");
ylabel("RDiff Axis: PC2 (Hz): "+round(eigVls_RDiff(2)/sum(eigVls_RDiff)*100, 2)+"%");
legend(h, ["Tiny", "Huge"], Location="best"); set(gcf,'position',[0,0,550,550]);
saveas(gcf, outputFolder+"2Dprojection-"+date+".jpg");

%% plot reward PC1 and difficulty PC1
rAxes = ["Tiny", "Huge", "All"];
for k=1:3 % 1:2
    Y = zeros(3,2); % reward x difficulty
    Yerr = zeros(3,2);
    if k==1
        neuralData_onDPC = (neuralData - mean(neuralData, 1)) * wRTiny(:,1);
    elseif k==2
        neuralData_onDPC = (neuralData - mean(neuralData, 1)) * wRHuge(:,1);
    else
        neuralData_onDPC = (neuralData - mean(neuralData, 1)) * wR(:,1);
    end

    VS_reward_size(neuralData_onDPC', rewardLabels, difficultyLabels, "Label", "Reward Axis (spikes/s)", "OutputFolder", outputFolderProjection+"Axis-in-"+rAxes(k)+"-"+date);
end


neuralData_onDPC = (neuralData - mean(neuralData, 1)) * wR(:,1);


%% Reward Axis
VS_reward_size(neuralData_onDPC', rewardLabels, difficultyLabels, "Label", "Reward Axis (spikes/s)", "OutputFolder", outputFolderParams+date+"-"+"RewardAxis", "Ylim", [-25 20])
VS_delay_time(neuralData_onDPC', rewardLabels, difficultyLabels, delayTimes, "Label", "Reward Axis (spikes/s)", "OutputFolder", outputFolderParams+date+"-"+"RewardAxis")
VS_direction(neuralData_onDPC', rewardLabels, difficultyLabels, directionLabels, "Label", "Reward Axis (spikes/s)", "OutputFolder", outputFolderParams+date+"-"+"RewardAxis")
VS_trial(trialNums, neuralData_onDPC', rewardLabels, difficultyLabels, "Label", "Reward Axis (spikes/s)", "OutputFolder", outputFolderParams+date+"-"+"RewardAxis", ...
            "IsMultipleDays",false)
VS_focus_difficulty(neuralData_onDPC', rewardLabels, difficultyLabels, ones(1, ntrials), focusDifficultyLabels, ...
"Label", "Reward Axis (spikes/s)", "OutputFolder", outputFolderParams+date+"RewardAxis")

%% vs time
% addpath(genpath("./utils/processing/"));
% [neuralActivity_smoothed, ~, ~, ~] = kernelSmooth(double(neuralActivity));
% for k=1:nneurons
%     PSTHs = squeeze(neuralActivity_smoothed(k, :, :));
%     VS_time(PSTHs, rewardLabels, difficultyLabels, "Label", "firing rate (spikes/s)", "OutputFolder", outputFolderEachNeuron+"neuron"+k, ...
%         "Timeperiod", "TO", "Xlim", [-200 600], "DiffLegendPos", "north", "RewLegendPos", "northwest");
% end

%% 
rewardAxisName = "rewardAxisBoth";
neuralData_onDPC = (neuralData - mean(neuralData, 1)) * wR(:,1);
for i=1:nrewards
    for j=1:ndifficulties
        for k=1:ndirections
            curInds = rewardLabels==rewards(i) & difficultyLabels==difficulties(j) & directionLabels==directions(k);
            Yall{d, i, j, k} = neuralData_onDPC(curInds, 1);
        end
    end
end

end

rewardAxisRaw = Yall;
rewardAxisMean = zeros(size(Yall));
rewardAxisStd  = zeros(size(Yall));
for d=1:ndates
    for j=1:ndifficulties
        for i=1:nrewards
            for k=1:ndirections
                Ys = Yall{d,i,j,k};
                rewardAxisMean(d,i,j,k) = mean(Ys);
                rewardAxisStd(d,i,j,k)  = std(Ys) / sqrt(length(Ys));
            end
        end
    end
end


comprisons = [3 1];
comprisonsName = ["L_vs_S"];
for k = 1:1
    X = zeros(ndates,ndirections);
    Y = zeros(ndates,ndirections);
    for i = 1:ndates
        for j = 1:ndirections
            X(i,j) = rewardAxisMean(i,comprisons(k, 1),1,j) - rewardAxisMean(i,comprisons(k, 2),1,j);
            Y(i,j) = rewardAxisMean(i,comprisons(k, 1),2,j) - rewardAxisMean(i,comprisons(k, 2),2,j);
        end
    end
    rankTest2d(X(:), Y(:), "XLabel", "Reward axis range, Tiny targets (Spikes/s)", ...
                "YLabel", "Reward axis range, Huge targets (Spikes/s)", ...
                "OutputFolder", outputFolderStretching+rewardAxisName+"-stretch-"+comprisonsName(k))
end