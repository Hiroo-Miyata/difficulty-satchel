close all; clear all; addpath(genpath("./utils/plotting/"));
dates = ["0417", "0419", "0420"];
ndates = length(dates);
folderName = "202309-summary/neuralData_GoCue";

Yall = cell(ndates, 2, 2, 8);
for d = 1:ndates
date = dates(d);
rootDir = "../";
load(rootDir+"data/processed/"+date+"_GC_200_100_11159.mat");
outputFolder = "../results/"+folderName+"/projection-2D/"; makeDir(outputFolder); 
outputFolderProjection = "../results/"+folderName+"/rewardAxis-projection/"; makeDir(outputFolderProjection);
outputFolderParams= "../results/"+folderName+"/rewardAxis-params/"; makeDir(outputFolderParams);
outputFolderStretching = "../results/"+folderName+"/stretching/"; makeDir(outputFolderStretching);
analysisBin = (50:250); % Pay Attention!! HT=(350:550), GC=(50:250), TO=(400:600)
axisOutputFileName = rootDir+"/interim/nonstitched_GC_rewardAxis_"+date+".mat";

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
meanNeuralData = squeeze(mean(neuralData_byParameters_mean,[1 3]));
% meanNeuralData = squeeze(meanNeuralData(:, 2, :));
[wR,zR,eigVls_R] = pca(meanNeuralData-mean(meanNeuralData, 1),'numcomponents',nrewards-1);
if zR(nrewards, 1) < zR(1, 1)
    wR(:, 1) = -wR(:, 1);
end

% reward axis of Tiny
meanNeuralData = squeeze(mean(neuralData_byParameters_mean,[1]));
meanNeuralData = squeeze(meanNeuralData(:, 1, :));
[wRTiny,zRTiny,eigVls_RTiny] = pca(meanNeuralData-mean(meanNeuralData, 1),'numcomponents',nrewards-1);
if zRTiny(nrewards, 1) < zRTiny(1, 1)
    wRTiny(:, 1) = -wRTiny(:, 1);
end

% reward axis of Huge
meanNeuralData = squeeze(mean(neuralData_byParameters_mean,[1]));
meanNeuralData = squeeze(meanNeuralData(:, 2, :));
[wRHuge,zRHuge,eigVls_RHuge] = pca(meanNeuralData-mean(meanNeuralData, 1),'numcomponents',nrewards-1);
if zRHuge(nrewards, 1) < zRHuge(1, 1)
    wRHuge(:, 1) = -wRHuge(:, 1);
end

% difficulty axis
meanNeuralData = squeeze(mean(neuralData_byParameters_mean,[1 2]));
[wDif,zDif,eigVls_Dif]  = pca(meanNeuralData-mean(meanNeuralData, 1),'numcomponents',ndifficulties-1);
if zDif(ndifficulties, 1) < zDif(1, 1)
    wDif(:, 1) = -wDif(:, 1);
end

% R-Diff axis
meanNeuralData = squeeze(mean(neuralData_byParameters_mean,[1]));
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

%% angle of wDif and wR
dot_product = dot(wRTiny(:,1), wRHuge(:,1));
angle = acos(dot_product) * 180 / pi;
if angle > 90
    angle = 180 - angle;
end
disp(angle);

%% plot reward PC1 and difficulty PC1 by 2D scatter plot
% neuralData_onDPC = (neuralData - mean(neuralData, 1)) * wRDiff(:,1:2);
% figure;
% hold on;
% for i=1:ndifficulties
%     for j=1:nrewards
%         curInds = rewardLabels==rewards(j) & difficultyLabels==difficulties(i);
%         X(i,j) = mean(neuralData_onDPC(curInds, 1));
%         Xerr(i,j) = std(neuralData_onDPC(curInds, 1))/sqrt(sum(curInds));
%         Y(i,j) = mean(neuralData_onDPC(curInds, 2));
%         Yerr(i,j) = std(neuralData_onDPC(curInds, 2))/sqrt(sum(curInds));
%     end
%     h(i) = errorbar(X(i,:), Y(i,:), Yerr(i,:), Yerr(i,:), Xerr(i,:), Xerr(i,:), "Color", diffColors(i,:), 'LineWidth', 2.5); hold on;
%     scatter(X(i,:), Y(i,:), 100, rewColors(rewards, :), "filled"); hold on;
% end
% set(gca, 'fontsize', 20, 'fontname', 'arial', 'tickdir', 'out', 'fontweight', 'bold');
% maxlim = ceil(max(X(:)) / 10)*10;
% minlim = floor(min(X(:)) / 10)*10;
% ylim([minlim maxlim]); xlim([minlim maxlim]);
% % ylim([-30 20]); xlim([-30 20]);
% xlabel("RDiff Axis: PC1 (Hz): "+round(eigVls_RDiff(1)/sum(eigVls_RDiff)*100, 2)+"%");
% ylabel("RDiff Axis: PC2 (Hz): "+round(eigVls_RDiff(2)/sum(eigVls_RDiff)*100, 2)+"%");
% legend(h, ["Tiny", "Huge"], Location="best"); set(gcf,'position',[0,0,550,550]);
% saveas(gcf, outputFolder+"2Dprojection-"+date+".jpg");


%% plot reward PC1 and difficulty PC1
% rAxes = ["Tiny", "Huge", "All"];
% for k=1:3 % 1:2
%     Y = zeros(3,2); % reward x difficulty
%     Yerr = zeros(3,2);
%     if k==1
%         neuralData_onDPC = (neuralData - mean(neuralData, 1)) * wRTiny(:,1);
%     elseif k==2
%         neuralData_onDPC = (neuralData - mean(neuralData, 1)) * wRHuge(:,1);
%     else
%         neuralData_onDPC = (neuralData - mean(neuralData, 1)) * wR(:,1);
%     end
% 
%     for i=1:nrewards
%         for j=1:ndifficulties
%             curInds = rewardLabels==rewards(i) & difficultyLabels==difficulties(j);
%             Y(i,j) = mean(neuralData_onDPC(curInds));
%             Yerr(i,j) = std(neuralData_onDPC(curInds))/sqrt(sum(curInds));
%         end
%     end
%     figure;
%     for j=1:ndifficulties
%         errorbar(1:nrewards, Y(:,j), Yerr(:,j), "Color", diffColors(j,:), 'LineWidth', 2); hold on;
%     end
%     set(gca, 'fontsize', 20, 'fontname', 'arial', 'tickdir', 'out', 'fontweight', 'bold');
%     xticks([1 2 3]); xticklabels(["S", "M", "L"]); xlim([0.7 3.3]); legend(["Tiny", "Huge"], Location="best")
%     set(gcf,'position',[100,100,400,650]);
%     ylabel("Reward Axis (Hz)");
%     saveas(gcf, outputFolderProjection+date+"-Axis-in-"+rAxes(k)+".jpg"); close all;
% end


close all;
neuralData_onDPC = (neuralData - mean(neuralData, 1)) * wR(:,1);
data = neuralData_onDPC';
[rewardNames, rewardLegends, rewColors, diffColors, direColors, DiffStyle, DelayTimes, nDelayTimes] = getExperimentConstants();

figure; hold on;
for j=1:nfocusBlocks
    Y = zeros(nrewards,1);
    Yerr = zeros(nrewards,1);
    for i=1:nrewards
        curInds = rewardLabels == rewards(i) & focusBlockLabels == focusBlocks(j) & ~isnan(data);
        Y(i) = mean(data(curInds));
        Yerr(i) = std(data(curInds)) / sqrt(sum(curInds));
    end
    errorbar(1:nrewards, Y,Yerr, Color="k", LineWidth=2, LineStyle=DiffStyle{j});
end
set(gca, 'fontsize', 20, 'fontname', 'arial', 'tickdir', 'out', 'fontweight', 'bold');
xlim([0.7, nrewards+0.3]); xticks(1:nrewards); xticklabels(rewardNames); ylabel("Reward Axis (Spike/s)"); legend(["Hard", "Easy"], Location="best");
set(gcf,'position',[0,0,550,550]);
saveas(gcf, outputFolderParams+date+"-vs-FocusBlock.jpg"); close all;

end
% rewardAxisRaw = Yall;
% rewardAxisMean = zeros(size(Yall));
% rewardAxisStd  = zeros(size(Yall));
% for d=1:size(Yall, 1)
%     for j=1:ndifficulties
%         for i=1:nrewards
%             for k=1:ndirections
%                 Ys = Yall{d,i,j,k};
%                 rewardAxisMean(d,i,j,k) = mean(Ys);
%                 rewardAxisStd(d,i,j,k)  = std(Ys) / sqrt(length(Ys));
%             end
%         end
%     end
% end
% 
% 
% comprisons = [2 1];
% comprisonsName = ["L_vs_S"];
% for k = 1
% 
%     X = zeros(ndates,ndirections);
%     Y = zeros(ndates,ndirections);
%     for i = 1:ndates
%         for j = 1:ndirections
%             X(i,j) = rewardAxisMean(i,comprisons(k, 1),1,j) - rewardAxisMean(i,comprisons(k, 2),1,j);
%             Y(i,j) = rewardAxisMean(i,comprisons(k, 1),2,j) - rewardAxisMean(i,comprisons(k, 2),2,j);
%         end
%     end
%     rankTest2d(X(:), Y(:), "XLabel", "Reward axis range, Tiny targets (Spikes/s)", ...
%                 "YLabel", "Reward axis range, Huge targets (Spikes/s)", ...
%                 "OutputFolder", outputFolderStretching+rewardAxisName+"-stretch-"+comprisonsName(k))
% end