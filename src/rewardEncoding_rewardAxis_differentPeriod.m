close all; clear all; addpath(genpath("./utils/plotting/"));
dates = ["0417", "0419", "0420"];
ndates = length(dates);
folderName = "202307w1-all/neuralData_TargetHold";

for d = 1:ndates
date = dates(d);
rootDir = "../";
load(rootDir+"data/processed/"+date+"_HT_200_375_11159.mat");
outputFolder = "../results/"+folderName+"/projection-2D_GCrewardAxis/"; makeDir(outputFolder); 
outputFolderProjection = "../results/"+folderName+"/rewardAxis-projection_GCrewardAxis/"; makeDir(outputFolderProjection);
outputFolderStretching = "../results/"+folderName+"/stretching_GCrewardAxis/"; makeDir(outputFolderStretching);
analysisBin = (350:550); % Pay Attention!! HT=(350:550), GC=(50:250), TO=(400:600)
axisOutputFileName = rootDir+"/interim/nonstitched_GC_rewardAxis_"+date+".mat";

% remove calibration part 
trialNums = [trialData.trial];
trialData = trialData(trialNums > 24);
% Get labels
directionLabels = [trialData.directionLabel]; directions = unique(directionLabels); ndirections = length(directions);
rewardLabels = [trialData.rewardLabel]; rewards = unique(rewardLabels); nrewards = length(rewards);
difficultyLabels =[trialData.targetSizeLabel]; difficulties = unique(difficultyLabels); ndifficulties = length(difficulties);
focusDifficultyLabels = [trialData.focusDifficultyLabel]; focusDifficulties = unique(focusDifficultyLabels); nfocusDifficulties = length(focusDifficulties);
rewColors = [1 0 0; 1 0.6470 0; 0 0 1];
diffColors  = [0 0.447 0.741; 0.466 0.674 0.188]; %tiny huge: blue and green
direColors = {[1 .5 .5],[.75 .75 .5],[.5 1 .5],[.25 .75 .5],[0 .5 .5],[0.25 0.25 .5],[0.5 0 .5],[0.75 0.25 .5]};
trialNums =  [trialData.trial]; 

% Get data and avg within dir x rew
neuralActivity = cat(3, trialData(:).firingRates);
neuralData = mean(neuralActivity(:, analysisBin, :), 2);
neuralData = squeeze(neuralData)';
[ntrials,nneurons] = size(neuralData);

load(axisOutputFileName)

%% plot reward PC1 and difficulty PC1 by 2D scatter plot
neuralData_onDPC = (neuralData - mean(neuralData, 1)) * wRDiff(:,1:2);
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
    scatter(X(i,:), Y(i,:), 100, rewColors(rewards, :), "filled"); hold on;
end
set(gca, 'fontsize', 20, 'fontname', 'arial', 'tickdir', 'out', 'fontweight', 'bold');
maxlim = ceil(max(X(:)) / 10)*10;
minlim = floor(min(X(:)) / 10)*10;
ylim([minlim maxlim]); xlim([minlim maxlim]);
% ylim([-30 20]); xlim([-30 20]);
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

    for i=1:nrewards
        for j=1:ndifficulties
            curInds = rewardLabels==rewards(i) & difficultyLabels==difficulties(j);
            Y(i,j) = mean(neuralData_onDPC(curInds));
            Yerr(i,j) = std(neuralData_onDPC(curInds))/sqrt(sum(curInds));
        end
    end
    figure;
    for j=1:ndifficulties
        errorbar(1:nrewards, Y(:,j), Yerr(:,j), "Color", diffColors(j,:), 'LineWidth', 2); hold on;
    end
    set(gca, 'fontsize', 20, 'fontname', 'arial', 'tickdir', 'out', 'fontweight', 'bold');
    xticks([1 2 3]); xticklabels(["S", "M", "L"]); xlim([0.7 3.3]); legend(["Tiny", "Huge"], Location="best")
    set(gcf,'position',[100,100,400,650]);
    ylabel("Reward Axis (Hz)");
    saveas(gcf, outputFolderProjection+date+"-Axis-in-"+rAxes(k)+".jpg"); close all;
end

close all;
end
