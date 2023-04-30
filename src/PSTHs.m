close all; clear all;
dates = ["0216", "0405", "0414"];

for day = 1:length(dates)
date = dates(day);
rootDir = "../";
load(rootDir+"data/processed/"+date+"_GC_300_100_11_159.mat");
% load(rootDir+"data/processed/"+date+"_TO_200_500_70.mat");
outputFolder = "../results/202304w2/PSTHs-GC/"+date;

% remove calibration part 
trialNums = [trialData.trial];
trialData = trialData(trialNums > 24);
% Get labels
directionLabels = [trialData.directionLabel]; directions = unique(directionLabels); ndirections = length(directions);
rewardLabels = [trialData.rewardLabel]; rewards = unique(rewardLabels); nrewards = length(rewards);
difficultyLabels =[trialData.targetSizeLabel]; difficulties = unique(difficultyLabels); ndifficulties = length(difficulties);
rewColors = [1 0 0; 1 0.6470 0; 0 0 1];
diffColors  = [0 0.447 0.741; 0.466 0.674 0.188]; %tiny huge: blue and green
direColors = {[1 .5 .5],[.75 .75 .5],[.5 1 .5],[.25 .75 .5],[0 .5 .5],[0.25 0.25 .5],[0.5 0 .5],[0.75 0.25 .5]};
trialNums =  [trialData.trial];

% Get data and avg within dir x rew
analysisBin = (50:350);
neuralActivity = cat(3, trialData(:).firingRates);
[nneurons, ntimeBin ,ntrials] = size(neuralActivity);

% save PSTHs of all neurons
neuralData_byParameters_mean = nan(ndirections,nrewards,ndifficulties,nneurons, ntimeBin);
for l = 1:ndirections
    for r = 1:nrewards
        for j = 1:ndifficulties
            curInds = directionLabels==directions(l) & rewardLabels==rewards(r) & difficultyLabels==difficulties(j);
            neuralData_byParameters_mean(l,r,j,:, :) = mean(neuralActivity(:, :, curInds), 3);
        end; clear j
    end; clear r
end; clear l

meanPSTHs = squeeze(mean(neuralData_byParameters_mean,[1 3])); % nreward x nneuron x ntime

% plot PSTHs of all reward condition in each neuron.

for n = 1:nneurons
    figure; hold on
    for r = 1:nrewards
        Y = squeeze(meanPSTHs(r,n,:));
        sigma = 25;
        kernel = normpdf(-3*sigma:3*sigma,0,sigma);
        Y = conv(Y, kernel, "same");
        Bin = (50:350);
        plot(Bin, Y(Bin), 'color', rewColors(r,:), 'linewidth', 2)
    end; clear r
    legend(["Small", "Medium", "Large"])
    xlabel("Time (ms)")
    ylabel("Firing Rate (Hz)")
    xticks([50 300 350])
    xticklabels(["-250", "GC", "+50"])
    set(gca, 'fontsize', 14, 'fontname', 'arial', 'tickdir', 'out','box','off');
    saveas(gcf, outputFolder+"-neuron"+n+".png")
    close all;
end; clear n

% find the neurons which modulate their firing rate as a function of reward
% At first, get the firing rate of each neuron during the analysis window in each trial
% in each reward condition, calculate the mean and std of firing rate of each neuron
% Then , in each neuron, find the neuron which modulate their firing rate as a function of reward
% To do this, we compare the mean firing rate in Small and Medium reward condition at first, then
% compare the mean firing rate in Medium and Large reward condition.
% if it show monotonic encoding, classify the neuron based on the criteria below
% type1 : significant diff b/w Small and Medium, and significant diff b/w Medium and Large
% type2 : significant diff b/w Small and Medium, but no significant diff b/w Medium and Large
% type3 : no significant diff b/w Small and Medium, but significant diff b/w Medium and Large
% type4 : no significant diff b/w Small and Medium, and no significant diff b/w Medium and Large

FR_mean = zeros(nneurons, nrewards);
FR_std = zeros(nneurons, nrewards);
FRs = squeeze(mean(neuralActivity(:, analysisBin, :), 2)); % nneuron x ntrial
for r = 1:nrewards
    FR_mean(:,r) = mean(FRs(:,rewardLabels==rewards(r)),2);
    FR_std(:,r) = std(FRs(:,rewardLabels==rewards(r)),0,2);
end; clear r

rewardNeuronsCount = zeros(4,1);
rewardNeurons = cell(4,1);

for n = 1:nneurons
    dif_SM = FR_mean(n,1) - FR_mean(n,2);
    dif_ML = FR_mean(n,2) - FR_mean(n,3);
    if (dif_SM > 0 && dif_ML > 0) || (dif_SM < 0 && dif_ML < 0)
        % get p-value the data size is different in each reward condition
        [~,p_SM] = ttest2(FRs(n,rewardLabels==rewards(1)), FRs(n,rewardLabels==rewards(2)));
        [~,p_ML] = ttest2(FRs(n,rewardLabels==rewards(2)), FRs(n,rewardLabels==rewards(3)));
        if p_SM < 0.05 && p_ML < 0.05
            rewardNeuronsCount(1) = rewardNeuronsCount(1) + 1;
            rewardNeurons{1} = [rewardNeurons{1} n];
        elseif p_SM < 0.05 && p_ML >= 0.05
            rewardNeuronsCount(2) = rewardNeuronsCount(2) + 1;
            rewardNeurons{2} = [rewardNeurons{2} n];
        elseif p_SM >= 0.05 && p_ML < 0.05
            rewardNeuronsCount(3) = rewardNeuronsCount(3) + 1;
            rewardNeurons{3} = [rewardNeurons{3} n];
        else
            rewardNeuronsCount(4) = rewardNeuronsCount(4) + 1;
            rewardNeurons{4} = [rewardNeurons{4} n];
        end
    end
end; clear n
save(outputFolder+"-rewardNeurons.mat", "rewardNeurons", "rewardNeuronsCount");

disp(rewardNeuronsCount)
disp(nneurons)

end