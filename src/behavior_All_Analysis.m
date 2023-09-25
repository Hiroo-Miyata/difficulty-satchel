close all; clear all; addpath(genpath("./utils/plotting/"));
dates = ["0417", "0419", "0420"];
fileName = "202309-summary";
outputFolderSR = "../results/"+ fileName + "/behavior_success-rate/"; makeDir(outputFolderSR);
outputFolderVigor = "../results/"+ fileName + "/behavior_vigors/"; makeDir(outputFolderVigor);

%% conbine across days
trialDataAll = struct.empty;
behavDataAll = struct.empty;
trialNumBegin = 0;
for d = 1:length(dates)
    date = dates(d);
    dataFolder = "../data/preprocessed/"+date+"/";
    taskfile = dir(dataFolder+ "*" +date+ "*_TSK_*.mat");
    load(taskfile.folder +"/"+ taskfile.name);
    trialNum = [trialData.trial];
    for i=1:length(trialData)
        trialData(i).dayLabel = d;
        trialData(i).newTrial = trialData(i).trial + trialNumBegin;
    end

    trialDataAll = cat(2, trialDataAll, trialData); clear trialData
    trialNumBegin = trialNumBegin + max(trialNum);
end; clear taskfile date dataFolder
trialData = trialDataAll; clear trialDataAll

catchLabels = [trialData.catchLabel];
trialData = trialData(catchLabels==0);

%% get labels
directionLabels = [trialData.directionLabel]; directions = unique(directionLabels); ndirections = length(directions);
rewardLabels = [trialData.rewardLabel]; rewards = unique(rewardLabels); nrewards = length(rewards);
difficultyLabels =[trialData.targetSizeLabel]; difficulties = unique(difficultyLabels); ndifficulties = length(difficulties);
delayTimes = [trialData.delayLength]; delayTimes = round(delayTimes / 50) * 50;
dayLabels = [trialData.dayLabel];
trialNums = [trialData.newTrial];
rewColors = [1 0 0; 1 0.6470 0; 0 0 1]; diffColors = [0 0.447 0.741; 0.466 0.674 0.188];
direColors = {[1 .5 .5],[.75 .75 .5],[.5 1 .5],[.25 .75 .5],[0 .5 .5],[0.25 0.25 .5],[0.5 0 .5],[0.75 0.25 .5]};
DiffStyle = ["-", ":"];
DelayTimes = 500:50:2000; nDelayTimes = length(DelayTimes);

%% get success rate and movement time
ntrials = length(trialData);
stateLabels = nan(1, ntrials);
movementTimes = nan(1, ntrials);
for i=1:ntrials
    stateTransition = trialData(i).stateTable;
    if all(ismember([70], stateTransition(1,:))) == 1
        if all(ismember([11], stateTransition(1,:))) == 1
            if all(ismember([159], stateTransition(1,:))) == 1
                if all(ismember([150], stateTransition(1,:))) == 1
                    stateLabel = 1;
                else
                    stateLabel = -13;
                end
                stateTime = stateTransition(2, find(stateTransition(1, :)==11));
                endTime = stateTransition(2, find(stateTransition(1, :)==159));
                movementTimes(i) = endTime - stateTime;
            else
                stateLabel = -12;
            end
        else
            stateLabel = -11;
        end
    else
        stateLabel = -10;
    end
    stateLabels(i) = stateLabel;
end



%% SUCCESS RATE
failedLabels = stateLabels ~= 1;
VS_reward_size_success_rate(failedLabels, rewardLabels, difficultyLabels, "Label", "success rate (%)", "OutputFolder", outputFolderSR+"SR")
% VS_delay_time_success_rate(failedLabels, rewardLabels, difficultyLabels, delayTimes, "Label", "success rate (%)", "OutputFolder", outputFolderSR+"SR")
VS_direction_success_rate(failedLabels, rewardLabels, difficultyLabels, directionLabels, "Label", "success rate (%)", "OutputFolder", outputFolderSR+"SR")
% VS_trial(trialNums, 100 * int16(~failedLabels), rewardLabels, difficultyLabels, "Label", "success rate (%)", "OutputFolder", outputFolderSR+"SR", ...
%         "IsMultipleDays",true, "DayLabels", dayLabels, "FigWidth", 1600)

%% success in preparation rate
delayFailure = stateLabels == -11;
VS_reward_size_success_rate(delayFailure, rewardLabels, difficultyLabels, "Label", "success in preparation rate (%)", "OutputFolder", outputFolderSR+"SRInPreparation", "Ylim", [60 100])
% VS_delay_time_success_rate(delayFailure, rewardLabels, difficultyLabels, delayTimes, "Label", "success in preparation rate (%)", "OutputFolder", outputFolderSR+"SRInPreparation", "Ylim", [60 100])
VS_direction_success_rate(delayFailure, rewardLabels, difficultyLabels, directionLabels, "Label", "success in preparation rate (%)", "OutputFolder", outputFolderSR+"SRInPreparation", "Ylim", [60 100])
% VS_trial(trialNums, 100 * int16(~delayFailure), rewardLabels, difficultyLabels, "Label", "success in preparation rate (%)", "OutputFolder", outputFolderSR+"SRInPreparation", ...
%         "IsMultipleDays",true, "DayLabels", dayLabels, "FigWidth", 1600)

%% success in reaching rate
reachFailure = stateLabels == -12;
VS_reward_size_success_rate(reachFailure, rewardLabels, difficultyLabels, "Label", "success in reaching rate (%)", "OutputFolder", outputFolderSR+"SRInReaching", "Ylim", [60 100])
% VS_delay_time_success_rate(reachFailure, rewardLabels, difficultyLabels, delayTimes, "Label", "success in reaching rate (%)", "OutputFolder", outputFolderSR+"SRInReaching", "Ylim", [60 100])
VS_direction_success_rate(reachFailure, rewardLabels, difficultyLabels, directionLabels, "Label", "success in reaching rate (%)", "OutputFolder", outputFolderSR+"SRInReaching", "Ylim", [60 100])
% VS_trial(trialNums, 100 * int16(~reachFailure), rewardLabels, difficultyLabels, "Label", "success in reaching rate (%)", "OutputFolder", outputFolderSR+"SRInReaching", ...
%         "IsMultipleDays",true, "DayLabels", dayLabels, "FigWidth", 1600)

%% success in holding rate
holdFailure = stateLabels == -13;
VS_reward_size_success_rate(holdFailure, rewardLabels, difficultyLabels, "Label", "success in holding rate (%)", "OutputFolder", outputFolderSR+"SRInHolding", "Ylim", [60 100])
% VS_delay_time_success_rate(holdFailure, rewardLabels, difficultyLabels, delayTimes, "Label", "success in holding rate (%)", "OutputFolder", outputFolderSR+"SRInHolding", "Ylim", [60 100])
VS_direction_success_rate(holdFailure, rewardLabels, difficultyLabels, directionLabels, "Label", "success in holding rate (%)", "OutputFolder", outputFolderSR+"SRInHolding", "Ylim", [60 100])
% VS_trial(trialNums, 100 * int16(~holdFailure), rewardLabels, difficultyLabels, "Label", "success in holding rate (%)", "OutputFolder", outputFolderSR+"SRInHolding", ...
%         "IsMultipleDays",true, "DayLabels", dayLabels, "FigWidth", 1600)

%% REACTION TIME
VS_reward_size(movementTimes, rewardLabels, difficultyLabels, "Label", "Movement Time (ms)", "OutputFolder", outputFolderVigor+"MovementTime")
% VS_delay_time(reactionTimes, rewardLabels, difficultyLabels, delayTimes, "Label", "Reaction Time (ms)", "OutputFolder", outputFolderVigor+"MovementTime")
VS_direction(movementTimes, rewardLabels, difficultyLabels, directionLabels, "Label", "Movement Time (ms)", "OutputFolder", outputFolderVigor+"MovementTime")
% VS_trial(trialNums, reactionTimes, rewardLabels, difficultyLabels, "Label", "Reaction Time (ms)", "OutputFolder", outputFolderVigor+"MovementTime", ...
%         "IsMultipleDays",true, "DayLabels", dayLabels, "FigWidth", 1600)
