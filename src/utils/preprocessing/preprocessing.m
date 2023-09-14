
clear all; close all;  addpath(genpath("../plotting/"));
dates = ["0417", "0419", "0420"];
ndates = length(dates);
% Global variables
dataDir = "../../../data/raw/";
processedDataDir = "../../../data/preprocessed/";

% load the data
for day = 1:ndates
date = dates(day);
rawFile = dir(dataDir + "*" + date + "*.mat");
load(dataDir + rawFile.name); 
outputFolder = processedDataDir + date + "/"; makeDir(outputFolder);
datOutM1 = dat; clear dat;

% the data structure
% datOutM1: 1 * ntrials : struct array with the fields:
%   block: 1 * 1 : double
%   channels: nchannels * 2 : int : the first column is the channel
%   number and the second column is unit number( there are only 1 or 255)
%   time: 1 * 2 : double : the first element is the start time and the
%   second element is the end time 
%   trialcodes: nevents * 3 : int : the first column is always 0 (ignore)
%   the second column is the event code and the third column is the time
%   there are event code and hand kinematics 142 in the second column, 
%   that is the indicator that the next two rows are the joystick position
%   row+1 has the x position, and row+2 has the y position.
%   event:  nevents * 3 : int : the first column is always 0 (ignore)
%   the second column is the event code and the third column is the timestamp
%   of sampling rate. it stores all events including 3XX events.
%   firstspike: double : the first spike time (time stamp of sampling rate)
%   spiketimesdiff: nspikes-1 * 1 : int : the spike times difference
%   (time stamp of sampling rate)
%   spikeinfo: nspikes * 2 : int : the idx is corresponding to the spiketimesdiff 
%   and the first column is the channel number and the second column is unit number(1 or 255)
%   results: 1 * 1 : int : sometimes it has NaN
%   params: 1 * 1 : struct array with the fields:
%     block: 1 * 1 : struct array with the fields:
%     trial: 1 * 1 : struct array with the fields:
%       focusDifficulty: 1 * 1 : double
%       targetAngle: 1 * 1 : double
%       variableRewardIdx: 1 * 1 : double
%       heJoystickMinHoldAngle: 1 * 1 : double
%       currBlock: 1 * 1 : double
%       targRad: 1 * 1 : double


%% convert the data to the format
% the data structure
% trialData: ntrials * 1 : struct array with the fields:
%   trial: 1 * 1 : int : the trial number
%   rewardLabel : 1 * 1 : double : the reward amount
%   directionLabel : 1 * 1 : double : the target direction
%   targetSizeLabel : 1 * 1 : double : the target size
%   focusDifficultyLabel : 1 * 1 : double : the difficulty of the task
%   catchLabel : 1 * 1 : int : whether it is a catch trial
%   startTime : 1 * 1 : double : the start time of the trial
%   endTime : 1 * 1 : double : the end time of the trial
%   stateTable: 2 * nevents : int : the first row is the event code and the second row is the time
%   firingRate : nchannels * ntimes : int : the firing rate of each channel (sampled at 1000Hz)
%   time : 1 * ntimes : double : the time of each sample (sampled at 1000Hz) start from 0
%   (time stamp of new sampleing rate: 1000Hz)
%   handKinematics: struct array with the fields:
%     position: 2 * ntimes : double : the first row is the x position and the second row is the y position
%     velocity: 2 * ntimes : double : the first row is the x velocity and the second row is the y velocity

fs = 30000;
newFs = 1000;

% get each labels
ntrials = length(datOutM1);
rewardLabels = zeros(ntrials,1);
directionLabels = zeros(ntrials,1);
targetSizeLabels = zeros(ntrials,1);
focusDifficultyLabels = zeros(ntrials,1);
for i = 1:ntrials
    rewardLabels(i,:) = datOutM1(i).params.trial.variableRewardIdx;
    directionLabels(i,:) = datOutM1(i).params.trial.targetAngle;
    targetSizeLabels(i,:) = datOutM1(i).params.trial.targRad;
    focusDifficultyLabels(i,:) = datOutM1(i).params.trial.focusDifficulty;
end

taskInfo = datOutM1(1).params;
taskInfo.rewards = unique(rewardLabels);
taskInfo.directions = unique(directionLabels);
taskInfo.targetSizes = unique(targetSizeLabels);
taskInfo.focusDifficulties = unique(focusDifficultyLabels);

channelLabels = datOutM1(1).channels(:,1);
taskInfo.channels = unique(channelLabels);

%%%%%
%for me
stateDict = [];
stateTransit = ["0-0"];
stateTransitLength = cell(1,1);
%%%%%


trialData = struct.empty(ntrials,0);
neuralData = struct.empty(ntrials,0);
kinematicData = struct.empty(ntrials,0);
XPOSData = [];
YPOSData = [];

spikeInfo.badTrials = false(ntrials, 1);
for i = 1:ntrials
    reward = datOutM1(i).params.trial.variableRewardIdx;
    rewardLabel = find(taskInfo.rewards == reward);
    direction = datOutM1(i).params.trial.targetAngle;
    directionLabel = find(taskInfo.directions == direction);
    targetSize = datOutM1(i).params.trial.targRad;
    targetSizeLabel = find(taskInfo.targetSizes == targetSize);
    focusDifficulty = datOutM1(i).params.trial.focusDifficulty;
    focusDifficultyLabel = find(taskInfo.focusDifficulties == focusDifficulty);
    block = datOutM1(i).params.trial.currBlock;

    % make a matrix of firing rate
    % put the spikes based on the spikesTimeDiff and spikeinfo

    ntime = ceil((datOutM1(i).time(2)-datOutM1(i).time(1))*newFs);
    nchannel = length(taskInfo.channels);
    firingRate = logical(zeros(nchannel,ntime));
    startTimeFs = ceil(datOutM1(i).time(1)*fs);

    % make a spike time array based on the spikeinfo
    spikeTimes = zeros(size(datOutM1(i).spikeinfo,1),1);
    spikeTimes(1) = datOutM1(i).firstspike + 1;
    for j = 2:size(datOutM1(i).spikeinfo,1)
        spikeTimes(j) = spikeTimes(j-1) + double(datOutM1(i).spiketimesdiff(j-1));
    end

    for j = 1:size(datOutM1(i).spikeinfo,1)
        spikeTime = ceil((spikeTimes(j) - startTimeFs)/fs*newFs);
        % ignore the spike if the spikeinfo(j, 2) is 255 (noise)
        if datOutM1(i).spikeinfo(j,2) == 255
            continue;
        end

        channel = datOutM1(i).spikeinfo(j,1);
        channelLabel = find(taskInfo.channels == channel);
        if spikeTime > 0 && spikeTime <= ntime
            % if there is a spike in the same time, it will be ignored but output a warning
            if firingRate(channelLabel,spikeTime) == 1
                fprintf("there is a spike in the same time. trial: " + i + ", channel: " + channel + ", spike time: " + spikeTime + "\n");
            end
            firingRate(channelLabel,spikeTime) = 1;
        else
            if spikeInfo.badTrials(i) == false
                fprintf("spike time is out of range. trial: " + i + ", channel: " + channel + ", spike time: " + spikeTime + "\n");
                spikeInfo.badTrials(i) = true;
            end
        end
    end

    startTime = datOutM1(i).time(1);
    endTime = datOutM1(i).time(2);
    trialCode = datOutM1(i).trialcodes(:, 2:3);

    % extract state table and hand kinematics from the trialCode
    % handposition indicater is 142

    stateTable = zeros(2,0);
    handKinematics.position = zeros(ntime, 2);
    catchLabel = 0;
    delayLength = NaN;
    j = 1;
    while j < size(trialCode,1)
        if trialCode(j, 1) == 142
            % extract hand kinematics
            xPos = 10000 - trialCode(j+1, 1);
            yPos = 10000 - trialCode(j+2, 1);
            % get the time
            timeIdx = ceil((trialCode(j, 2) - startTime)*newFs);
            handKinematics.position(timeIdx, 1) = xPos;
            handKinematics.position(timeIdx, 2) = yPos;

            XPOSData = cat(1, XPOSData, xPos);
            YPOSData = cat(1, YPOSData, yPos);
            
            j = j + 3;
        elseif trialCode(j,1) == 7777
            % which means huge target trial
            % targetSizeLabel = 1;
            stateTable(1, end+1) = trialCode(j, 1);
            timeIdx = ceil((trialCode(j, 2) - startTime)*newFs);
            stateTable(2, end) = timeIdx;
            j = j + 1;
        elseif trialCode(j,1) == 5555
            % which means catch trial
            catchLabel = 1;
            stateTable(1, end+1) = trialCode(j, 1);
            timeIdx = ceil((trialCode(j, 2) - startTime)*newFs);
            stateTable(2, end) = timeIdx;
            j = j + 1;
        elseif trialCode(j,1) > 2000 && trialCode(j,1) < 4001
            % which means delay length
            delayLength = trialCode(j,1) - 2000;
            j = j + 1;
        else
            % extract state table
            stateTable(1, end+1) = trialCode(j, 1);
            timeIdx = ceil((trialCode(j, 2) - startTime)*newFs);
            stateTable(2, end) = timeIdx;

            %%%%%
            % for me 
            if ~ismember(trialCode(j, 1), stateDict)
                stateDict = [stateDict, trialCode(j, 1)];
            end
            %%%%%

            j = j + 1;
        end
    end

%     if targetSizeLabel == 1 & all(ismember([11 159], trialCode(:, 1))) == 1
%         disp("this huge target trial have 11");
%     end

    trialData(i).trial = i;
    trialData(i).block = block;
    trialData(i).rewardLabel = rewardLabel;
    trialData(i).directionLabel = directionLabel;
    trialData(i).targetSizeLabel = targetSizeLabel;
    trialData(i).focusDifficultyLabel = focusDifficultyLabel;
    trialData(i).catchLabel = catchLabel;
    trialData(i).delayLength = delayLength;
    trialData(i).startTime = startTime;
    trialData(i).endTime = endTime;
    trialData(i).stateTable = stateTable;
    trialData(i).time = linspace(0, ntime*1/newFs, ntime);
    kinematicData(i).position = handKinematics.position;

    neuralData(i).spikeMatrix = firingRate;


    %%% for me
    for st = 2:size(stateTable, 2)
        transition = num2str(stateTable(1, st-1)) + "-" + num2str(stateTable(1, st));
        if ~ismember(transition, stateTransit)
            stateTransit = cat(1, stateTransit, transition);
            stateTransitLength{end+1} = [stateTable(2, st) - stateTable(2, st-1)];
        else
            idx = find(strcmp(stateTransit, transition));
            stateTransitLength{idx} = cat(1, stateTransitLength{idx}, stateTable(2, st) - stateTable(2, st-1));
        end
    end

end

stateDict = sort(stateDict);
% for i = 2:length(stateTransit)
%     % print the state transition name and mean + std of the length
%     disp(stateTransit(i) + ": " + mean(stateTransitLength{i}) + " +/- " + std(stateTransitLength{i}));
% end

% save the trialData
% the name of the file is the first two components (separated by '_') of the file name
save(outputFolder+rawFile.name + "_TSK_" + ".mat", 'trialData', 'taskInfo');
save(outputFolder+rawFile.name + "_NER_" + ".mat", 'neuralData', 'spikeInfo');
save(outputFolder+rawFile.name + "_KIN_" + ".mat", 'kinematicData');

disp("finished one day preprocessing")
figure; hold on;
scatter(XPOSData, YPOSData, 1, 'filled');
saveas(gcf, "results/handPosition_" + date + ".png"); close all;
end
