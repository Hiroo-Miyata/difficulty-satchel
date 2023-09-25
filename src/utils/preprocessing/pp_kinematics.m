clear all; close all;  addpath(genpath("../plotting/"));
dates = ["0417", "0419", "0420"];
ndates = length(dates);
% Global variables
dataDir = "../../../data/preprocessed/";
processedDataDir = "../../../data/preprocessed/";

% load the data
for day = 1:ndates
date = dates(day);
fileDir = dataDir + date + "/";
rawFile = dir(fileDir + "*_KIN_*.mat");
load(fileDir + rawFile.name);

fs = 1000;
% input
% kinematicData (* x 1 struct array)
% - position (ntime x 2 double array)

% output
% kinematicData (* x 1 struct array)
% - position (ntime x 2 double array)
% - velocity (ntime x 2 double array)
% - acceleration (ntime x 2 double array)

% calculate velocity and acceleration
initialFs = 60; % nominal frequency of the hand data
finalFs = 1000; % final desired frequency to interpolate to
fc = 10;        % cutoff frequency of LPF for kinematics
order = 4/2;    % order (/2) of LPF for kinematics; filtfilt means the order is applied twice, hence /2
[bLPF,aLPF] = butter(order,fc/(initialFs/2),'low'); % make the filter
for i = 1:length(kinematicData)
    % at first update kinematicData(i).position
    ntimes = size(kinematicData(i).position, 1);

     % Our sampling is sometimes uneven, which we can't filter on; but, the
    % data are pretty smooth as is, so we can confidently interpolate to a
    % valid sampling rate. We'll use 60Hz, the nominal resolution
    dt = 1000/initialFs;
    curKinTime = find(~ismember(kinematicData(i).position,[0,0],'rows'));
    if length(curKinTime) > 10
        posTime = (curKinTime(1):dt:curKinTime(end))';
        curPos = kinematicData(i).position(curKinTime,:);
        curPos_rs = interp1(curKinTime,curPos,posTime,'spline');
        
        % Filter the data and calculate velocity and acceleration
        curPos_filt = filtfilt(bLPF,aLPF,double(curPos_rs)); % in mm
        curVel_filt = diff(curPos_filt)./diff(posTime); % mm/ms = m/s
        velTime = (posTime(1:end-1)+posTime(2:end))/2; % this velocity approximates that between time indices
        curAcc_filt = diff(curVel_filt)./(diff(velTime)/1000); % (m/s)/(ms/1000) = (m/s)/s = m/s^2
        accTime = posTime(2:end-1); % same as doing the averages for the velocity time
        
        % Interpolate to the new fs and align to center target appearance
        time = (ceil(min(posTime)):1000/finalFs:floor(max(posTime)))';
        kinematicData(i).position(time(1):time(end), :) = interp1(posTime,curPos_filt,time,'spline');
        kinematicData(i).velocity(time(1):time(end), :) = interp1(velTime,curVel_filt,time,'spline');
        kinematicData(i).acceleration(time(1):time(end), :) = interp1(accTime,curAcc_filt,time,'spline');
    else
        kinematicData(i).position = nan;
        kinematicData(i).velocity = nan;
        kinematicData(i).acceleration = nan;
    end
end

% save the data
% change the filename to *_KIN2_*.mat
filename = rawFile.name;
filename = strrep(filename, "_KIN_", "_KIN2_");
save(fileDir + filename, 'kinematicData');

end