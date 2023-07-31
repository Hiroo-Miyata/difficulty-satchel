clear all; close all;  addpath(genpath("../plotting/"));
dates = ["0216", "0405", "0414", "0417", "0418", "0419", "0420", "0425", "0426", "0428"];
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

for i = 1:length(kinematicData)
    % at first update kinematicData(i).position
    ntimes = size(kinematicData(i).position, 1);
    kinematicData(i).velocity = zeros(ntimes, 2);
    kinematicData(i).acceleration = zeros(ntimes, 2);

    % dataIdx has the data. other index is filled by zero. so ignore it
    dataIdx = ~ismember(kinematicData(i).position, [0, 0], 'rows');
    dataIdx = find(dataIdx(:, 1) == 1);
    ndata = length(dataIdx);
    for j = 2:ndata
        timeLength = dataIdx(j) - dataIdx(j - 1);
        position = kinematicData(i).position(dataIdx(j), :);
        kinematicData(i).position(dataIdx(j-1)+1:dataIdx(j), :) = ones(timeLength, 1) * position;
        prevPosition = kinematicData(i).position(dataIdx(j - 1), :);
        velocity = (position - prevPosition) / (dataIdx(j) - dataIdx(j - 1)) * fs;
        kinematicData(i).velocity(dataIdx(j-1)+1:dataIdx(j), :) = ones(timeLength, 1) * velocity;
        if j > 2
            prevVelocity = kinematicData(i).velocity(dataIdx(j - 1), :);
            acceleration = (velocity - prevVelocity) / (dataIdx(j) - dataIdx(j - 1)) * fs;
            kinematicData(i).acceleration(dataIdx(j-1)+1:dataIdx(j), :) = ones(timeLength, 1) * acceleration;
        end
    end
    
    for d = 1:2
        kinematicData(i).velocity(:, d) = movmean(kinematicData(i).velocity(:, d), 100);
        kinematicData(i).acceleration(:, d) = movmean(kinematicData(i).acceleration(:, d), 100);
    end
end

% save the data
% change the filename to *_KIN2_*.mat
filename = rawFile.name;
filename = strrep(filename, "_KIN_", "_KIN2_");
save(fileDir + filename, 'kinematicData');

end