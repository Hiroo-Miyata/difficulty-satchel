clear all;
%% Input 
dates = ["0216", "0405", "0414"]; %, "0405", "0414"
for day = 1:length(dates)
date = dates(day);
rootDir = "../";
target = "TO";
beforeT = 200;
afterT = 500;
conditionState = [70];
outputFile = rootDir + "data/processed/" + date + "_"+target+"_"+num2str(beforeT)...
    +"_"+num2str(afterT)+"_"+"70";

%% load
preprocessedFolder = rootDir + "data/preprocessed/";
preprocessedfile = dir(preprocessedFolder+"*"+date+"*.mat");
load(preprocessedfile.folder +"/"+ preprocessedfile.name);

%% data around GoCue -200~+600 ms
if target == "TO"
    state = 70;
elseif target == "GC"
    state = 11;
elseif target == "HT"
    state = 159;
end

% array indicate extracted data
validTrials = false(length(trialData), 1);
%% extract data corresponding to the inputs
for i = 1:length(trialData)
    stateTransition = trialData(i).stateTable;
    catchLabel = trialData(i).catchLabel;
    %% if the trial has the state of interest: [3]
    if all(ismember(conditionState, stateTransition(1,:))) == 1 && catchLabel == 0
        stateTime = stateTransition(2, find(stateTransition(1, :)==state));
        newFR = 1000 * trialData(i).firingRates(:, stateTime-beforeT:stateTime+afterT);
        trialData(i).firingRates = newFR;

        validTrials(i) = true;
    end
end

trialData = trialData(validTrials);
% trialData = rmfield(trialData,["firingRate"]);
save(outputFile, "trialData")

end