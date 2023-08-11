clear all;
%% Input 
dates = ["0216", "0405", "0414", "0417", "0418", "0419", "0420", "0425", "0426", "0428"];
for day = 1:length(dates)
date = dates(day);
rootDir = "../";
target = "HT";
beforeT = 200;
afterT = 375;
conditionState = [11 159];
outputFile = rootDir + "data/processed/" + date + "_"+target+"_"+num2str(beforeT)...
    +"_"+num2str(afterT)+"_"+"11159";

%% load
preprocessedFolder = rootDir + "data/preprocessed/"+date + "/";
synchfile = dir(preprocessedFolder+'/*_TSK_*.mat');
neurofile = dir(preprocessedFolder+'/*_NER_*.mat');
load(synchfile.folder +"/"+ synchfile.name);
load(neurofile.folder +"/"+ neurofile.name);

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
    if all(ismember(conditionState, stateTransition(1,:))) == 1 && catchLabel == 0 && spikeInfo.badTrials(i) == 0
        stateTime = stateTransition(2, find(stateTransition(1, :)==state));
        
        firingRates = zeros(size(neuralData(i).spikeMatrix)); % nneurons * times
        firingRates = 1000 * neuralData(i).spikeMatrix;
        trialData(i).firingRates = firingRates(:, stateTime-beforeT:stateTime+afterT);

        validTrials(i) = true;
    end
end

trialData = trialData(validTrials);
% trialData = rmfield(trialData,["firingRate"]);
save(outputFile, "trialData")

end