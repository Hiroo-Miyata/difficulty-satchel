clear all;
%% Input 
dates = ["0417", "0419", "0420"];
for day = 1:length(dates)
date = dates(day);
rootDir = "../";
target = "DF";
beforeT = 200;
afterT = 1000;
conditionState = [11 154];
outputFile = rootDir + "data/processed/non_stitched/" + date + "_"+target+"_"+num2str(beforeT)...
    +"_"+num2str(afterT)+"_"+"11154";

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
elseif target == "DF"
    state = 154;
elseif target == "RF"
    state = 153;
end

% array indicate extracted data
validTrials = false(length(trialData), 1);
%% extract data corresponding to the inputs
for i = 1:length(trialData)
    stateTransition = trialData(i).stateTable;
    catchLabel = trialData(i).catchLabel;
    trialNum = trialData(i).trial;
    %% if the trial has the state of interest: [3]
    if all(ismember(conditionState, stateTransition(1,:))) == 1 && catchLabel == 0 && spikeInfo.badTrials(i) == 0
        stateTime = stateTransition(2, find(stateTransition(1, :)==state));
        lastT = max(int16(trialData(i).time * 1000));
        
        if (stateTime + afterT) < lastT
            trialData(i).firingRates = 1000 * neuralData(i).spikeMatrix(:, stateTime-beforeT:stateTime+afterT);          
            validTrials(i) = true;
        else
            if trialData(i+1).trial == (trialNum + 1)

                firingRate = neuralData(i).spikeMatrix;
                [nneurons, ntimeBins] = size(firingRate);
                combinedFR = nan(nneurons, afterT+beforeT+1);
                dataLimit = lastT-stateTime;
                combinedFR(:, 1:(dataLimit+beforeT+1)) = 1000 * firingRate(:, stateTime-beforeT:end);
                firingRate = neuralData(i+1).spikeMatrix;
                combinedFR(:, (dataLimit+beforeT+2):end) = 1000 * firingRate(:, 1:afterT-dataLimit);
                trialData(i).firingRates = combinedFR;

                validTrials(i) = true;
            end

        end
      
    end
end

trialData = trialData(validTrials);
% trialData = rmfield(trialData,["firingRate"]);
save(outputFile, "trialData")

end