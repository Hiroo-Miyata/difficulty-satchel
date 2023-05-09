close all; clear all;
dates = ["0216", "0405", "0414"]; ndate = length(dates);
Yall = cell(length(dates), 3, 2, 8);

for day = 1:length(dates)
date = dates(day);
rootDir = "../";
load(rootDir+"data/processed/"+date+"_GC_300_100_11_159.mat");
% load(rootDir+"data/processed/"+date+"_TO_200_500_70.mat");
outputFolder = "../results/202305w2/rewardAxis-GC-2D/";

% remove calibration part 
% trialNums = [trialData.trial];
% trialData = trialData(trialNums > 24);
focusDifficultyLabels = [trialData.focusDifficultyLabel]; focusDifficulties = unique(focusDifficultyLabels); nfocusDifficulties = length(focusDifficulties);
trialData = trialData(focusDifficultyLabels==focusDifficulties(3) | focusDifficultyLabels==focusDifficulties(4));


% Get labels
directionLabels = [trialData.directionLabel]; directions = unique(directionLabels); ndirections = length(directions);
rewardLabels = [trialData.rewardLabel]; rewards = unique(rewardLabels); nrewards = length(rewards);
difficultyLabels =[trialData.targetSizeLabel]; difficulties = unique(difficultyLabels); ndifficulties = length(difficulties);
rewColors = [1 0 0; 1 0.6470 0; 0 0 1];
diffColors  = [0 0.447 0.741; 0.466 0.674 0.188]; %tiny huge: blue and green
direColors = {[1 .5 .5],[.75 .75 .5],[.5 1 .5],[.25 .75 .5],[0 .5 .5],[0.25 0.25 .5],[0.5 0 .5],[0.75 0.25 .5]};
trialNums =  [trialData.trial];

% Get data and avg within dir x rew
analysisBin = (150:350);
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

% Look at reward pts similarly
% meanNeuralData = squeeze(mean(neuralData_byParameters_mean,[1 3]));
% % meanNeuralData = squeeze(meanNeuralData(:, 2, :));
% [wR,zR,eigVls_R] = pca(meanNeuralData-mean(meanNeuralData, 1),'numcomponents',nrewards-1);
% % check the score of high reward condition has the largest value on the first PC respectively
% % if not, flip the eigen vector of the first two principal components
% if zR(3, 1) < zR(1, 1)
%     wR(:, 1) = -wR(:, 1);
% end
% % % Look at reward pts similarly
% meanNeuralData = squeeze(mean(neuralData_byParameters_mean,[1 2]));
% [wDif,zDif,eigVls_Dif]  = pca(meanNeuralData-mean(meanNeuralData, 1),'numcomponents',ndifficulties-1);
% % check the score of high reward condition has the largest value on the first PC respectively
% % if not, flip the eigen vector of the first two principal components
% if zDif(2, 1) < zDif(1, 1)
%     wDif(:, 1) = -wDif(:, 1);
% end

% % Look at reward pts similarly
meanNeuralData = squeeze(mean(neuralData_byParameters_mean,[1]));
meanNeuralData = reshape(meanNeuralData, [nrewards*ndifficulties, nneurons]);
[wR,zR,eigVls_R] = pca(meanNeuralData-mean(meanNeuralData, 1),'numcomponents',nrewards*ndifficulties-1);
% check the score of high reward condition has the largest value on the first PC respectively
% if not, flip the eigen vector of the first two principal components
if zR(3, 1) < zR(1, 1)
    wR(:, 1) = -wR(:, 1);
end

if zR(4, 2) < zR(1, 2)
    wR(:, 2) = -wR(:, 2);
end


%% angle of wDif and wR
% dot_product = dot(wR(:,1), wDif(:,1));
% angle = acos(dot_product) * 180 / pi;
% disp(angle);

%% plot reward PC1 and difficulty PC1 by 2D scatter plot
neuralData_onDPC = (neuralData - mean(neuralData, 1)) * wR(:,1:2);
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
    h(i) = errorbar(X(i,:), Y(i,:), Yerr(i,:), Yerr(i,:), Xerr(i,:), Xerr(i,:), "Color", diffColors(i,:), 'LineWidth', 2); hold on;
    scatter(X(i,:), Y(i,:), 100, rewColors, "filled"); hold on;
end
set(gca, 'fontsize', 14, 'fontname', 'arial', 'tickdir', 'out');
xlabel("RDiff Axis: PC1 (Hz): "+round(eigVls_R(1)/sum(eigVls_R)*100, 2)+"%");
ylabel("RDiff Axis: PC2 (Hz): "+round(eigVls_R(2)/sum(eigVls_R)*100, 2)+"%");
legend(h, ["Tiny", "Huge"], Location="best")
saveas(gcf, outputFolder+"rewardAxis-2D-"+date+".jpg");

% %% plot reward PC1 and difficulty PC1
% for pc=1:2 % 1:2
%     Y = zeros(3,2); % reward x difficulty
%     Yerr = zeros(3,2);
%     if pc==1
%         neuralData_onDPC = (neuralData - mean(neuralData, 1)) * wR(:,1);
%     else
%         neuralData_onDPC = (neuralData - mean(neuralData, 1)) * wDif(:,1);
%     end

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
%     set(gca, 'fontsize', 14, 'fontname', 'arial', 'tickdir', 'out');
%     xticks([1 2 3])
%     xticklabels(["S", "M", "L"])
%     xlim([0.7 3.3])
%     legend(["Tiny", "Huge"], Location="best")
%     set(gcf,'position',[100,100,300,650]);
%     if pc == 1
%         ylabel("Reward Axis (Hz)");
%         saveas(gcf, outputFolder+"rewardAxis-"+date+".jpg");
%     else
%         ylabel("Difficulty Axis (Hz)");
%         saveas(gcf, outputFolder+"difficultyAxis-"+date+".jpg");
%     end
% end


% neuralData_onDPC = (neuralData - mean(neuralData, 1)) * wR(:,1);
% for i=1:nrewards
%     for j=1:ndifficulties
%         for k=1:ndirections
%             curInds = rewardLabels==rewards(i) & difficultyLabels==difficulties(j) & directionLabels==directions(k);
%             Yall{day, i, j, k} = neuralData_onDPC(curInds);
%         end
%     end
% end

%% time trajectory 
% figure;
% neuralData_onDPC = (neuralData - mean(neuralData, 1)) * wR(:,1);
% for r=1:nrewards
%     for diff=1:ndifficulties
%         curInds = rewardLabels == rewards(r) & difficultyLabels == difficulties(diff);
%         X = trialNums(curInds);
%         Y = movmean(neuralData_onDPC(curInds), 50);
%         p(diff) = plot(X, Y , Color=diffColors(diff,:), LineWidth=2);hold on;
%         scatter(X([1 end]), Y([1 end]), 100, rewColors(r, :), "filled"); hold on;
%     end
% end
% set(gca, 'fontsize', 14, 'fontname', 'arial', 'tickdir', 'out');
% legend(p, ["Tiny", "Huge"], Location="best")
% xlabel("Trials")
% ylabel("Reward Axis: PC1 (Hz)");
% saveas(gcf, outputFolder+"trajectory-"+date+".jpg");

close all;
end

%% plot Yall

% figure;
% for i=1:3
%     % scatter plot
%     % the marker is star
%     % the color is black
%     % the median is a red star
%     scatter(ones(size(Yall, 1), 1) * i, Yall(:, i), 30, 'k', "marker", "*"); hold on;
%     scatter(i, median(Yall(:, i)), 100, 'r', "marker", "*"); hold on;
% end
% % plot a black line on y = 0
% plot([0 4], [0 0], 'k', 'LineWidth', 2); hold on;
% xticks(1:3)
% xticklabels(["S", "M", "L"])
% xlim([0.7 3.3])
% set(gca, 'fontsize', 14, 'fontname', 'arial', 'tickdir', 'out');
% ylabel("Tiny - Huge on Reward Axis (Hz)");
% lim = round(max(Yall(:))) + 3;
% ylim([-1*lim lim])
% set(gcf,'position',[100,100,300,650]);
% saveas(gcf, outputFolder+"rewardAxis-GC-400-600.jpg");


%% plot Yall
% figure;
% for j=1:ndifficulties
%     Y = zeros(1,nrewards);
%     Yerr = zeros(1, nrewards);
%     for i=1:nrewards
%         Ys = cat(1, Yall{:, i, j});
%         Y(i) = mean(Ys);
%         Yerr(i) = std(Ys) / sqrt(length(Ys));
%     end
%     errorbar(1:nrewards, Y, Yerr, "Color", diffColors(j,:), 'LineWidth', 2); hold on;
% end
% set(gca, 'fontsize', 14, 'fontname', 'arial', 'tickdir', 'out');
% xticks([1 2 3])
% xticklabels(["S", "M", "L"])
% xlim([0.7 3.3])
% legend(["Tiny", "Huge"], Location="best")
% set(gcf,'position',[100,100,300,650]);
% ylabel("Reward Axis (Hz)");
% saveas(gcf, outputFolder+"rewardAxis-all.jpg");

rewardAxisRaw = Yall;
rewardAxisMean = zeros(ndate,3,2,8);
rewardAxisStd  = zeros(ndate,3,2,8);
for d=1:ndate
    for j=1:ndifficulties
        for i=1:nrewards
            for k=1:ndirections
                Ys = Yall{d,i,j,k};
                rewardAxisMean(d,i,j,k) = mean(Ys);
                rewardAxisStd(d,i,j,k)  = std(Ys) / sqrt(length(Ys));
            end
        end
    end
end
% save("rewardAxisAroundGCForAdam2.mat", "rewardAxisRaw", "rewardAxisMean", "rewardAxisStd");


%% plot
% X axis: the difference between reward = 3 and reward = 1 in difficulty = 0 condition
% Y axis: the difference between reward = 3 and reward = 1 in difficulty = 1 condition
% plot each direction and day separately (8*8 = 64 plots)
% also plot the line x = y
% then in the different plot, show the histogram of the difference between xaxis and yaxis

% figure;
% hold on;
% eachxy = zeros(ndate,8, 2);
% diffxy = zeros(ndate,8);
% for i = 1:ndate
%     for j = 1:ndirections
%         x = rewardAxisMean(i,3,1,j) - rewardAxisMean(i,1,1,j);
%         y = rewardAxisMean(i,3,2,j) - rewardAxisMean(i,1,2,j);
%         scatter(x, y, 100, 'k', "marker", "*"); hold on;
%         diffxy(i,j) = y - x;
%         eachxy(i,j,1) = x;
%         eachxy(i,j,2) = y;
%     end
% end
% plot([0 ylim+10], [0 ylim+10], 'k', 'LineWidth', 2); hold on;
% scatter(mean(eachxy(:,:,1), [1 2]), mean(eachxy(:,:,2), [1 2]), 250, 'r', 'marker', '*', 'LineWidth', 2);
% lx = xlabel("Reward axis range, Tiny targets (Spikes/s)");
% ly = ylabel("Reward axis range, Huge targets (Spikes/s)");
% lx.Color = diffColors(1,:);
% ly.Color = diffColors(2,:);
% set(gca, 'fontsize', 14, 'fontname', 'arial', 'tickdir', 'out');
% saveas(gcf, outputFolder+"rewardAxis-stretch.jpg");


% figure;
% hold on;
% Y = diffxy(:);
% % calculate p value by wilcoxon rank sum test
% p = signrank(Y);

% histogram(Y, 'FaceColor', 'none', 'EdgeColor', 'k', 'LineWidth', 2);
% % add median line with red color
% hold on;
% line([median(Y(:)) median(Y(:))], ylim, 'Color', 'r', 'LineWidth', 2);
% % show the median correlation coefficient
% text(median(Y(:)), 0.8*max(ylim), "Median: "+median(Y(:)), 'Color', 'r');
% % show the black dot line of correlation coefficient = 0
% line([1 1], ylim, 'Color', 'k', 'LineWidth', 2, 'LineStyle', '--');

% xlabel("Huge range - Tiny range");
% ylabel("Count");
% set(gca, 'fontsize', 14, 'fontname', 'arial', 'tickdir', 'out');
% title("p = "+p);
% saveas(gcf, outputFolder+"rewardAxis-stretch-hist.jpg");
% 
% 
% 
