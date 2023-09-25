function VS_focus_difficulty(data, rewardLabels, difficultyLabels, focusBlockLabels, focusDifficultyLabels, options)

    arguments
        data (:,1) double
        rewardLabels (:,1) double
        difficultyLabels (:,1) double
        focusBlockLabels (:,1) double
        focusDifficultyLabels (:,1) double
        options.OutputFolder (1,:) char
        options.Label (1,:) char
        options.Legend = true
        options.DiffLegendPos = "southwest"
        options.RewLegendPos = "south"
    end

    [rewardNames, rewardLegends, diffLegends, rewColors, diffColors, direColors, DiffStyle, DelayTimes, nDelayTimes] = getExperimentConstants();
    rewards = unique(rewardLabels); nrewards = length(rewards);
    difficulties = unique(difficultyLabels); ndifficulties = length(difficulties);
    focusBlocks = unique(focusBlockLabels); nfocusBlocks = length(focusBlocks);
    focusDifficulties = unique(focusDifficultyLabels); nfocusDifficulties = length(focusDifficulties);

    figure; hold on;

    for h=1:nfocusBlocks
        for j=1:ndifficulties
            for i=1:nrewards
                curInds = focusBlockLabels == focusBlocks(h) & rewardLabels == rewards(i) & difficultyLabels == difficulties(j);
                tmp = focusDifficultyLabels(curInds);
                tmp_focusDifficulties = unique(tmp);
                X = find(ismember(focusDifficulties, unique(tmp)));
                Y = zeros(length(tmp_focusDifficulties), 1);
                for k=1:length(tmp_focusDifficulties)
                    curInds = focusBlockLabels == focusBlocks(h) & rewardLabels == rewards(i) & difficultyLabels == difficulties(j) ...
                             & focusDifficultyLabels == tmp_focusDifficulties(k) & ~isnan(data);
                    Y(k) = mean(data(curInds));
                    Yerr(k) = std(data(curInds)) / sqrt(sum(curInds));
                end
                errorbar(X, Y,Yerr, 'Color', rewColors(rewards(i), :), 'LineWidth', 2.5,  'LineStyle', DiffStyle(j));
                lh(j) = plot(0, 0, 'Color', 'k', 'LineWidth', 2.5,  'LineStyle', DiffStyle(j));
                lr(i) = plot(0, 0, 'Color', rewColors(rewards(i), :), 'LineWidth', 2.5);
            end
        end
    end
    set(gca, 'fontsize', 20, 'fontname', 'arial', 'tickdir', 'out', 'fontweight', 'bold');
    xlim([0.7 nfocusDifficulties+.3]); xticks(1:nfocusDifficulties);
    xticklabels(focusDifficulties);
    xlabel("Focus Difficulty"); ylabel(options.Label); set(gcf,'position',[0,0,550,550]);
    if options.Legend
        ah1=axes('position',get(gca,'position'),'visible','off');
        legend(ah1,lr,rewardLegends, Location=options.RewLegendPos);
        legend(lh, ["Tiny", "Huge"], Location=options.DiffLegendPos);
    end
    set(gca, 'fontsize', 20, 'fontname', 'arial', 'tickdir', 'out', 'fontweight', 'bold');
    saveas(gcf, options.OutputFolder+"-vs-focus-difficulty.jpg");
    close all;
end