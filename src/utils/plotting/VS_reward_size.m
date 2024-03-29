function VS_reward_size(data, rewardLabels, difficultyLabels, options)
    
    arguments
        data
        rewardLabels
        difficultyLabels
        options.Label
        options.OutputFolder
        options.Ylim = nan
        options.Size = nan
    end

    [rewardNames, rewardLegends, diffLegends, rewColors, diffColors, direColors, DiffStyle, DelayTimes, nDelayTimes] = getExperimentConstants();
    rewards = unique(rewardLabels); nrewards = length(rewards);
    difficulties = unique(difficultyLabels); ndifficulties = length(difficulties);

    figure; hold on;
    for j=1:ndifficulties
        Y = zeros(nrewards,1);
        Yerr = zeros(nrewards,1);
        for i=1:nrewards
            curInds = rewardLabels == rewards(i) & difficultyLabels == difficulties(j) & ~isnan(data);
            Y(i) = mean(data(curInds));
            Yerr(i) = std(data(curInds)) / sqrt(sum(curInds));
        end
        errorbar(1:nrewards, Y,Yerr, Color=diffColors(j, :), LineWidth=4);
    end
    set(gca, 'fontsize', 20, 'fontname', 'arial', 'tickdir', 'out', 'fontweight', 'bold');
    xlim([0.7, nrewards+0.3]); xticks(1:nrewards); xticklabels(rewardNames); ylabel(options.Label); legend(diffLegends, Location="best");
    if ~isnan(options.Ylim); ylim(options.Ylim); end;
    if isnan(options.Size)
        set(gcf,'position',[0,0,550,550]);
    else
        set(gcf,'position',[0,0,options.Size(1), options.Size(2)]);
    end
    saveas(gcf, options.OutputFolder+"-vs-Reward.jpg"); close all;
end