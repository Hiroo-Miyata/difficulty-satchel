function VS_time(data, rewardLabels, difficultyLabels, options)
    
    arguments
        data (:,:) double
        rewardLabels
        difficultyLabels
        options.Label
        options.OutputFolder
        options.Xlim (1,2) double = nan
        options.Ylim = nan
        options.Timeperiod = nan
    end

    [rewardNames, rewardLegends, rewColors, diffColors, direColors, DiffStyle, DelayTimes, nDelayTimes] = getExperimentConstants();
    rewards = unique(rewardLabels); nrewards = length(rewards);
    difficulties = unique(difficultyLabels); ndifficulties = length(difficulties);

    [timeBin, ntrials] = size(data);
    if isnan(options.Xlim)
        options.Xlim = [0, timeBin];
        Xticks = [0, timeBin];
    else
        if options.Xlim(1) < 0
            Xticks = []
        else
            Xticks = [0, timeBin(end)];
        end
    end


    figure; hold on;
    for j=1:ndifficulties
        for i=1:nrewards
            Y = zeros(1, timeBin); Yerr = zeros(1, timeBin);
            curInds = rewardLabels == rewards(i) & difficultyLabels == difficulties(j) & ~isnan(data);
            Y = mean(data(:, curInds), 2);
            % calculate 95%CI
            Yerr = 1.96*std(data(:, curInds), 0, 2)/sqrt(sum(curInds));

            % plot
            plot(timeBin, Y, 'Color', rewColors(rewards(i), :), 'LineWidth', 2.5,  'LineStyle', DiffStyle(j));
            fill([1:timeBin, timeBin:-1:1], [Y+Yerr, Y(end:-1:1)-Yerr(end:-1:1)], rewColors(rewards(i), :), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        end
    end
    set(gca, 'fontsize', 20, 'fontname', 'arial', 'tickdir', 'out', 'fontweight', 'bold');
    xlim(options.Xlim); xlabel("Time (ms)"); 
    ylabel(options.Label); legend(["Tiny", "Huge"], Location="best");
    if ~isnan(options.Ylim); ylim(options.Ylim); end;
    set(gcf,'position',[0,0,550,550]);
    saveas(gcf, options.OutputFolder+"-vs-Reward.jpg"); close all;
end