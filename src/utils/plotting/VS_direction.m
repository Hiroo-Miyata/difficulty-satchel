function VS_direction(data, rewardLabels, difficultyLabels, directionLabels, options)
    
    arguments
        data
        rewardLabels
        difficultyLabels
        directionLabels
        options.Label
        options.OutputFolder
        options.DiffLegendPos = "southwest"
        options.RewLegendPos = "south"
    end

    [rewardNames, rewardLegends, diffLegends, rewColors, diffColors, direColors, DiffStyle, DelayTimes, nDelayTimes] = getExperimentConstants();
    rewards = unique(rewardLabels); nrewards = length(rewards);
    difficulties = unique(difficultyLabels); ndifficulties = length(difficulties);
    directions = unique(directionLabels); ndirections = length(directions);

    figure; hold on;
    for j=1:ndifficulties
        for i=1:nrewards
            Y = zeros(ndirections, 1); Yerr = zeros(ndirections, 1);
            for k=1:ndirections
                curInds = rewardLabels == rewards(i) & difficultyLabels == difficulties(j) & directionLabels == directions(k) & ~isnan(data);
                Y(k) = mean(data(curInds));
                Yerr(k) = std(data(curInds)) / sqrt(sum(curInds));
            end
            errorbar(1:ndirections, Y,Yerr, 'Color', rewColors(rewards(i), :), 'LineWidth', 2.5,  'LineStyle', DiffStyle(j));
            lh(j) = plot(0, 0, 'Color', 'k', 'LineWidth', 2.5,  'LineStyle', DiffStyle(j));
            lr(i) = plot(0, 0, 'Color', rewColors(rewards(i), :), 'LineWidth', 2.5);
        end
    end
    set(gca, 'fontsize', 20, 'fontname', 'arial', 'tickdir', 'out', 'fontweight', 'bold');
    xlim([0.7 ndirections+.3]); xticks(1:ndirections); xticklabels([0 45 90 135 180 225 270 315]);
    xlabel("Direction"); ylabel(options.Label); set(gcf,'position',[0,0,550,550]);
    ah1=axes('position',get(gca,'position'),'visible','off');
    legend(ah1,lr,rewardLegends, Location=options.RewLegendPos);
    legend(lh, diffLegends, Location=options.DiffLegendPos); 
    set(gca, 'fontsize', 20, 'fontname', 'arial', 'tickdir', 'out', 'fontweight', 'bold');
    saveas(gcf, options.OutputFolder+"-vs-Direction.jpg");
    close all;
end