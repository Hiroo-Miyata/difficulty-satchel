function VS_reward_size(data, rewardLabels, difficultyLabels, label, outputFolder)
    nrewards = 2; rewards = [1 3]; rewardNames = ["S", "L"];
    ndifficulties = 2; difficulties = [0 1];
    rewColors = [1 0 0; 1 0.6470 0; 0 0 1]; diffColors = [0 0.447 0.741; 0.466 0.674 0.188];
    direColors = {[1 .5 .5],[.75 .75 .5],[.5 1 .5],[.25 .75 .5],[0 .5 .5],[0.25 0.25 .5],[0.5 0 .5],[0.75 0.25 .5]};
    DiffStyle = ["-", ":"];

    figure; hold on;
    for j=1:ndifficulties
        Y = zeros(nrewards,1);
        Yerr = zeros(nrewards,1);
        for i=1:nrewards
            curInds = rewardLabels == rewards(i) & difficultyLabels == difficulties(j) & ~isnan(data);
            Y(i) = mean(data(curInds));
            Yerr(i) = std(data(curInds)) / sqrt(sum(curInds));
        end
        errorbar(1:nrewards, Y,Yerr, Color=diffColors(j, :), LineWidth=2);
    end
    set(gca, 'fontsize', 20, 'fontname', 'arial', 'tickdir', 'out', 'fontweight', 'bold');
    xlim([0.7, nrewards+0.3]); xticks(1:nrewards); xticklabels(rewardNames); ylabel(label); legend(["Tiny", "Huge"], Location="best");
    set(gcf,'position',[0,0,550,550]);
    saveas(gcf, outputFolder+"-vs-Reward.jpg"); close all;
end