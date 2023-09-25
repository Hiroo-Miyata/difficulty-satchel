function [neuralActivity, nneurons, ntimebins, ntrials] = kernelSmooth(neuralActivityRaw)

    arguments
        neuralActivityRaw (:, :, :) double
    end
    [nneurons, ntimebins, ntrials] = size(neuralActivityRaw);

    neuralActivity = zeros(nneurons, ntimebins, ntrials);
    % apply guassian kernel
    for i = 1:nneurons
        for j = 1:ntrials
            sigma = 25;
            kernel = normpdf(-3*sigma:3*sigma,0,sigma);
            neuralActivity(i, :, j) = conv(squeeze(neuralActivityRaw(i, :, j)), kernel, "same");
        end
    end
end