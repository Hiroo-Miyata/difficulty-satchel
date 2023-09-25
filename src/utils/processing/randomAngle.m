function [angles]  = randomAngle(ndim, nrep)
    % Generate a random n-dimensional vector with normal distribution
    eigenVec = randn(ndim, 1);

    % Normalize it to have length 1 (making it a vector on the n-sphere)
    eigenVec = eigenVec / norm(eigenVec);

    angles = zeros(nrep, 1);
    randVec = randn(ndim, nrep);
    for i = 1:nrep
        randVec(:, i) = randVec(:, i) / norm(randVec(:, i));
        angles(i) = acos(randVec(:, i)' * eigenVec);
        angles(i) = angles(i) * 180 / pi;
        if angles(i) > 90
            angles(i) = 180 - angles(i);
        end
    end

    % Plot the histogram
    % with 95% confidence interval
    figure;
    histogram(angles, 20);
    title('Histogram of angles between random vectors and eigenvector');
    xlabel('Angle (degrees)');
    ylabel('Frequency');
    hold on;
    xline(mean(angles), 'r', 'LineWidth', 2);
    xline(mean(angles) + 2 * std(angles), 'r--', 'LineWidth', 2);
    xline(mean(angles) - 2 * std(angles), 'r--', 'LineWidth', 2);
    hold off;
end