function [W, W0, r, best_lambda, FitInfo] = findVigorAxisCV(neuralData, vigorData, kfolds, outputFileName)

    % Lasso regression with cross-validation
    lambdas = logspace(-4, 1, 240);
    [B, FitInfo] = lasso(neuralData, vigorData, 'CV', kfolds, 'Lambda', lambdas, 'alpha', 0.0001);

    % Best lambda
    % best_lambda = FitInfo.LambdaMinMSE;
    best_lambda = 0.5;
    [~, index] = min(abs(lambdas - 0.01));

    % Intercept and coefficients corresponding to best lambda
    
    % index = FitInfo.IndexMinMSE;
    W0 = FitInfo.Intercept(index);
    W = B(:,index);

    % Residuals
    r = vigorData - (neuralData * W + repmat(W0, size(vigorData, 1), 1));

    % calculate R square
    R2 = 1 - sum(r.^2) / sum((vigorData - mean(vigorData)).^2);

    % Plot lambda vs MSE
    loglog(FitInfo.Lambda, FitInfo.MSE, 'o-');
    hold on;
    % loglog(best_lambda, FitInfo.MSE(index), 'ro');
    xlabel('Lambda');
    ylabel('MSE');
    title('Cross-validation result');
    % legend('MSE','Lambda with minimum MSE','Location','best');
    hold off;
    saveas(gcf, outputFileName + "_lambda_vs_MSE.png");

    % Plot Vigor Axis Projection VS residuals
    figure;
    plot(neuralData * W + repmat(W0, size(vigorData, 1), 1), r, 'o');
    xlabel('Vigor Axis Projection');
    ylabel('Residuals');
    title('Vigor Axis Projection VS residuals: R^2 = ' + string(R2));
    saveas(gcf, outputFileName + "_VigorAxisProjection_vs_residuals.png");
    close all;
end

% function [W, W0, r, best_lambda, FitInfo] = findVigorAxisCV(neuralData, vigorData, kfolds)
% this function is used to check the overfitting of the model
% do cross-validation and see the Square error of the model
% in the training set and the test set
% Then find the best regularization parameter based on the test error
% across kfolds
% Input:
% - neuralData: (M by N) matrix, M is the number of trials,
% N is the number of neurons
% neuralData is already averaged across trials
% - vigorData: (M by 1) vector, M is the number of trials
% - kfolds: (1 by 1) scalar, the number of folds
% Output: 
% - W: (N by 1) vector, the axis which is correlated to vigor
% - W0: (1 by 1) scalar, the intercept of the linear model
% - r: (M by 1) vector, the residual of the linear model
% - best_lambda: (1 by 1) scalar, the best regularization parameter
% - FitInfo: (1 by 1) struct, the result of cross-validation
% Algorithm:
% - split the data into kfolds (randomly)
% - for each fold, fit a linear model on the training set
% - loss function: sum of square error + regularization term
% - regularization term: Lasso (L1) regularization
% - calculate the square error of the model on the training set
% - calculate the square error of the model on the test set
% - average the square error across kfolds
% - find the best regularization parameter based on the test error