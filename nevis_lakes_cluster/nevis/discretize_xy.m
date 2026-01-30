function [xd, yd, xs, ys] = discretize_xy(X, Y, E)
% Calculate binned Y values based on a discretized x variable
% Laura A. Stevens 05 Apr 2017
% INPUTS
% X             X variable, 2-D matrix
% Y             Y variable, 2-D matrix
% E             Number of edges for your bins (=number of bins+1)
% OUTPUTS
% xd            mean of binned x-values
% xs            standard deviation of binned x-values
% yd            mean of y-values that fall within the x bins
% ys            standard deviation of y-values that fall within x bins

N = linspace(min(min(X)), max(max(X)), E);          % assign bin edges
I = discretize2(X, N);                               % discretize x-axis

for i=1:(length(N)-1);
    I2 = I; I2 (I2 ~= i) = NaN; I3 = I2./I2;
    interestX = (I3).*X;                            % select X of interest
    xd(i,:) = nanmean(nanmean(interestX));
    xs(i,:) = nanstd(nanstd(interestX));
    interestY = (I3).*Y;                            % select Y of interest
    yd(i,:) = nanmean(nanmean(interestY));
    ys(i,:) = nanstd(nanstd(interestY));
end
end
