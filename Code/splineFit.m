function [mse_x, mse_y, rmse, coeffs_x, coeffs_y] = splineFit(numPieces, tVec, x, y, plotFlagSpline)

t0=tVec(1);
tf=tVec(end);
% tVec;

breaks = linspace(t0, tf, numPieces + 1);
knotVector = augknt(breaks, 4);

% Fit spline to x(t)
sp_x = spap2(knotVector, 4, tVec, x);

% Fit spline to y(t)
sp_y = spap2(knotVector, 4, tVec, y);

% Coefficients for x(t)
coeffs_x = sp_x.coefs.';

% Coefficients for y(t)
coeffs_y = sp_y.coefs.';

numCoeffs = length(coeffs_x);  % Should be the same for x and y
% fprintf('Number of coefficients: %d\n', numCoeffs);

sp_x_pred = spmak(knotVector, coeffs_x.');
sp_y_pred = spmak(knotVector, coeffs_y.');

% Evaluate the spline at tVec
% x_fit = fnval(sp_x, tVec).';
% y_fit = fnval(sp_y, tVec).';
x_fit = fnval(sp_x_pred, tVec).';
y_fit = fnval(sp_y_pred, tVec).';
% Calculate MSE
mse_x = mean((x - x_fit).^2);
mse_y = mean((y - y_fit).^2);

% Calculate RMSE
rmse = sqrt(mse_x + mse_y);
% fprintf('RMSE: %f\n', rmse);

if plotFlagSpline == 1
    if length(tVec) > 99
        figure(3);
    else
        figure(10);
    end
    hold on
    % plot(x, y, 'o', 'DisplayName', 'Original Data');
    hold on;
    plot(x_fit, y_fit, 'r-', 'LineWidth', 3, 'DisplayName', 'Spline Fit');
    legend('show');
    xlabel('x');
    ylabel('y');
    title('Spline Fit vs. Original Data');
    grid on;
    axis equal
end