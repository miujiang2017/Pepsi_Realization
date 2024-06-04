function [fitresult, gof] = createFit_stdR_res(a, g)
%CREATEFIT(A,G)
%  Create a fit.
%
%  Data for 'untitled fit 1' fit:
%      X Input : a
%      Y Output: g
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  另请参阅 FIT, CFIT, SFIT.

%  由 MATLAB 于 23-Apr-2023 20:01:59 自动生成


%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( a, g );

% Set up fittype and options.
ft = fittype( 'poly6' );

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft );

% Plot fit with data.
% figure( 'Name', 'untitled fit 1' );
% h = plot( fitresult, xData, yData );
% legend( h, 'g vs. a', 'untitled fit 1', 'Location', 'NorthEast', 'Interpreter', 'none' );
% % Label axes
% xlabel( 'a', 'Interpreter', 'none' );
% ylabel( 'g', 'Interpreter', 'none' );
% grid on


