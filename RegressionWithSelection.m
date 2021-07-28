% This script performs regression with selection as described in Section S2
% of Lobanov et al., bioRxiv 2021.07.16.452643.

clc
clear
close all

%% Parameters

% Precision (accuracy) of numerical integration
Prec = 1e-5;

% Left boundary (minimum values) of the parameters space v = [sigma; b0; b]
vMin = [5; -5; -5];

% Right boundary (maximum values) of the parameters space v = [sigma;b0;b]
vMax = [15; 5; 5];

% Selection threshold
Rthr = 17.832094282432730;

% Selection width
Delta = 3.285861129837341;

% Full path to the 'Additional file 6.xlsx'. This file will be available
% once publication is accepted and deposited on the journal website.
% Note, bioRxiv version does not contain this file.
FileLoad = 'Additional file 6.xlsx';

%% Read data

IHD = readtable(FileLoad);
IHD.RL = IHD.dNmax_QTR + IHD.dNmin_QTR;
IHD.IsLR = strcmp(IHD.Study,'MGH');
IHD(isnan(IHD.RL),:) = [];

%% Perform regression with selection

% Linear regression (LR)
mdl = fitlm(IHD,'RAAO ~ RL');
sigmaLR = mdl.RMSE;
b0LR = mdl.Coefficients.Estimate(1);
bLR = mdl.Coefficients.Estimate(2);
pLR = mdl.Coefficients.pValue(2);

% H1
[sigma,b0,b,m2l] = SelectionFindML(IHD.RAAO,IHD.RL,IHD.IsLR,Prec, ...
    [],vMin,vMax,Rthr,Delta);

% H0
vMin(3) = 0;
vMax(3) = 0;
[sigmaH0,b0H0,~,m2lH0] = SelectionFindML(IHD.RAAO,IHD.RL,IHD.IsLR, ...
    Prec,[],vMin,vMax,Rthr,Delta);

% Find p-value using likelihood ration test
p = chi2cdf(m2lH0-m2l,1,'upper');

%% Create summary table T

T = table({'sigma';'b0';'b';'p'},[sigmaLR;b0LR;bLR;pLR], ...
    [sigmaH0;b0H0;NaN;NaN],[sigma;b0;b;p],'VariableNames', ...
    {'Variable';'LR';'H0';'H1'});