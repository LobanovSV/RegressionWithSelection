% This script estimates parameters of the selection function as described
% in Section S2 of Lobanov et al., bioRxiv 2021.07.16.452643.
% Note, the function SelectionFitCP performs replications until you close
% the figure. The longer you wait, the more chance that you find global
% minimum.

clc
clear
close all

%% Parameters

% Precision (accuracy) of numerical integration
Prec = 1e-5;

% Left boundary (minimum values) of the parameters space v = [sigma; Rthr;
% Delta]
vMin = [6.7; 17.3; 3.2];

% Right boundary (maximum values) of the parameters space v = [sigma; Rthr;
% Delta]
vMax = [7.4; 18.5; 3.35];

% Full path to the 'Additional file 6.xlsx'. This file will be available
% once publication is accepted and deposited on the journal website.
% Note, bioRxiv version does not contain this file.
FileLoad = 'Additional file 6.xlsx';

%% Read data

IHD = readtable(FileLoad);
IHD = IHD(strcmp(IHD.Study,'Cardiff'),:);

%% Estimate parameters of the selection function

[sigma,Rthr,Delta,p] = SelectionFitCP(IHD.RAAO,Prec,vMin,vMax);
% sigma - standard deviation
% Rthr - selection threshold
% Delta - selection width
% p - one-sample Kolmogorov–Smirnov p-value