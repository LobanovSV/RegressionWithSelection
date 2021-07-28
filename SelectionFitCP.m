function [sigma,Rthr,Delta,p] = SelectionFitCP(R,Prec,vMin,vMax)
% Find optimal parameters sigma, Rthr and Delta, which minimise maximum 
% absolute difference D between the expected Fe(|R|) and observed Fo(|R|)
% cumulative probabilities. Here, sigma is the standard deviation of the
% expected probability density of the initial HD population (normal 
% distribution), Rthr is the selection threshold, and Delta is the
% selection width.
%
% Syntax:
%    [sigma,Rthr,Delta,p] = SelectionFitCP(R);
%    [sigma,Rthr,Delta,p] = SelectionFitCP(R,Prec);
%    [sigma,Rthr,Delta,p] = SelectionFitCP(R,Prec,vMin,vMax);
% Note: replace Prec, vMin, and/or vMax by [] to use default values.
%
% Input:
%    R    - residual ages at onset in the selected sample
%   Prec  - precision (accuracy) of numerical integration
%   vMin  - left boundary (minimum values) of the parameters space
%           v = [sigma; Rthr; Delta]
%   vMax  - right boundary (maximum values) of the parameters space
%           v = [sigma;Rthr; Delta]
% Output:
%  sigma  - standard deviation
%   Rthr  - selection threshold
%  Delta  - selection width
%    p    - one-sample Kolmogorov–Smirnov p-value
%
% Please cite the following paper if you use this code:
% Lobanov et al., bioRxiv 2021.07.16.452643. 
%
% Function is created by S. Lobanov in 2021.

%% Default parameters

if nargin < 2 || isempty(Prec)
    Prec = 1e-6;
end
if nargin < 3 || isempty(vMin)
    vMin = [6.7;17.3;3.2];
end
if nargin < 4 || isempty(vMax)
    vMax = [7.4;18.5;3.35];
end

%% Parameters

PrecNiter = 20;

%% Find cumulative probability of the residual ages at onset

R = abs(R);
r = reshape(sort(R),1,[]);
r = [r;r];
r = r(:);
Fo = linspace(0,1,length(R)+1);
Fo = [Fo(1:end-1);Fo(2:end)];
Fo = Fo(:);

%% Main loop

Dv = vMax - vMin;
D = inf;
F = figure;
repl = 0;
Legend = cell(1e3,3);
iL = 0;
Quantity = {'\sigma','R_{thr}','\Delta'};
while ishandle(F)
    repl = repl + 1;
    vi = vMin + rand(3,1) .* Dv;
    Di = FindErr(r,Fo,Prec,vi);
    for i = 4:3+PrecNiter
        dv = 2^(-i) * Dv .* [-1,1];
        Modified = true;
        while Modified
            Modified = false;
            for j = 1:3
                if dv(j,1) == 0
                    continue
                end
                vj = vi;
                for k = 1:2
                    vj(j) = max(vMin(j),min(vMax(j),vi(j) + dv(j,k)));
                    
                    Dj = FindErr(r,Fo,Prec,vj);
                    
                    if Dj < Di
                        Di = Dj;
                        vi = vj;
                        if k == 2
                            dv(j,:) = dv(j,[2,1]);
                        end
                        Modified = true;
%                         FindErr(r,Fo,Prec,vi,R,true);
%                         drawnow
                        break
                    end
                end
            end
        end
    end
    if Di < D
        D = Di;
        v = vi;
        iL = iL + 1;
        for j = 1:3
            subplot(1,3,j)
            plot(v(j),D,'.')
            hold on
            xlabel(Quantity{j})
            ylabel('D_{min}')
            Legend{iL,j} = [int2str(repl) ': ' Quantity{j} ' = ' ...
                num2str(v(j),3)];
            legend(Legend(1:iL,j))
        end
    end
    title(['Replication ' int2str(repl) ...
        '. Close figure to stop computation'])
    drawnow
end

%% Find p-value

[~,p] = FindErr(r,Fo,Prec,v,R,true);

%% Extract parameters

sigma = v(1);
Rthr = v(2);
Delta = v(3);

%% Check that parameters are not on the boundary

if any((v == vMin | v == vMax) & vMin ~= vMax)
    warning('Some estimates are on the boundary. Move this boundary.')
end

end

function [D,p] = FindErr(r,Fo,Prec,v,R,Draw)
if nargin < 6 || isempty(Draw)
    Draw = false;
end

sigma = v(1);
Rthr = v(2) / sigma;
Delta = v(3) / sigma;
r = r / sigma;

Zmax = 2 * sqrt(-log(Prec));
Z = (0:Prec:Zmax).';
ya = 2 / sqrt(2*pi) * exp(-1/2 * Z.^2);
S = 1 ./ (1 + exp((Rthr-Z)/Delta));
pa = ya.*S;
pa = pa / trapz(pa);
Fe = cumsum(pa) - pa(1) / 2;
if any(r > Zmax)
    Z = [Z;unique(r>Zmax)];
    Zmax = max(Z);
    Fe = [Fe;ones(length(Z)-length(Fe),1)];
end
D = max(abs(interp1(Z,Fe,r,'linear')-Fo));

if Draw
    Z = Z * sigma;
    r = r * sigma;
    Zmax = Zmax * sigma;
    [~,p] = kstest(R,'CDF',[Z,Fe]);
    if Draw
        plot(r,Fo,'b.',Z,Fe,'r')
        grid on
        ylim([0,1])
        xlim([0,Zmax])
        xlabel('Absolute value of the residual age at onset (years)')
        ylabel('Cumulative probability')
        legend('Observed','expected')
    end
end
end