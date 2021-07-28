function [sigma,b0,b,m2l] = SelectionFindML(R,RL,IsLR,Prec,v,vMin,vMax, ...
    Rthr,Delta)
% Find maximum likelihood L and estimates of the unknown standard
% deviation sigma, intercept b0 and effect size b. Note, m2l = -2*log(L).
%
% Syntax:
%    [sigma,b0,b,m2l] = SelectionFindML(R,RL,IsLR);
%    [sigma,b0,b,m2l] = SelectionFindML(R,RL,IsLR,Prec);
%    [sigma,b0,b,m2l] = SelectionFindML(R,RL,IsLR,Prec,v);
%    [sigma,b0,b,m2l] = SelectionFindML(R,RL,IsLR,Prec,v,vMin,vMax);
% Note: replace Prec, v, vMin, and/or vMax by [] to use default values.
% Input:
%    R    - residual ages at onset
%    RL   - repeat lengths
%   IsLR  - logical value indicating whether linear regression must be
%           applied to corresponding residual age at onset
%   Prec  - precision (accuracy) of numerical integration
%    v    - initial estimates: v = [sigma;b0;b]
%   vMin  - left boundary (minimum values) of the parameters space
%   vMax  - right boundary (maximum values) of the parameters space
%   Rthr  - selection threshold
%  Delta  - selection width
% Output:
%  sigma  - standard deviation
%    b0   - intercept
%    b    - effect size
%   m2l   - -2*log(L), where L is likelihood
%
% Please cite the following paper if you use this code:
% Lobanov et al., bioRxiv 2021.07.16.452643.
%
% Function is created by S. Lobanov in 2021.

%% Default parameters

if nargin < 4 || isempty(Prec)
    Prec = 1e-6;
end
if nargin < 6 || isempty(vMin)
    vMin = [5;-5;-5];
end
if nargin < 7 || isempty(vMax)
    vMax = [10;5;5];
end
if nargin < 5 || isempty(v)
    v = (vMin + vMax) / 2;
end
if nargin < 8 || isempty(Rthr)
    Rthr = 17.832094282432730;
end
if nargin < 9 || isempty(Delta)
    Delta = 3.285861129837341;
end

%% Parameters

PrecNiter = 20;

%% Find unique RL

% Remove NaNs
ii = isnan(R) | isnan(RL);
if any(ii)
    RL(ii) = [];
    R(ii) = [];
    IsLR(ii) = [];
end

[RL,~,ii] = unique(RL);
n = zeros(size(RL));
for i = 1:length(RL)
    n(i) = sum(ii == i);
end
[~,ii] = sort(ii);
R = R(ii);
IsLR = IsLR(ii);
R = mat2cell(R,n);
IsLR = mat2cell(IsLR,n);

clear n ii i

%% Main loop

m2l = SelectedRAAO_m2l(R,RL,IsLR,v,Prec,Rthr,Delta);
for i = 4:3+PrecNiter
    dv = 2^(-i) * (vMax-vMin) .* [-1,1];
    Modified = true;
    while Modified
        Modified = false;
        for j = 1:3
            if dv(j,1) == 0
                continue
            end
            vj = v;
            for k = 1:2
                vj(j) = max(vMin(j),min(vMax(j),v(j) + dv(j,k)));
                
                m2lj = SelectedRAAO_m2l(R,RL,IsLR,vj,Prec,Rthr,Delta);
                
                if m2lj < m2l
                    m2l = m2lj;
                    v = vj;
                    if k == 2
                        dv(j,:) = dv(j,[2,1]);
                    end
                    Modified = true;
                    break
                end
            end
        end
    end
end

%% Extract parameters

sigma = v(1);
b0 = v(2);
b = v(3);

%% Check that parameters are not on the boundary

if any((v == vMin | v == vMax) & vMin ~= vMax)
    warning('Some estimates are on the boundary. Move this boundary.')
end

end

function m2l = SelectedRAAO_m2l(R,RL,IsLR,v,Prec,Rthr,Delta)
% Find m2l = -2 * log(L) for the HD cohort with selection.
% RL is vector, RAAO is cell of the same size as RL.

%% Parameters of the selection function

% Rthr = 17.793524124464870;
% Delta = 3.288084607124466;

%% Calculate

sigma = v(1);
b0 = v(2);
b = v(3);

Zmax = 2 * sqrt(-log(Prec));
NZ = ceil(Zmax/Prec);
Z = linspace(-Zmax,Zmax,2*NZ+1);
N = Zmax / NZ / sqrt(2*pi) * exp(-0.5 * Z.^2);

m2l = 0;
for i = 1:length(RL)
    delta = b0 + b * RL(i);
    if i== 1 || b~=0
        Norm = sum(N ./ (1 + exp((Rthr-abs(sigma * Z + delta))/Delta)));
    end
    
    f = 1 ./ (1 + exp((Rthr-abs(R{i}))/Delta));
    f(IsLR{i}) = Norm;
    
    L = 1/sqrt(2*pi)/sigma/Norm * exp(-0.5/sigma^2 * (R{i}-delta).^2) .* f;
    
    m2l = m2l - 2 * sum(log(L));
end
end