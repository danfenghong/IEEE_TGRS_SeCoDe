% Demo_Samson
close all;
clc, clear;

addpath(genpath(pwd));

load samson_1.mat;
Y = V;

load end3.mat;
h = nRow; w = nCol;
l = nBand;
clear nRow nCol nBand SlectBands maxValue

k = 3;
S = A; %abundance
A = M; %mixing matrix
clear M

%% endmember initialization
load end3_extraction.mat;
A_0 = EM;
clear EM
s=1-pdist2(A_0',A','cosine');
[~,in]=sort(s,'descend');
index=in(1,:);
A_0 = A_0(:,index);

%% abundance initialization
tic
[~, X_0] = SPCLSU_ADMM(Y,A_0,1000);
X_0=X_0./repmat(sum(X_0),size(X_0,1),1);
toc

%% SeCoDE
% Set up cbpdndl parameters
lambda = 1e-2;
opt = [];
opt.Verbose = 1;
opt.MaxMainIter = 20;
opt.rho = 100*lambda + 0.5;
opt.sigma = 5;
opt.AutoRho = 1;
opt.AutoRhoPeriod = 10;
opt.AutoSigma = 1;
opt.AutoSigmaPeriod = 10;
opt.XRelaxParam = 1.8;
opt.DRelaxParam = 1.8;
opt.NonNegCoef = 1;

% initialize F
% with pretrained D
load('ConvDict.mat');
dmap = containers.Map(ConvDict.Label, ConvDict.Dict);
D = dmap('12x12x36');
p = 12;

mu = 1e1;
F = D;

delta = 1e-1;
Y_a = [Y;delta*ones(1,h*w)];

for iter = 1:6

    % MC-CSC
    X_3d = reshape(X_0', [h,w,k]);
    npd = p/2;
    fltlmbd = 10;
    [Xl, Xh] = lowpass(X_3d, fltlmbd, npd);

    [F, M, FMh, FM, optinf] = cbpdndl_unmixing(F, Xh, Xl, S, A_0, Y, lambda, opt);
    Z = reshape(FM, size(X_0'))'; % k-N
    Z = double(Z);
    
    % CSU
    [A_0, X_0] = NMF_ADMM(Y, Z, 1e-1, 0.5, 1000);
    
end

Z1 = Z./repmat(sum(Z,1),[k,1]); 
[ARMSE(Z1, S), SadEval(A_0, A, k), OA_CM(Z1, S)]
