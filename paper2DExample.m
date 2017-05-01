%example with rho = 0.792888603593367, 5 modes, n = 4
%clear all;
%close all;
%clc;
addpath(genpath('YALMIP'));
addpath(genpath('JSR'));
addpath(genpath('sdpt3'));

%options for the solver
ops = sdpsettings('solver','sdpt3');

% dimension of the state space
n = 2;

%number of modes
m = 4; 

A{1}=[0.6498 0.0576;
    0.0576  0.6313];

A{2}=[0.1795 0.1092;
    0.1092 0.2093];

A{3}=  [0.6427 0.0016;
    0.0016 0.6335];

A{4}=[0.1062 0.0729;
    0.0729 0.7647];
%Campi's certainty
beta = 0.92;
N = 1000;

[lowerBound, upperBound] = computeRhoSequentialSampling(A, beta, N);