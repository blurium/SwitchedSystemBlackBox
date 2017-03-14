%example with rho = 0.958096329922833
%5 modes, n = 4
clear all;
close all;
clc;
addpath(genpath('YALMIP'));
addpath(genpath('JSR'));
addpath(genpath('sdpt3'));

% dimension of the state space
n = 4;

%number of modes
m = 5; 

A{1}=[-0.1197   -0.4672    0.2844   -0.0927;
   -0.4672    0.2110    0.1060   -0.1289;
    0.2844    0.1060    0.4887   -0.2726;
   -0.0927   -0.1289   -0.2726   -0.4656];

A{2}=[0.5977   -0.1345   -0.5015    0.2847;
   -0.3523    0.2387   -0.2387   -0.5257;
    0.1443   -0.6623    0.1582   -0.2652;
   -0.4535    0.0277   -0.4704    0.4822];

A{3}=[0.1465    0.4447   -0.2367    0.0089;
   -0.1938    0.4535    0.5118   -0.2639;
   -0.4632   -0.1213   -0.0623   -0.3205;
    0.0413    0.3958   -0.1186    0.2844];

A{4}=[-0.7974    0.1273   -0.0250   -0.1665;
   -0.0708   -0.7925    0.1024    0.1787;
    0.1820   -0.1764   -0.7208   -0.0708;
    0.0799    0.0098   -0.2412   -0.6642];

A{5}=[0.2753   -0.1374    0.2820   -0.4957;
    0.2781    0.3644    0.4157    0.2455;
   -0.2612   -0.4200    0.4171    0.1872;
   -0.4456    0.3394    0.1653   -0.1120];

%bounds to generate the sampled points x_i
uBound = 4;
lBound = -4;

%options for the solver
ops = sdpsettings('solver','sdpt3');

%number of points in the sampling
for N=20:20:1000


%generate the N points
for i=1:N
    X{i} = lBound + (uBound-lBound)*rand(1,n);
    X{i}=X{i}';
end

%apply randomly one mode to every point, store the result in array Y
for i=1:N
    k=unidrnd(m);
    Y{i}=A{k}*X{i};
end

%computation of lower bound by bisection, we start with lambda = 1, and we stop after 50 iterations
lambdaNext = 1;
lambdaU = 1;
lambdaL = 0;
epsilon = 0.005;

P_var = sdpvar(n,n); %optimization variable
Objective = 1;
counter = 0;
lambdaUfound = false;

while(counter<50 && lambdaUfound==false)

    Constraints = [];
    lambda = lambdaNext;
    Constraints = Constraints + (P_var > 0);
    for i=1:N
    Constraints = Constraints + (Y{i}'*P_var*Y{i} <= lambda^2*X{i}'*P_var*X{i});
    end
    sol = optimize(Constraints,Objective,ops);
    if sol.problem == 0
       lambdaU = lambda;
       lambdaL = lambda/2;
       lambdaUfound=true;
    else
        lambdaNext=2*lambda;
    end
    counter=counter+1;
    display('try again!');
end

if lambdaUfound==false
    display('lower bound too high')
end

    feasible = true;
lambdaNext = (lambdaU + lambdaL)/2;
if lambdaUfound==true
    while(lambdaU - lambdaL > epsilon)
    Constraints = [];
    Constraints = Constraints + (P_var >= eye(n)*0.000005);
    lambda =  lambdaNext;
    for i=1:N
        Constraints = Constraints + (Y{i}'*P_var*Y{i} <= lambda^2*X{i}'*P_var*X{i});
    end
    sol = optimize(Constraints,Objective,ops);
    if sol.problem==0
        feasible=true;
        lambdaU=lambda;
        lambdaNext=(lambdaU+lambdaL)/2;
    else
        feasible=false;
        lambdaL=lambda;
        lambdaNext=(lambdaU+lambdaL)/2;
    end
    end
if feasible==false
    lambda=lambdaU;
end

lowerBound{N/20}=lambda/sqrt(n);
fprintf('The lower bound is %4f',lowerBound{N/20})
end
end

lb=[];
for j=1:50
 lb=[lb lowerBound{j}];    
end

lb


