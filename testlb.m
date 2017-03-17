%clear all;
%close all;
%clc;
addpath(genpath('YALMIP'));
addpath(genpath('JSR'));
addpath(genpath('sdpt3'));

clear A;
clear X;
clear Y;
% dimension of the state space
%n = 4;
n=1+unidrnd(6);
%number of modes
%m = 5; 
m=1+unidrnd(4);
%number of points in the sampling
%N = 400;
N= 30 + unidrnd(800);
%bounds to generate the sampled points x_i
uBound = 5;
lBound = -5;

%options for the solver
ops = sdpsettings('solver','sdpt3');

l=unidrnd(2);

%generate the m modes
for i = 1:m
    v=rand(1,n);
    D = diag(0.85*l*v); % Random eigenvalues almost in the unit circle
    V = orth(randn(n)); 
    E = V*D*V';
    A{i} = E;
end

%generate the N points
for i=1:N
    X{i} = lBound + (uBound-lBound)*rand(1,n);
    X{i}=X{i}';
    X{N+i}=-X{i};
end

%apply randomly one mode to every point, store the result in array Y
for i=1:N
    k=unidrnd(m);
    Y{i}=A{k}*X{i};
    Y{N+i}=-Y{i};
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
    Constraints = Constraints + (P_var >= 0.000005*eye(n));
    for i=1:(2*N)
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
    for i=1:(2*N)
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

lowerBound=lambda/sqrt(n);
fprintf('The lower bound is %4f',lowerBound)

a=max(rho(A));
b=jsr_lift_semidefinite(A);
c=jsr_prod_bruteForce(A);
end