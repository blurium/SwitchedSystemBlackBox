%clear all;
%close all;
%clc;
addpath(genpath('YALMIP'));
addpath(genpath('JSR'));
addpath(genpath('sdpt3'));

beta=0.92;

clear A;
clear X;
clear Y;
clear n;
clear m;
clear N;
% dimension of the state space
%n = 4;
n=1+unidrnd(6);
%number of modes
%m = 5; 
m=1+unidrnd(4);
%number of points in the sampling
%N = 400;
N= 30 + unidrnd(800);

d=n*(n+1)/2+1;
epsilon=1-betaincinv(1-beta,N-d,d+1);

%options for the solver
ops = sdpsettings('solver','sdpt3','verbose','0');

l=unidrnd(2);

%generate the m modes
for i = 1:m
    v=rand(1,n);
    D = diag(0.85*l*v); % Random eigenvalues almost in the unit circle
    V = orth(randn(n)); 
    E = V*D*V';
    A{i} = E;
end

c=jsr_prod_bruteForce(A);

for i=1:N %generate uniformly points of the unit sphere
      v=randn(n,1);
      X{i}=v/sqrt(sum(v.^2));
end

%apply randomly one mode to every point, store the result in array Y
for i=1:N
    k=unidrnd(m);
    Y{i}=A{k}*X{i};
end


lambdaNext = 1;
lambdU = 1;
lambdaL = 0;


P_var = sdpvar(n,n); %optimization variable
Objective = 1;
counter = 0;
lambdaUfound = false;

while(counter<50 && lambdaUfound==false)
    Constraints = [];
    lambda = lambdaNext;
    Constraints = Constraints + (P_var >= 0.00005*eye(n));
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
    while((lambdaU - lambdaL) > 0.005)
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
end

if feasible==false
    lambda=lambdaU;
end

gammaStar=lambda; %this is the value of gamma*(\omega_N) we will use for the upper bound

%find P
Constraints = [];
Constraints = Constraints + (P_var >= eye(n));
Objective=lambda_max(P_var);
    for i=1:N
        Constraints = Constraints + (Y{i}'*P_var*Y{i} <= gammaStar^2*X{i}'*P_var*X{i});
    end
    sol = optimize(Constraints,Objective,ops);
if sol.problem~=0
    display('Problem with gammaStar')
end

P=value(P_var);
    
    
lambdaMax=max(eig(P));
dP=det(P);
epsilon1 = min(1,m*epsilon/2*sqrt(lambdaMax^n/dP));

alpha1=min(1,epsilon1*gamma(d/2)/(pi^(d/2)));
alpha=betaincinv(alpha1,(d-1)/2,1/2);

delta=sqrt(1-alpha);
upperBound=gammaStar/delta;

if (upperBound >=c(2))
   counterOK=counterOK+1;
end