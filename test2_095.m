%example with rho = 0.958096329922833, 5 modes, n = 4
%clear all;
%close all;
%clc;
addpath(genpath('YALMIP'));
addpath(genpath('sdpt3'));

%options for the solver
ops = sdpsettings('solver','sdpt3','verbose','0');

% dimension of the state space
n = 4;

%number of modes
m = 5; 
A = cell(1,m); %A stores the modes of the system

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


beta = 0.92; %certainty required by the user
d = n*(n+1)/2+1;

ini = d+1; %size of the first sampling, taken as the minimum possible size
inc = 1000; %number of points added to the sample at each iteration
iter = 10; %number of iterations after the initial sample
final = ini+iter*inc; %size of the last sample considered

dArray = zeros(1,iter+1); %array of values of delta computed
lbArray = zeros(1,iter+1); %array of values of the lower bound computed
ubArray = zeros(1,iter+1); %array of values of the upper bound computed

Nsamples=[0 ini:inc:final]; %cardinality of the samples, we start with empty set of cardinality zero
X = cell(1,final); %cell where sampled vectors of the unit sphere are stored
Y = cell(1,final); %cell where output vectors are stored

for i=1:length(Nsamples)

N=Nsamples(i+1); %size of the current sample

%epsilon as function of beta and N ; epsilon = 1 - I^{-1}(beta, N-d,d+1)
epsilon=1-betaincinv(1-beta,N-d,d+1);

for j=(Nsamples(i)+1):Nsamples(i+1) %generate uniformly random points of the unit sphere
      v=randn(n,1);
      X{j}=v/sqrt(sum(v.^2));
end

%apply randomly one mode to every point, store the result in array Y
for j=(Nsamples(i)+1):Nsamples(i+1)
    k=unidrnd(m); %random uniform generation of the index of the mode applied to sampled point X{j}
    Y{j}=A{k}*X{j};
end


%computation of gamma* for the optimization problem by bisection, we start
%by looking for a valid upper bound, lambdaU, of gamma* before starting the bisection
%itself

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
    Constraints = Constraints + (P_var >= eye(n));
        for j=1:Nsamples(i+1)
            Constraints = Constraints + (Y{j}'*P_var*Y{j} <= lambda^2*X{j}'*P_var*X{j});
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
    display('rho bigger!');
end

if lambdaUfound==false
    display('rho too high!')
end

%we proceed to the bisection to find gamma* (we get in fact an
%overapproximation of gamma*, ensuring feasibility of the next optimization
%problem)
feasible = true;
lambdaNext = (lambdaU + lambdaL)/2;

if lambdaUfound==true

    while(lambdaU - lambdaL > 0.005)
        Constraints = [];
        Constraints = Constraints + (P_var >= eye(n));
        lambda =  lambdaNext;
        for j=1:Nsamples(i+1)
             Constraints = Constraints + (Y{j}'*P_var*Y{j} <= lambda^2*X{j}'*P_var*X{j});
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
lbArray(i)=lowerBound;

gammaStar=lambda; %this is the value of gamma*(\omega_N) we will use for to find P, delta and the upper bound

%find P
Constraints = [];
Constraints = Constraints + (P_var >= eye(n));
Objective=lambda_max(P_var); %we want to minize lambda_max(P), while keeping lambda_min(P)>=1
    for j=1:Nsamples(i+1)
        Constraints = Constraints + (Y{j}'*P_var*Y{j} <= gammaStar^2*X{j}'*P_var*X{j});
    end
    sol = optimize(Constraints,Objective,ops);
if sol.problem~=0
    display('Feasibility problem')
end

P=value(P_var);
    
%computation of delta and the upper bound    
lambdaMax=max(eig(P));
dP=det(P);
epsilon1 = min(1,m*(epsilon/2)*sqrt(lambdaMax^n/dP));

alpha=betaincinv(epsilon1,(n-1)/2,1/2);

delta=sqrt(1-alpha);
upperBound=gammaStar/delta;
dArray(i)=delta;
ubArray(i)=upperBound;
end
end