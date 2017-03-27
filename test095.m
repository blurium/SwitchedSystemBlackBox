%example with rho = 0.958096329922833, 5 modes, n = 4
%clear all;
%close all;
%clc;
addpath(genpath('YALMIP'));
addpath(genpath('JSR'));
addpath(genpath('sdpt3'));

%options for the solver
ops = sdpsettings('solver','sdpt3');

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
%Campi's certainty
beta = 0.92;
d=n*(n+1)/2+1; %t minized 

deltaarray=zeros(1,26);
lbarray=zeros(1,26);
ubarray=zeros(1,26);
ini = 15;
inc = 200;
final = ini+10*inc; %we osbserve 11 samplings

for i=1:ini %generate uniformly 20 more random points of the unit sphere
      v=randn(n,1);
      X{i}=v/sqrt(sum(v.^2));
end

%apply randomly one mode to every point, store the result in array Y
for i=1:ini
    k=unidrnd(m);
    Y{i}=A{k}*X{i};
end

epsilon=1-betaincinv(1-beta,ini-d,d+1);


for i=1:ini
   W{i}=X{i};
   W{i+ini}=-X{i};
end

for i=1:ini
   Z{i}=Y{i};
   Z{i+ini}=-Y{i};
end
%computation of gamma* for the opt problem of the lower bound by bisection, we start with lambda = 1, and we stop after 50 iterations
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
    Constraints = Constraints + (P_var >= 0.000005*eye(n));
    for i=1:(2*ini)
    Constraints = Constraints + (Z{i}'*P_var*Z{i} <= lambda^2*W{i}'*P_var*W{i});
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
    while(lambdaU - lambdaL > 0.01)
    Constraints = [];
    Constraints = Constraints + (P_var >= eye(n));
    lambda =  lambdaNext;
    for i=1:(2*ini)
        Constraints = Constraints + (Z{i}'*P_var*Z{i} <= lambda^2*W{i}'*P_var*W{i});
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
lbarray(1)=lowerBound;
end

lambdaNext = 1;
lambdU = 1;
lambdaL = 0;
counter = 0;
lambdaUfound = false;
while(counter<50 && lambdaUfound==false)
    Constraints = [];
    lambda = lambdaNext;
    Constraints = Constraints + (P_var >= eye(n));
    for i=1:ini
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
    while((lambdaU - lambdaL) > 0.01)
    Constraints = [];
    Constraints = Constraints + (P_var >= eye(n));
    lambda =  lambdaNext;
        for i=1:ini
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
    for i=1:ini
        Constraints = Constraints + (Y{i}'*P_var*Y{i} <= gammaStar^2*X{i}'*P_var*X{i});
    end
    sol = optimize(Constraints,Objective,ops);
if sol.problem~=0
    display('Problem with gammaStar')
end

P=value(P_var);
    
    
lambdaMax=max(eig(P));
dP=det(P);
epsilon1 = min(1,m*(epsilon/2)*sqrt((lambdaMax^n)/dP));

alpha1=min(1,epsilon1);
alpha=betaincinv(alpha1,(n-1)/2,1/2);

delta=sqrt(1-alpha);
upperBound=gammaStar/delta;
deltaarray(1)=delta;
ubarray(1)=upperBound;



%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%
%loop for all other cases

for N=ini:inc:(final-inc)

%epsilon as function of beta and N
%epsilon = 1 - I^{-1}(beta, N-d,d+1)
epsilon=1-betaincinv(1-beta,N-d,d+1);

for i=(N+1):(N+inc) %generate uniformly 20 more random points of the unit sphere
      v=randn(n,1);
      X{i}=v/sqrt(sum(v.^2));
end

%apply randomly one mode to every point, store the result in array Y
for i=(N+1):(N+inc)
    k=unidrnd(m);
    Y{i}=A{k}*X{i};
end

%we use the points W and Z{i} for the lower bound which are the points sampled
%and their opposite, W is the set of points X extended, Z is the set of
%points Y (the images) extended

for i=1:(N+inc)
   W{i}=X{i};
   W{i+N+inc}=-X{i};
end

for i=1:(N+inc)
   Z{i}=Y{i};
   Z{i+N+inc}=-Y{i};
end
%computation of gamma* for the opt problem of the lower bound by bisection, we start with lambda = 1, and we stop after 50 iterations
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
    for i=1:(2*(N+inc))
    Constraints = Constraints + (Z{i}'*P_var*Z{i} <= lambda^2*W{i}'*P_var*W{i});
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
    while(lambdaU - lambdaL > 0.01)
    Constraints = [];
    Constraints = Constraints + (P_var >= eye(n));
    lambda =  lambdaNext;
    for i=1:(2*(N+inc))
        Constraints = Constraints + (Z{i}'*P_var*Z{i} <= lambda^2*W{i}'*P_var*W{i});
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
lbarray(1+((N+inc-ini)/inc))=lowerBound;
end
%%%%%%%%%%%%%%%%%%%%%%%UPPER BOUND%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%we compute the gamma* for the problem of the upper bound
lambdaNext = 1;
lambdU = 1;
lambdaL = 0;
counter = 0;
lambdaUfound = false;
while(counter<50 && lambdaUfound==false)
    Constraints = [];
    lambda = lambdaNext;
    Constraints = Constraints + (P_var >= eye(n));
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
    while((lambdaU - lambdaL) > 0.01)
    Constraints = [];
    Constraints = Constraints + (P_var >= eye(n));
    lambda =  lambdaNext;
        for i=1:(N+inc)
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
    for i=1:(N+inc)
        Constraints = Constraints + (Y{i}'*P_var*Y{i} <= gammaStar^2*X{i}'*P_var*X{i});
    end
    sol = optimize(Constraints,Objective,ops);
if sol.problem~=0
    display('Problem with gammaStar')
end

P=value(P_var);
    
    
lambdaMax=max(eig(P));
dP=det(P);
epsilon1 = min(1,m*(epsilon/2)*sqrt(lambdaMax^n/dP));

alpha1=min(1,epsilon1);
alpha=betaincinv(alpha1,(n-1)/2,1/2);

delta=sqrt(1-alpha);
upperBound=gammaStar/delta;
deltaarray(1+((N+inc-ini)/inc))=delta;
ubarray(1+((N+inc-ini)/inc))=upperBound;

end