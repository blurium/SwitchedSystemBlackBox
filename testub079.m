%example with rho = 0.792888603593367, 5 modes, n = 4
%clear all;
%close all;
%clc;
addpath(genpath('YALMIP'));
addpath(genpath('JSR'));
addpath(genpath('sdpt3'));

%options for the solver
ops = sdpsettings('solver','sdpt3','verbose','0');

% dimension of the state space
n = 4;

%number of modes
m = 5; 

A{1}=[0.3908   -0.1289    0.1059   -0.0352;
   -0.1289    0.4131    0.1630    0.0636;
    0.1059    0.1630    0.1677    0.0249;
   -0.0352    0.0636    0.0249    0.0226];

A{2}=[0.5015   -0.0902    0.0964   -0.1331;
   -0.0902    0.6505   -0.0615   -0.1678;
    0.0964   -0.0615    0.6865   -0.0510;
   -0.1331   -0.1678   -0.0510    0.5698];

A{3}=[0.4777    0.1141    0.1424   -0.2475
    0.1141    0.6455   -0.0432    0.0335
    0.1424   -0.0432    0.5784    0.1376
   -0.2475    0.0335    0.1376    0.4474];

A{4}=[0.3340   -0.1382    0.2024   -0.1523;
   -0.1382    0.5044    0.0449    0.0595;
    0.2024    0.0449    0.5062   -0.1514;
   -0.1523    0.0595   -0.1514    0.4033];

A{5}=[0.7129   -0.0105    0.0018    0.0353;
   -0.0105    0.6289   -0.0526   -0.1887;
    0.0018   -0.0526    0.7330   -0.0293;
    0.0353   -0.1887   -0.0293    0.4932];


c=jsr_prod_bruteForce(A);

%Campi's certainty
beta = 0.92;
d=n*(n+1)/2+1; %t minized 

deltaarray=zeros(1,11);
lbarray=zeros(1,11);
ubarray=zeros(1,11);

for N=5000:5000:10000

%epsilon as function of beta and N
%epsilon = 1 - I^{-1}(beta, N-d,d+1)
epsilon=1-betaincinv(1-beta,N-d,d+1);

for i=(N-499):N %generate uniformly 20 more random points of the unit sphere
      v=randn(n,1);
      X{i}=v/sqrt(sum(v.^2));
end

%apply randomly one mode to every point, store the result in array Y
for i=(N-499):N
    k=unidrnd(m);
    Y{i}=A{k}*X{i};
end

%we use the points W and Z{i} for the lower bound which are the points sampled
%and their opposite, W is the set of points X extended, Z is the set of
%points Y (the images) extended

for i=1:N
   W{i}=X{i};
   W{i+N}=-X{i};
end

for i=1:N
   Z{i}=Y{i};
   Z{i+N}=-Y{i};
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
    for i=1:(2*N)
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
    while(lambdaU - lambdaL > 0.005)
    Constraints = [];
    Constraints = Constraints + (P_var >= eye(n)*0.000005);
    lambda =  lambdaNext;
    for i=1:(2*N)
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
lbarray((N)/500)=lowerBound;
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
epsilon1 = min(1,m*(epsilon/2)*sqrt(lambdaMax^n/dP));

alpha1=min(1,epsilon1*gamma(d/2)/(pi^(d/2)));
alpha=betaincinv(alpha1,(d-1)/2,1/2);

delta=sqrt(1-alpha);
upperBound=gammaStar/delta;
deltaarray(N/500)=delta;
ubarray(N/500)=upperBound;

end
display('FINI!!!!')