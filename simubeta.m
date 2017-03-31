%clear all;
%close all;
%clc;
addpath(genpath('YALMIP'));
addpath(genpath('JSR'));
addpath(genpath('sdpt3'));

%options for the solver
ops = sdpsettings('solver','sdpt3','verbose','0');

beta=0.92; %certainty required by the user
counterOK = 0; %number of times we have an upper bound for the JSR

Nsystems = 500; %number of tests: each test we take a different system and a different sample

for i = 1:Nsystems
    % dimension of the state space, between 2 and 7
    n = 1+unidrnd(6);
    %number of modes, between 2 and 5
    m = 1+unidrnd(4);
    A = cell(1,m); %cell A stores the modes of the system
    d = n*(n+1)/2+1;

    %number of points in the sampling, between d+1 and 1000
    N = d + unidrnd(800);
    X = cell(1,N);
    Y = cell(1,N);

    epsilon=1-betaincinv(1-beta,N-d,d+1);

    l=unidrnd(2); %if l=1, we have rho < 1, if l=2, we can have rho>1

    %generate the m modes
    for j = 1:m
        v=rand(1,n);
        D = diag(0.85*l*v); % random eigenvalues of the mode i 
        V = orth(randn(n)); 
        E = V*D*V';
        A{j} = E;
    end

    c=jsr_prod_bruteForce(A); %computation of the JSR with brute force method provided by the JSR Toolbox

    for j = 1:N %generate uniformly points of the unit sphere
        v=randn(n,1);
        X{j}=v/sqrt(sum(v.^2));
    end

    %apply randomly one mode to every point, store the result in array Y
    for j = 1:N
        k = unidrnd(m); %pick uniformely random index of the mode apply to the sampled point X{i}
        Y{j} = A{k}*X{j}; %output obtained
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
        Constraints = Constraints + (P_var >= eye(n));
        for j=1:N
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

    feasible = true;
    lambdaNext = (lambdaU + lambdaL)/2;
    if lambdaUfound==true
        while((lambdaU - lambdaL) > 0.005)
            Constraints = [];
            Constraints = Constraints + (P_var >= eye(n));
            lambda =  lambdaNext;
                for j=1:N
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
    end

    if feasible==false
        lambda=lambdaU;
    end

    gammaStar=lambda; %this is the value of gamma*(\omega_N) we will use for the upper bound

    %find P
    Constraints = [];
    Constraints = Constraints + (P_var >= eye(n));
    Objective=lambda_max(P_var);
        for j=1:N
            Constraints = Constraints + (Y{j}'*P_var*Y{j} <= gammaStar^2*X{j}'*P_var*X{j});
        end
    sol = optimize(Constraints,Objective,ops);
        if sol.problem~=0
            display('Feasibility problem')
        end

    P=value(P_var);
    
    lambdaMax=max(eig(P));
    dP=det(P);
    epsilon1 = min(1,m*epsilon/2*sqrt(lambdaMax^n/dP));


    alpha=betaincinv(epsilon1,(n-1)/2,1/2);

    delta=sqrt(1-alpha);
    upperBound=gammaStar/delta;

    if (upperBound >=c(2))
        counterOK=counterOK+1;
    end
end
proba = counterOK/Nsystems;
fprintf('%f\n',proba)