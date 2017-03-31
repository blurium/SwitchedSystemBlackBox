function [lowerBound, upperBound] = computeRhoBlackboxGivenX(A, beta, X)
%inputs : A matrices of switched system
%beta : desired confidence factor
%X = cell(1,N); %cell where sampled vectors of the unit sphere are stored
ops = sdpsettings('solver','sdpt3','verbose','0');

m = length(A); % number of modes
n = size(A{1},1); % number of states
d = n*(n+1)/2+1;
N = length(X);

Y = cell(1,N); %cell where output vectors are stored


%epsilon as function of beta and N ; epsilon = 1 - I^{-1}(beta, N-d,d+1)
epsilon=1-betaincinv(1-beta,N-d,d+1);

for j=1:N %generate uniformly random points of the unit sphere
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
    if feasible==false
        lambda=lambdaU;
    end
    
    
    lowerBound=lambda/sqrt(n);
    
    gammaStar=lambda; %this is the value of gamma*(\omega_N) we will use for to find P, delta and the upper bound
    
    %find P
    Constraints = [];
    Constraints = Constraints + (P_var >= eye(n));
    Objective=lambda_max(P_var); %we want to minize lambda_max(P), while keeping lambda_min(P)>=1
    for j=1:N
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
end

end