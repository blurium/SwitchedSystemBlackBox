function [lowerBound, upperBound] = computeRhoBlackboxGivenXY(A, beta, X, Y)
%inputs : A matrices of switched system
%beta : desired confidence factor
%X = cell(1,N); %cell where sampled vectors of the unit sphere are stored
ops = sdpsettings('solver','sdpt3','verbose','0');

m = length(A); % number of modes
n = size(A{1},1); % number of states
d = n*(n+1)/2+1;
N = length(X);
%epsilon as function of beta and N ; epsilon = 1 - I^{-1}(beta, N-d,d+1)
epsilon=1-betaincinv(1-beta,N-d,d+1);
if (epsilon > 2 / m)
    lowerBound = 0;
    upperBound = inf;
    return;
end
%computation of gamma* for the optimization problem by bisection, we start
%by looking for a valid upper bound, lambdaU, of gamma* before starting the bisection
%itself

lambdaL = 0;

P_var = sdpvar(n,n); %optimization variable
Objective = 1;
flag = 0;
lambdaVar = 0;

for j=1:N
   lambdaVar = max((Y{j}'*Y{j})/(X{j}'*X{j}), lambdaVar);
end
%we proceed to the bisection to find gamma* (we get in fact an
%overapproximation of gamma*, ensuring feasibility of the next optimization
%problem)
lambdaU = lambdaVar;
%we proceed to the bisection to find gamma* (we get in fact an
%overapproximation of gamma*, ensuring feasibility of the next optimization
%problem)
lambdaNext = (lambdaU + lambdaL)/2;
feasibleLast = false;

while(lambdaU - lambdaL > 0.1 || feasibleLast ~= true)
    Constraints = [];
    Constraints = Constraints + (P_var >= eye(n));
    lambda =  lambdaNext;
    for j=1:N
        Constraints = Constraints + (Y{j}'*P_var*Y{j} <= lambda^2*X{j}'*P_var*X{j});
    end
    sol = optimize(Constraints,Objective,ops);
    if sol.problem==0
        flag = 1;
        lambdaU=lambda;
        lambdaNext=(lambdaU+lambdaL)/2;
        feasibleLast = true;
    else
        if (flag == 0)
            lambdaU = 2*lambdaU;
            lambdaNext=(lambdaU+lambdaL)/2;
            feasibleLast = false;
        else
            lambdaL=lambda;
            lambdaNext=(lambdaU+lambdaL)/2;
            lambda = lambdaU;
            feasibleLast = false;
        end
        
    end    
end

if (flag == 0)
    error('couldnt compute gammastar\n');
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
    Constraints = [];
    Constraints = Constraints + (P_var >= eye(n));
    for j=1:N
        Constraints = Constraints + (Y{j}'*P_var*Y{j} <= (gammaStar+0.1)^2*X{j}'*P_var*X{j});
    end
    sol = optimize(Constraints,Objective,ops);
    if sol.problem~=0
        disp('Feasibility problem')
    end
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