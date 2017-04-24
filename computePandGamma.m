function [gammaStar, P]=computePandGamma(X, Y)
lambdaL = 0;
ops = sdpsettings('solver','sdpt3','verbose','0');

Objective = 1;
flag = 0;
lambdaVar = 0;
N = length(X);
n = length(X{1});
P_var = sdpvar(n,n); %optimization variable

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
end
