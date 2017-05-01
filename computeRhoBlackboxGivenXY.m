function [lowerBound, upperBound] = computeRhoBlackboxGivenXY(A, beta, X, Y)
%inputs : A matrices of switched system
%beta : desired confidence factor
%X = cell(1,N); %cell where sampled vectors of the unit sphere are stored

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

[gammaStar, P, lowerBound]=computePandGamma(X, Y);

%computation of delta and the upper bound
lambdaMax=max(eig(P));
dP=det(P);
epsilon1 = min(1/2,m*(epsilon/2)*sqrt(lambdaMax^n/dP));;

alpha=betaincinv(2*epsilon1,(n-1)/2,1/2);

delta=sqrt(1-alpha);
upperBound=gammaStar/delta;

end