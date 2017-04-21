function plotRhoForRandomSystems(m,n, beta, Nmax)
close all;
%this function takes m and n and plots the evaluation of rho (and delta)
%with N for a random system
addpath(genpath('YALMIP'));
addpath(genpath('JSR'));
addpath(genpath('sdpt3'));

increment = 120;
for i = 1:m
    sys = drss(n);
    A{i} = sys.a;
end

jsrRaphael = jsr_prod_bruteForce(A);
jsrMean = mean(jsrRaphael);

d = (n*(n+1))/2;
counter  = 1;
for i = d+1 : increment : Nmax    
    [lowerBound(counter), upperBound(counter)] = computeRhoBlackbox(A, beta, i);
    counter = counter +1;
end

indexArray = d+1 : increment : Nmax;
figure;
plot(indexArray, lowerBound, '-r*');
hold on;
plot(indexArray, upperBound, '-b*');
hold on;
plot(indexArray, jsrMean*ones(length(indexArray)),'k');

end