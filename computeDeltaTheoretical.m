function [delta, N]=  computeDeltaTheoretical(n,m,beta,Nmax)
%this function takes n,m and beta and computes delta(N)
%from n(n+1)/2 to Nmax
d = n*(n+1)/2+1;
counter = 1;
for i = (n*(n+1))/2+1:1:Nmax
    epsilon = 1-betaincinv(1-beta,i-d,d+1);
    epsilon1 = min(1/2,m*(epsilon/2));  
    epsilon2 = min(1/2,(1-(1-epsilon*m)));
    epsilonMin = min(epsilon1, epsilon2);
    delta(counter) = sqrt(1-betaincinv(2*epsilonMin,(n-1)/2,1/2));
    N(counter) = i;
    counter = counter +1;
end





end