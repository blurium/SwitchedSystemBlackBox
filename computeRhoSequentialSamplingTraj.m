function [lowerBound, upperBound, delta] = computeRhoSequentialSamplingTraj(A, beta, N, increment, l)
m = length(A); % number of modes
n = size(A{1},1); % number of states
d = n*(n+1)/2+1;

X = cell(1,N); %cell where sampled vectors of the unit sphere are stored
Y = cell(1,N); %cell where output vectors are stored


%epsilon as function of beta and N ; epsilon = 1 - I^{-1}(beta, N-d,d+1)
epsilon=1-betaincinv(1-beta,N-d,d+1);
if (epsilon > 2 / m)
    lowerBound = 0;
    upperBound = inf;
    return;
end


for j=1:N %generate uniformly random points of the unit sphere
    v=randn(n,1);
    X{j}=v/sqrt(sum(v.^2));
    k=unidrnd(m); %random uniform generation of the index of the mode applied to sampled point X{j}
    Y{j}=(A{k}^l)*X{j};
end
counter = 1;
for i = d*2:increment:N
    [lowerBound(counter), upperBound(counter), delta(counter)] = computeRhoBlackboxGivenXYTraj(A, beta, X(1:i), Y(1:i), l);
    counter = counter +1;
end

% close all;
% c=jsr_prod_bruteForce(A); 
% plot(d*2:increment:N, lowerBound);
% hold on;
% plot(d*2:increment:N, upperBound);
% hold on;
% plot(d*2:increment:N, c(2)*ones(length(d*2:increment:N),1),'k-.','Linewidth',0.75);
% hold on;
% plot(d*2:increment:N, (c(2)/sqrt(n))*ones(length(d*2:increment:N),1), 'k-.','Linewidth', 0.75);
% legend('Lower bound','Upper bound','\rho','\rho / \surd n','Interpreter','latex');
% xlabel('Number of samples (N)')
% 
% 
% grid on;
% set(gcf,'paperunits','centimeters','papersize',[15 10],'paperposition',[0 0 15 10])
% print -dpdf experiment2
% 
% figure;
% plot(d*2:increment:N, delta)
% grid on;
% set(gcf,'paperunits','centimeters','papersize',[15 10],'paperposition',[0 0 15 10])
% xlabel('Number of samples (N)')
% ylabel('$\delta(0.92,\omega_N)$','Interpreter','latex');
% print -dpdf delta2


end