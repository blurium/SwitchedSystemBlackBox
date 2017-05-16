%example with rho = 0.792888603593367, 5 modes, n = 4
clear all;
close all;
clc;
rng(2011);
addpath(genpath('YALMIP'));
addpath(genpath('JSR'));
addpath(genpath('sdpt3'));

%options for the solver
ops = sdpsettings('solver','sdpt3');

% dimension of the state space
n = 2;

%number of modes
m = 4; 
%Campi's certainty
beta = 0.95;
N = 500;
increment = 25;
maxJSR = 0.99;
[A, jsrRaphael] = createRandomSystem(n,m,maxJSR);
jsr_prod_pruningAlgorithm(A)

return;

for i = 1:5
    [lowerBound{i}, upperBound{i}, delta{i}] = computeRhoSequentialSampling(A, beta, N,increment);
end

for i = 1:length(lowerBound{1})   
    lowerBoundMean(i) = mean(cellfun(@(v) v(i), lowerBound));
    upperBoundMean(i) = mean(cellfun(@(v) v(i), upperBound));
    deltaMean(i) = mean(cellfun(@(v) v(i), delta));
end
%%
%save paper2D
%%
close all; 
d = n*(n+1)/2+1;
figure;
fig1 = gcf;
plot(d*2:increment:N, lowerBoundMean,'LineWidth',1.5);
hold on;
plot(d*2:increment:N, upperBoundMean,'LineWidth',1.5);
hold on;
plot(d*2:increment:N, mean(jsrRaphael)*ones(length(d*2:increment:N),1),'k-.','Linewidth',0.75);
hold on;
plot(d*2:increment:N, (mean(jsrRaphael)/sqrt(n))*ones(length(d*2:increment:N),1), 'k-.','Linewidth', 0.75);
legend('Lower bound','Upper bound','\rho','\rho / \surd n','Interpreter','latex');
xlabel('Number of samples (N)')

ylim([0 2.1])
grid on;
set(gcf,'paperunits','centimeters','papersize',[15 10],'paperposition',[0 0 15 10])
print -dpdf experiment1


figure;
plot(d*2:increment:N, deltaMean,'LineWidth',1.5)
grid on;
set(gcf,'paperunits','centimeters','papersize',[15 10],'paperposition',[0 0 15 10])
xlabel('Number of samples (N)')
ylabel('$\delta(0.95,\omega_N)$','Interpreter','latex');
print -dpdf delta1