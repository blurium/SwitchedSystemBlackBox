%example with rho = 0.792888603593367, 5 modes, n = 4
clear all;
close all;
clc;
%rng(2511);
addpath(genpath('YALMIP'));
addpath(genpath('JSR'));
addpath(genpath('sdpt3'));

%options for the solver
ops = sdpsettings('solver','sdpt3');

% dimension of the state space
n = 2;
l = 2;
noOfSystems = 1;
lMax = 6;
%number of modes
m = 4; 
%Campi's certainty
beta = 0.95;
N = 750;
increment = 50;
maxJSR = 1.2;
[A, jsrRaphael] = createRandomSystem(n,m,maxJSR);
jsr_prod_pruningAlgorithm(A)


for i = 1:lMax
    [lowerBound{i}, upperBound{i}, delta{i}] = computeRhoSequentialSamplingTraj(A, beta, N,increment,i);
end

% for i = 1:length(lowerBound{1})   
%     lowerBoundMean(i) = mean(cellfun(@(v) v(i), lowerBound));
%     upperBoundMean(i) = mean(cellfun(@(v) v(i), upperBound));
%     deltaMean(i) = mean(cellfun(@(v) v(i), delta));
% end
%%
%save paper2D
%%
close all; 
jsrRaphael = jsr_prod_bruteForce(A);
d = n*(n+1)/2+1;
figure;
fig1 = gcf;
for i = 1:5
    plot(d*2:increment:N, upperBound{i},'LineWidth',1.5);
    hold on;

    legend('l=1','l=2','l=3','l=4','l=5');
end
plot(d*2:increment:N, jsrRaphael(1)*ones(length(d*2:increment:N),1),'k-.','Linewidth',0.75);
hold on;
%plot(d*2:increment:N, (mean(jsrRaphael)/n^(1/(2*l)))*ones(length(d*2:increment:N),1), 'k-.','Linewidth', 0.75);
%legend('Lower bound','Upper bound','\rho','\rho / \surd n','Interpreter','latex');
xlabel('Number of samples (N)')
% 
% ylim([0 2.1])
% grid on;
% set(gcf,'paperunits','centimeters','papersize',[15 10],'paperposition',[0 0 15 10])
% print -dpdf experiment1
% 
figure;
for i = 1:5
plot(d*2:increment:N, lowerBound{i},'LineWidth',1.5);
hold on;
plot(d*2:increment:N, (jsrRaphael(2)/n^(1/(2*i)))*ones(length(d*2:increment:N),1), '-.','Linewidth', 0.75);
hold on;
  legend('l=1','l=1 (theoretical)','l=2','l=2 (theoretical))','l=3',...
      'l=3 (theoretical)','l=4','l=4 (theoretical)','l=5','l=5 (theoretical)')
end

figure;
for i = 1:5
plot(d*2:increment:N, lowerBound{i},'LineWidth',1.5);
hold on;
plot(d*2:increment:N, (jsrRaphael(2)/n^(1/(2*i)))*ones(length(d*2:increment:N),1), '-.','Linewidth', 0.75);
hold on;
  legend('l=1','l=1 (theoretical)','l=2','l=2 (theoretical))','l=3',...
      'l=3 (theoretical)','l=4','l=4 (theoretical)','l=5','l=5 (theoretical)')
end


%%
figure;
plot(d*2:increment:N, lowerBound{1},'r','LineWidth',1.5);
hold on;
plot(d*2:increment:N, lowerBound{2},'g','LineWidth',1.5);
hold on;
plot(d*2:increment:N, lowerBound{3},'b','LineWidth',1.5);
hold on;
plot(d*2:increment:N, upperBound{1},'r','LineWidth',1.5);
hold on;
plot(d*2:increment:N, upperBound{2},'g','LineWidth',1.5);
hold on;
plot(d*2:increment:N, upperBound{3},'b','LineWidth',1.5);
hold on;
plot(d*2:increment:N, (jsrRaphael(1))*ones(length(d*2:increment:N),1), 'k-.','Linewidth', 0.75);
legend('l=1','l=2','l=3');
title('n=2, m=4');
grid on;
%print -dpdf delta1