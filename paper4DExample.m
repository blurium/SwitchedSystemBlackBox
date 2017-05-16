%example with rho = 0.792888603593367, 5 modes, n = 4
%clear all;
%close all;
%clc;
close all;
clear;
rng(15);
addpath(genpath('YALMIP'));
addpath(genpath('JSR'));
addpath(genpath('sdpt3'));

%options for the solver
ops = sdpsettings('solver','sdpt3');

% dimension of the state space
n = 4;

%number of modes
m = 6; 


% ans =
% 
%     0.0550    0.2428   -0.2142   -0.3082
%     0.0971    0.2591    0.2548   -0.3065
%     0.0677    0.2342    0.0190    0.1161
%    -0.2873   -0.2054   -0.1838    0.1520
% 
% 
% ans =
% 
%     0.1436   -0.0645   -0.0147    0.0869
%    -0.0317    0.3489    0.2238    0.2809
%     0.0071    0.0623    0.3268    0.0205
%     0.1466    0.1293   -0.0046    0.0408
% 
% 
% ans =
% 
%    -0.1546   -0.0097    0.3450   -0.4638
%    -0.3075    0.0920   -0.4186    0.0632
%    -0.0141   -0.3167   -0.0505    0.1554
%    -0.3382    0.0483    0.0788    0.1706
% 
% 
% ans =
% 
%     0.1341    0.2431   -0.3785   -0.0886
%     0.0945   -0.0158   -0.2534    0.1543
%    -0.1038   -0.0470    0.2680   -0.1551
%     0.1920    0.1127    0.0155   -0.1870
% 
% 
% ans =
% 
%    -0.5914    0.4561    0.2015    0.1571
%     0.0692   -0.6361    0.1248    0.2781
%     0.2956    0.2182   -0.5056   -0.3107
%     0.2353    0.1327   -0.2034    0.6679
% 
% 
% ans =
% 
%    -0.5714   -0.0094    0.0349   -0.1995
%    -0.1435    0.0344    0.0185    0.4903
%     0.1328    0.0619   -0.3813    0.0460
%    -0.0280    0.4745    0.1192   -0.3247
maxJSR = 0.9;
[A, jsrRaphael] = createRandomSystem(n,m,maxJSR);

%Campi's certainty
beta = 0.95;
N = 5000;
increment = 200;
d = (n*(n+1))/2 +1;
for i = 1:5
    [lowerBound{i}, upperBound{i}, delta{i}] = computeRhoSequentialSampling(A, beta, N,increment);
end

for i = 1:length(lowerBound{1})   
    lowerBoundMean(i) = mean(cellfun(@(v) v(i), lowerBound));
    upperBoundMean(i) = mean(cellfun(@(v) v(i), upperBound));
    deltaMean(i) = mean(cellfun(@(v) v(i), delta));
end

%%
close all;
%save experiment4D

figure;
fig1 = gcf;
c = jsrRaphael;
plot(d*2:increment:N, lowerBoundMean,'LineWidth',1.5);
hold on;
plot(d*2:increment:N, upperBoundMean,'LineWidth',1.5);
hold on;
plot(d*2:increment:N, (mean(c))*ones(length(d*2:increment:N),1),'k-.','Linewidth',0.75);
hold on;
plot(d*2:increment:N, (mean(c)/sqrt(n))*ones(length(d*2:increment:N),1), 'k-.','Linewidth', 0.75);
legend('Lower bound','Upper bound','\rho','\rho / \surd n','Interpreter','latex');
xlabel('Number of samples (N)')
ylim([0 6])



grid on;
set(gcf,'paperunits','centimeters','papersize',[15 10],'paperposition',[0 0 15 10])
print -dpdf experiment2

figure;
fig2 = gcf;
plot(d*2:increment:N, deltaMean,'LineWidth',1.5)
hold on;
grid on;
set(gcf,'paperunits','centimeters','papersize',[15 10],'paperposition',[0 0 15 10])
xlabel('Number of samples (N)')
ylabel('$\delta(0.95,\omega_N)$','Interpreter','latex');
print -dpdf delta2

return;
load('dataSet1.mat','lowerBound','upperBound')
figure(fig1);

c = jsrRaphael;
plot(d*2:increment:N, upperBound);
hold on;
%legend('Lower bound','Upper bound','\rho','\rho / \surd n','Interpreter','latex');
%xlabel('Number of samples (N)')

%print -dpdf experiment2

figure(fig2);
plot(d*2:increment:N, delta)




load('dataLast.mat','lowerBound','upperBound')
beta = 0.7;
%[lowerBound, upperBound, delta] = computeRhoSequentialSampling(A, beta, N,increment);
%fig1 = gcf;
figure(fig1);
plot(d*2:increment:N, upperBound);



grid on;
set(gcf,'paperunits','centimeters','papersize',[15 10],'paperposition',[0 0 15 10])
%print -dpdf experiment2

figure(fig2);
plot(d*2:increment:N, delta)
beta = 0.5;
rng(15);
[lowerBound, upperBound, delta] = computeRhoSequentialSampling(A, beta, N,increment);
figure(fig1);
plot(d*2:increment:N, upperBound);



grid on;
%print -dpdf experiment2

figure(fig2);
plot(d*2:increment:N, delta)
