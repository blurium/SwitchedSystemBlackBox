%this function performs random experiments with 2D, 3D, 4D, 5D systems
close all;
clear;

rng(10);
Nsystems = 10;
n = 3;
m = 5;
maxJSR = 1.5;
beta = 0.9;
tic;
for i = 1:Nsystems
    [A, jsrRaphael{i}] = createRandomSystem(n,m,maxJSR);
    upperBoundRaphael(i) = jsrRaphael{i}(2);
    lowerBoundRaphael(i) = jsrRaphael{i}(1);
    N = 3000;
    increment = 200;
    d = (n*(n+1))/2 +1;
    [lowerBound{i}, upperBound{i}, delta{i}] = computeRhoSequentialSampling(A, beta, N,increment);  
end
toc;



Nsystems = 10;
n = 4;
m = 6;
tic;
for i = 1:Nsystems
    [A2, jsrRaphael2{i}] = createRandomSystem(n,m,maxJSR);
    upperBoundRaphael2(i) = jsrRaphael2{i}(2);
    lowerBoundRaphael2(i) = jsrRaphael2{i}(1);
    N = 3000;
    increment = 200;
    d = (n*(n+1))/2 +1;
    [lowerBound2{i}, upperBound2{i}, delta2{i}] = computeRhoSequentialSampling(A2, beta, N,increment);  
end
toc;

%save randomExperiments

%now compute averages
for i = 1:length(lowerBound2{1})   
    lowerBoundMean(i) = mean(cellfun(@(v) v(i), lowerBound)./lowerBoundRaphael);
    upperBoundMean(i) = mean(cellfun(@(v) v(i), upperBound)./upperBoundRaphael);
    lowerBoundMean2(i) = mean(cellfun(@(v) v(i), lowerBound2)./lowerBoundRaphael2);
    upperBoundMean2(i) = mean(cellfun(@(v) v(i), upperBound2)./upperBoundRaphael2);
end

%%


Nsystems = 10;
n = 4;
m = 5;
tic;
for i = 1:Nsystems
    [A3, jsrRaphael3{i}] = createRandomSystem(n,m,maxJSR);
    upperBoundRaphael3(i) = jsrRaphael3{i}(2);
    lowerBoundRaphael3(i) = jsrRaphael3{i}(1);
    N = 3000;
    increment = 200;
    d = (n*(n+1))/2 +1;
    [lowerBound3{i}, upperBound3{i}, delta3{i}] = computeRhoSequentialSampling(A3, beta, N,increment);  
end
toc;
%%
save randomExperiments2

%now compute averages
for i = 1:length(lowerBound2{1})   
    lowerBoundMean(i) = mean(cellfun(@(v) v(i), lowerBound)./lowerBoundRaphael);
    upperBoundMean(i) = mean(cellfun(@(v) v(i), upperBound)./upperBoundRaphael);
    lowerBoundMean2(i) = mean(cellfun(@(v) v(i), lowerBound2)./lowerBoundRaphael2);
    upperBoundMean2(i) = mean(cellfun(@(v) v(i), upperBound2)./upperBoundRaphael2);
    lowerBoundMean3(i) = mean(cellfun(@(v) v(i), lowerBound3)./lowerBoundRaphael3);
    upperBoundMean3(i) = mean(cellfun(@(v) v(i), upperBound3)./upperBoundRaphael3);
end




%%
figure;
fig1 = gcf;
c = jsrRaphael;
plot(d*2:increment:N, lowerBoundMean,'LineWidth',1.5);
hold on;
plot(d*2:increment:N, lowerBoundMean3,'LineWidth',1.5);
hold on;
plot(d*2:increment:N, lowerBoundMean2,'LineWidth',1.5);
hold on;
grid on;
legend('n=3, m=5', 'n=4, m=5','n=4, m=6');
xlabel('Number of samples (N)');
ylabel('Lower bound (black-box)/Lower bound (white-box)');
set(gcf,'paperunits','centimeters','papersize',[15 10],'paperposition',[0 0 15 10])
print -dpdf randomExperiment1


figure;
fig2 = gcf;
plot(d*2:increment:N, upperBoundMean,'LineWidth',1.5);
hold on;
plot(d*2:increment:N, upperBoundMean3,'LineWidth',1.5);
hold on;
plot(d*2:increment:N, upperBoundMean2,'LineWidth',1.5);
hold on;
grid on;
legend('n=3, m=5', 'n=4, m=5','n=4, m=6');
xlabel('Number of samples (N)');
ylabel('Upper bound (black-box)/Upper bound (white-box)');
set(gcf,'paperunits','centimeters','papersize',[15 10],'paperposition',[0 0 15 10])
print -dpdf randomExperiment2




%now plotting begins






