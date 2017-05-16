clear;
rng(75)
close all;
jsrRaphael = [inf inf];
while(jsrRaphael(2)>1 || jsrRaphael(2)<0.7) %0.7
Aini = rand(2,2);

B = [0.5 1]';

K = place(Aini,B,[0.8 -0.7]);


Abar = Aini - B*K;

A{1} = Abar^2*Abar^2*Aini^4;
A{2} = Aini*Aini*Abar^2*Aini^4;
A{3} = Aini*Abar*Abar^2*Aini^4;
A{4} = Abar*Aini*Abar^2*Aini^4;


jsrRaphael = jsr_prod_bruteForce(A)
end
n=2;
m=4;
beta = 0.9;
N = 900;
increment = 50;


[lowerBound, upperBound, delta] = computeRhoSequentialSampling(A, beta, N,increment);
fig1 = gcf;

%%
close all
jsrRaphael = jsr_prod_pruningAlgorithm(A)
d = (n*(n+1))/2 +1;
c = jsrRaphael;
plot(d*2:increment:N, lowerBound,'LineWidth',1.5);
hold on;
plot(d*2:increment:N, upperBound,'LineWidth',1.5);
hold on;
plot(d*2:increment:N, (mean(c))*ones(length(d*2:increment:N),1),'k-.','Linewidth',0.75);
hold on;
plot(d*2:increment:N, (mean(c)/sqrt(n))*ones(length(d*2:increment:N),1), 'k-.','Linewidth', 0.75);
legend('Lower bound','Upper bound','\rho','\rho / \surd n','Interpreter','latex');
xlabel('Number of samples (N)')
ylim([0 2.5])
xlim([0 750])
grid on;
set(gcf,'paperunits','centimeters','papersize',[15 10],'paperposition',[0 0 15 10])
print -dpdf networkControl