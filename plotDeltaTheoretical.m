clear;
m=2;
Nmax =  85000;
beta = 0.95;
close all;
figure;
grid on;
counter = 1;
for n=2:1:6
    [delta, N] = computeDeltaTheoretical(n,m,beta,Nmax);
    %plot(N, (delta(counter,:)));
    semilogx(N, delta,'LineWidth',1.5);
    %loglog(N, delta(counter,:));
    counter = counter +1;
    hold on;
    grid on;
end

xlabel('Number of Samples(N)');
ylabel('$\delta(0.95, N)$', 'Interpreter','latex')
legend('n=2','n=3','n=4','n=5','n=6');
set(gcf,'paperunits','centimeters','papersize',[15 10],'paperposition',[0 0 15 10])
%%
%print -dpdf deltaWrtN
figure;
n = 3;
for m = 2:1:6
    [delta, N] = computeDeltaTheoretical(n,m,beta,Nmax);
    %plot(N, (delta(counter,:)));
    semilogx(N, delta,'LineWidth',1.5);
    %loglog(N, delta(counter,:));
    counter = counter +1;
    hold on;
    grid on;
end

xlabel('Number of Samples(N)');
ylabel('$\delta(0.95, N)$', 'Interpreter','latex')
legend('m=2','m=3','m=4','m=5','m=6');
set(gcf,'paperunits','centimeters','papersize',[15 10],'paperposition',[0 0 15 10])