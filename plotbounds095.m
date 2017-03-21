close all;
delta=[0.8100    0.8100    0.8100    0.8192    0.8558    0.8675    0.9280    0.9358    0.9347    0.9317    0.9564];

lb=[0.3984    0.4375    0.4570    0.4609    0.4609    0.4609    0.4688    0.4688    0.4688    0.4688    0.4727];

ub=[0.9837    1.0802    1.1284    1.1253    1.0772    1.0627    1.0102    1.0018    1.0030    1.0062    0.9885];


points=[15:50:515];
set(gca,'FontSize',20)

a=0.958*ones(1,11);
b=1/sqrt(4)*a;
hold on
grid on
plot(points,lb,'b','LineWidth',2)
plot(points,ub,'r','LineWidth',2)
plot(points,a,'--k','LineWidth',1)
plot(points,b,'-.k','LineWidth',1)
ylim([0 1.2]);
xlim([0 520]);



xlabel('Number of samples (N)','FontWeight','bold')



h_legend = legend('Lower bound','Upper bound','\rho','\rho / \surd n','Interpreter','latex','FontSize')
set(h_legend,'FontSize',20);
%legend('\sqrt{n}','Interpreter','latex')
%ylabel('Bounds','FontWeight','bold')