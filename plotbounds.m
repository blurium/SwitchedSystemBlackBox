close all;
delta=[0         0    0.6863    0.8490    0.9120    0.9425    0.9596    0.9700    0.9769    0.9817    0.9851    0.9876    0.9896    0.9911    0.9923    0.9933    0.9941    0.9948    0.9953 0.9958    0.9962    0.9966    0.9969    0.9971];

lb=[0.4972    0.4972    0.5469    0.5469    0.5469    0.5469    0.5469    0.5469    0.5469    0.5469    0.5469    0.5469    0.5469    0.5469    0.5469    0.5469    0.5469    0.5469    0.5469  0.5469    0.5469    0.5469    0.5469    0.5469];

ub=[1.1270    0.9110    0.8481    0.8206    0.8060    0.7973    0.7917    0.7879    0.7851    0.7831    0.7816    0.7804    0.7794    0.7787    0.7780    0.7775    0.7771 0.7767    0.7764    0.7761    0.7759    0.7757];


points=[15:15:360];
points2=[45:15:360];
set(gca,'FontSize',20)

a=0.7727*ones(1,24);
b=1/sqrt(2)*a;
hold on
grid on
plot(points,lb,'b','LineWidth',2)
plot(points2,ub,'r','LineWidth',2)
plot(points,a,'--k','LineWidth',1)
plot(points,b,'-.k','LineWidth',1)
ylim([0.3 1.2]);
xlim([0 380]);



xlabel('Number of samples (N)','FontWeight','bold')



h_legend = legend('Lower bound','Upper bound','\rho','\rho / \surd n','Interpreter','latex','FontSize')
set(h_legend,'FontSize',20);
%legend('\sqrt{n}','Interpreter','latex')
%ylabel('Bounds','FontWeight','bold')