delta=[0.8036 0.8058 0.8900 0.9249 0.9493 0.9575 0.9695 0.9728 0.9754 0.9774 0.9790];

lb=[0.4751 0.5204 0.5370 0.5392 0.5458 0.5458 0.5469 0.5469 0.5469 0.5469 0.5469];

ub=[0.8944 0.9132 0.8547 0.8247 0.8133 0.8063 0.7977 0.7950 0.7930 0.7913 0.7900];

points=[15:50:515];
set(gca,'FontSize',15)

a=ones(1,11);
hold on
grid on
plot(points,delta,'b','LineWidth',2)
plot(points,a,'--r','LineWidth',1)
ylim([0 1.2]);
xlim([0 520]);



xlabel('Number of samples (N)','FontWeight','bold')

ylabel('$\delta(0.92,\omega_N)$','Interpreter','latex','FontWeight','bold')
