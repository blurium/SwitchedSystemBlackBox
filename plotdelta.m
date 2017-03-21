delta=[0  0    0.6806    0.6035    0.9876    0.9927    0.9952    0.9966    0.9974    0.9980    0.9984];

lb=[0.4917    0.5138    0.5138    0.5359    0.5469    0.5469    0.5469    0.5469    0.5469    0.5469    0.5469];

ub=[1.0675    1.2558    0.7831    0.7791    0.7772    0.7761    0.7754    0.7750    0.7747];

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