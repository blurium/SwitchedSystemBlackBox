close all;
delta=[    0.8100    0.8100    0.8307    0.8929    0.9265    0.9349    0.9419    0.9474    0.9517    0.9546    0.9576    0.9602    0.9624    0.9643    0.9660    0.9671];

lb=[   0.4375    0.4453    0.4609    0.4688    0.4766    0.4766    0.4766    0.4766    0.4766    0.4766    0.4766    0.4766    0.4766    0.4766    0.4766    0.4766];

ub=[  1.0802    1.0995    1.1098    1.0499    1.0288    1.0195    1.0119    1.0061    1.0015    0.9984    0.9953    0.9926    0.9904    0.9884    0.9867    0.9856];


points=[65:50:815];
set(gca,'FontSize',20)

a=0.958*ones(1,16);
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