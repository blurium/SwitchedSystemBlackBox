x=[0:pi/300:2*pi];
%set(0,'defaultAxesFontSize', 40);
set(gca,'FontSize',15)
hold on
grid on
axis([-1.1 1.1 -1.1 1.1])

x1=[1.20:pi/300:1.94];
x2=[4.34:pi/300:5.08];
plot(cos(x),sin(x),'k')
plot(cos(x1),sin(x1),'r','LineWidth',2)
plot(cos(x2),sin(x2),'r','LineWidth',2)

y=[0:pi/20:2*pi];
y1=cos(y);
y2=sin(y);
q=quiver(y1,y2,-2*y1,0.3*y2,'MaxHeadSize',0.1);

xlabel('x1','FontWeight','bold')
ylabel('x2','FontWeight','bold')

pbaspect([1 1 1])

U=q.UData;
V=q.VData;
X=q.XData;
Y=q.YData;

figure
set(gca,'FontSize',15)
hold on
grid on
axis([-1.1 1.1 -1.1 1.1])
plot(cos(x),sin(x),'k')
plot(cos(x1),sin(x1),'r','LineWidth',2)
plot(cos(x2),sin(x2),'r','LineWidth',2)

X=[X 0 0];
Y=[Y 1 -1];
U=[U 0 0];
V=[V 0.3 -0.3];
for i = 1:length(Y)
        ah = annotation('arrow','Color','b',...
            'headStyle','cback1','HeadLength',5,'HeadWidth',5);
        set(ah,'parent',gca);
        set(ah,'position',[X(i) Y(i) 0.2*U(i) 0.2*V(i)]);
end

pbaspect([1 1 1])
xlabel('x1','FontWeight','bold')
ylabel('x2','FontWeight','bold')


figure
set(gca,'FontSize',15)
hold on
grid on
axis([-1.1 1.1 -1.1 1.1])
x1=[1.50:pi/300:1.64];
x2=[4.64:pi/300:4.78];

plot(cos(x),sin(x),'k')
plot(cos(x1),sin(x1),'r','LineWidth',2)
plot(cos(x2),sin(x2),'r','LineWidth',2)

q1=quiver(y1,y2,-2*y1,0.01*y2);

xlabel('x1','FontWeight','bold')
ylabel('x2','FontWeight','bold')
pbaspect([1 1 1])

U1=q1.UData;
V1=q1.VData;
X1=q1.XData;
Y1=q1.YData;

figure
set(gca,'FontSize',15)
hold on
grid on
axis([-1.1 1.1 -1.1 1.1])
plot(cos(x),sin(x),'k')
plot(cos(x1),sin(x1),'r','LineWidth',2)
plot(cos(x2),sin(x2),'r','LineWidth',2)
X1=[X1 0 0];
Y1=[Y1 1 -1];
U1=[U1 0 0];
V1=[V1 0.01 -0.01];

for i = 1:length(Y1)
        ah = annotation('arrow','Color','b',...
            'headStyle','cback1','HeadLength',5,'HeadWidth',5);
        set(ah,'parent',gca);
        set(ah,'position',[X1(i) Y1(i) 0.2*U1(i) 0.2*V1(i)]);
end


xlabel('x1','FontWeight','bold')
ylabel('x2','FontWeight','bold')
pbaspect([1 1 1])
%hTitle = title('Title of the plot');
%set(hTitle,'FontSize',30)