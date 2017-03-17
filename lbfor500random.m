x=[1:500];
result=zeros(1,500);
ro=zeros(1,500);
for u=1:500
testlb;
result(u)=lowerBound;
ro(u)=c(1);
u
end
plot(x,result,x,ro)