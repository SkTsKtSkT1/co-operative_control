clear all;
clc;
A=[0 1 1 1 0;0 0 1 0 0;0 1 0 0 1;0 0 0 0 1;0 0 0 0 0];
X=[2;0.5;0;-0.5;0.5];
f=0.001;
yita=[0;0;0;0;0];
u=[0;0;0;0;0];
delta=[0;0;0;0;0];
gamma=1;
t=0:f:20;

  
for n=1:length(t)
    X(5,n+1)=sin((n-1)*f)*sin(2*X(5,n))*f+X(5,n); %求导数
    for i=1:length(A)-1
        u(i,n)=0;
        yita=[0;0;0;0;0];
        for K=1:length(A)
            yita(i)=yita(i)+A(i,K);
            delta(K,n)=1-K;
        end %计算yita
        
        for j=1:length(A)-1
            if n>=2
                test=(delta(j,n)-delta(j,n-1))/f;
                test2=X(i,n)-X(j,n);
                tset3=delta(i,n)-delta(j,n);
                u(i,n)=u(i,n)+(A(i,j)*((X(j,n)-X(j,n-1))/f-(delta(j,n)-delta(j,n-1))/f-gamma*((X(i,n)-X(j,n))-(delta(i,n)-delta(j,n)))))/yita(i);
            else
                u(i,n)=0;
            end
        end
            if n>=2
                test1=(delta(i,n)-delta(i,n-1))/f;
                u(i,n)=u(i,n)+(delta(i,n)-delta(i,n-1))/f+A(i,5)*(sin((n-1)*f)*sin(2*X(5,n))-gamma*(X(i,n)-delta(i,n)-X(5,n)))/yita(i);
            else
                u(i,n)=0;
            end
            X(i,n+1)=X(i,n)+u(i,n)*f;        
    end   
end
  

[row,col]=size(X);
figure(1);
for i=1:row
    plot(t,X(i,1:col-1));
    hold on;
end
hold off