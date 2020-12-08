A=[0 1 0 1;1 0 1 0;0 1 0 0;1 0 0 0];
X=[0;0.2;0.4;0.6];  %位置
Y=[0;0.1;0.2;0.3];  % 其他参数
%%观测器
% X_hat=rand([4 1]);
% Y_hat=zeros(4,1);
X_hat=X;
Y_hat=Y;
u=[0;0;0;0];
f=0.001;
t=1:f:20;

Q= [1.0000    0.5000    0.3333    0.2500;
    0.5000    1.0000    0.6667    0.5000;
    0.3333    0.6667    1.0000    0.7500;
    0.2500    0.5000    0.7500    1.0000];
F=[-1 0 0 0;0 -2 0 0;0 0 -3 0;0 0 0 -4]; %赫尔维茨矩阵->特征值均具有负实部->对角线元素为根
P = lyap(F', Q);
PF=P*F;
for n=1:length(t)
    sum=[0;0;0;0];
    x_hat_dot=[0;0;0;0];
    for i=1:length(A)
        u(i,n)=0;
        for j=1:length(A)
            sum(i)=sum(i)+A(i,j)*(X(i,n)-X(j,n));  %4*4*4*1=4*1
        end
        x_hat_dot(i)=F(i,:)*X_hat(:,n)+X(i,n);  %1*1+1*1+1*1
        X_hat(i,n+1)=X_hat(i,n)+x_hat_dot(i)*f;
        Y_hat(i,n)=PF(i,:)*X_hat(:,n)+det(P)*X(i,n);
        u(i,n)=-sum(i)-Y_hat(i);
        Y(i,n+1)=Y(i,n)+f*u(i,n);
        X(i,n+1)=X(i,n)+f*Y(i,n);
    end
  
end

[row,col]=size(X);
[row1,col1]=size(Y);
figure(1);
for i=1:row
    plot(t,X(i,1:col-1));
    hold on;
end
hold off
figure(2);
for i=1:row1
    plot(t,Y(i,1:col1-1));
    hold on;
end
hold off



