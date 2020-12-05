A=[0 0 1 0;1 0 0 0;0 1 0 0;1 0 0 0];
B=[0 0 1 0;1 0 0 0;0 1 0 0;1 0 0 0];
X=[0;0.2;0.4;0.6];  %位置
Y=[0;0.1;0.2;0.3];  % 其他参数
Kr=[1 0 0 0;0 1 0 0;0 0 1 0;0 0 0 1];
Kv=[1 0 0 1;0 1 0 0;1 0 1 0;0 0 0 1];
R=diag(Kr);
V=diag(Kv);
u=[0;0;0;0];
f=0.001;
t=0:f:20;

for n=1:length(t)
    for i=1:length(A)
        u(i,n)=0;
    for j=1:length(A)
       u(i,n)=u(i,n)-(A(i,j)*tanh(R(i)*(X(i,n)-X(j,n)))+B(i,j)*(V(i)*(Y(i,n)-Y(j,n))));
      % u(i,n)=u(i,n)-A(i,j)*tanh(R(i)*(X(i,n)-X(j,n)));
    end
    %u(i,n)=u(i,n)-tanh(V(i)*Y(i,n));
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
figure(3)
plot(t,u(1,1:col1-1),t,u(2,1:col1-1),t,u(3,1:col1-1),t,u(4,1:col1-1));
hold on;

        