A=[0 0 1 0 1;1 0 0 0 1;0 1 0 0 1;1 0 0 0 1;0 0 0 0 0];
X=[0;0.2;0.4;0.6;1];  %位置
Y=[0;0.1;0.2;0.3;1];  % 其他参数
u=[0;0;0;0];
gamma=1;
alpha=1;
f=0.001;
t=0:f:20;

for n=1:length(t)
    Y(5,n+1)=sin((n-1)*f)*sin(2*Y(5,n))*f+Y(5,n);
    for i=1:length(A)-1
        u(i,n)=0;
        for j=1:length(A)-1
            u(i,n)=u(i,n)-A(i,j)*((X(i,n)-X(j,n))+gamma*(Y(i,n)-Y(j,n)));
        end
        if n>=2
            u(i,n)=u(i,n)+(Y(5,n)-Y(5,n-1))/f-alpha*(Y(i,n)-Y(5,n));
        else
            u(i,n)=0;
        end
        Y(i,n+1)=Y(i,n)+f*u(i,n);
        X(i,n+1)=X(i,n)+f*Y(i,n);
    end
end

[row,col]=size(X);
[row1,col1]=size(Y);
figure(1);

for i=1:row-1
    plot(t,X(i,1:col-1));
    hold on;
end
title('\xi');
hold off
figure(2);

for i=1:row1
    plot(t,Y(i,1:col1-1));
    hold on;
end
title('\zeta');
hold off
