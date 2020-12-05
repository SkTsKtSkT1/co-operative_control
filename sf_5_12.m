A=[0 0 0 0 0 0 1;1 0 0 0 0 0 0;0 1 0 0 0 0 0;1 0 0 0 0 0 0;0 0 0 1 0 0 0;0 0 0 0 1 0 0;0 0 0 0 0 0 0];
X=[0;0.2;0.4;0.6;0.8;1.0;1];  %位置
Y=[0;0.1;0.2;0.3;0.4;0.5;1];  % 其他参数

Kr=[1 0 0 0 0 0 0;0 1 0 0 0 0 0;0 0 1 0 0 0 0;0 0 0 1 0 0 0;0 0 0 0 1 0 0;0 0 0 0 0 1 0;0 0 0 0 0 0 1];
Kv=[1 0 0 0 0 0 0;0 1 0 0 0 0 0;0 0 1 0 0 0 0;0 0 0 1 0 0 0;0 0 0 0 1 0 0;0 0 0 0 0 1 0;0 0 0 0 0 0 1];
R=diag(Kr);
V=diag(Kv);
f=0.001;
t=0:f:20;

for n=1:length(t)
    Y(5,n+1)=sin((n-1)*f)*sin(2*Y(7,n))*f+Y(7,n);
    for i=1:length(A)-1
        u(i,n)=0;
        if n>=2
            if i==1
                u(i,n)=(Y(7,n)-Y(7,n-1))/f-R(i)*(X(i,n)-X(7,n))-V(i)*(Y(i,n)-Y(7,n));
            else
                u(i,n)=(Y(1,n)-Y(1,n-1))/f-R(i)*(X(i,n)-X(1,n))-V(i)*(Y(i,n)-Y(1,n));
            end
        else
            u(i,n)=0;
        end
        Y(i,n+1)=Y(i,n)+f*u(i,n);
        X(i,n+1)=X(i,n)+f*Y(i,n);
        Y(7,n+1)=Y(7,n)+sin((n-1)*f)*sin(2*X(7,n))*sin(2*Y(7,n))*f;
        X(7,n+1)=X(7,n)+f*Y(7,n);
    end
end



[row,col]=size(X);
[row1,col1]=size(Y);
figure(1);

for i=1:row
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