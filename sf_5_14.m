A=[0 0 0 0 0 0 1;1 0 1 0 0 0 0;0 1 0 0 0 0 0 ;1 0 0 0 1 0 0;0 0 0 0 0 0 1;0 0 0 0 1 0 0;0 0 0 0 0 0 0];
X=[1;2;3;4;5;6;pi/2];  %位置
Y=[0.1;0.2;0.3;0.4;0.5;0.6;0];  % 其他参数

Kr=1*[1 0 0 0 0 0 0;0 1 0 0 0 0 0;0 0 1 0 0 0 0;0 0 0 1 0 0 0;0 0 0 0 1 0 0;0 0 0 0 0 1 0;0 0 0 0 0 0 1];
Kv=1*[1 0 0 0 0 0 0;0 1 0 0 0 0 0;0 0 1 0 0 0 0;0 0 0 1 0 0 0;0 0 0 0 1 0 0;0 0 0 0 0 1 0;0 0 0 0 0 0 1];
R=diag(Kr);
V=diag(Kv);
k=[];
f=0.001;
t=0:f:20;

for n=1:length(t)
    for i=1:length(A)-1
        u(i,n)=0;
        k(i)=0;
        for k_i=1:length(A)
            k(i)=k(i)+A(i,k_i);
        end
        for j=1:length(A)-1
            if n>=2
                u(i,n)=u(i,n)+A(i,j)*((Y(j,n)-Y(j,n-1))/f-R(i)*(X(i,n)-X(j,n))-V(i)*(Y(i,n)-Y(j,n)))/k(i);
            else
                u(i,n)=0;
            end
        end
            if n>=2
                u(i,n)=u(i,n)+A(i,7)*((Y(7,n)-Y(7,n-1))/f-R(i)*(X(i,n)-X(7,n))-V(i)*(Y(i,n)-Y(7,n)))/k(i);
            end
            Y(i,n+1)=Y(i,n)+f*u(i,n);
            X(i,n+1)=X(i,n)+f*Y(i,n);
            Y(7,n+1)=Y(7,n)+(-sin(X(7,n))/(1+exp(-(n-1)*f)))*f;
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