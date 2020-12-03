BB(:,:,1)=[0 1 0;0 0 0;0 0 0];
BB(:,:,2)=[0 1 0;0 0 1;2 0 0];
BB(:,:,3)=[0 0 0;0 0 1;0 0 0];


X=[-0.1; 0 ;0.5];  %位置
Y=[0.1;-0.2;0.3];  % 其他参数
u=[0;0;0;0];
f=0.001;
r=[1;1;1];  %r(t)>0
t=1:f:30;

for n=1:length(t)
    %     Q=round(10*rand(1,1));
    %     if(Q>9)
    %         A=BB(:,:,2);
    %     else
    %         A=BB(:,:,1);
    %     end
    if(mod((n-1)*f,1)<=1*0.9)
        A=BB(:,:,1);
    else
        A=BB(:,:,2);
    end

    for i=1:length(A)
        u(i,n)=0;
        for j=1:length(A)
            u(i,n)=u(i,n)-A(i,j)*((X(i,n)-X(j,n))+r(i)*(Y(i,n)-Y(j,n)));
        end
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
        