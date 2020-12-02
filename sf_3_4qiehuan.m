clear all;
clc;
BB(:,:,1)=[0 1 1 1 0;0 0 1 0 0;0 1 0 0 1;0 0 0 0 1;0 0 0 0 0]; 
BB(:,:,2)=[0 1 0 0 0;1 0 1 0 0;0 1 0 1 0;0 0 1 0 1;0 0 0 0 0]; 
BB(:,:,3)=[0 0 0 0 1;1 0 0 0 0;0 1 0 0 0;0 0 1 0 0;0 0 0 0 0];
BB(:,:,4)=[0 0 1 0 0;1 0 0 0 1;0 1 0 0 0;1 0 0 0 0;0 0 0 0 0];
X=[6;2.4;0;-4.2;1];
f=0.25;
gamma=0.4;
yita=[0;0;0;0;0];
u=[0;0;0;0;0];
t=0:f:20;

  
for n=1:length(t)
    m=round(rand(1,1)*3)+1;
    A=BB(:,:,m);
    X(5,n+1)=sin((n-1)*f)*sin(2*X(5,n))*f+X(5,n); %求导数
    for i=1:length(A)-1
        u(i,n)=0;
        yita=[0;0;0;0;0];
        for K=1:length(A)
            yita(i)=yita(i)+A(i,K);
        end %计算yita
        for j=1:length(A)-1
            if n>=2
                u(i,n)=u(i,n)+(A(i,j)*((X(j,n)-X(j,n-1))/f-gamma*(X(i,n)-X(j,n))))/yita(i);
            else
                u(i,n)=0;
            end
        end
            u(i,n)=u(i,n)+A(i,5)*(sin((n-1)*f)*sin(2*X(5,n))-gamma*(X(i,n)-X(5,n)))/yita(i);
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