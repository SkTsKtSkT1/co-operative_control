%A=[0 1 1 1 0;0 0 1 0 0;0 1 0 0 1;0 0 0 0 1;0 0 0 0 0]; %3.8
A=[0 1 0 0 0;1 0 1 0 0 ;0 1 0 1 0;0 0 1 0 1;0 0 0 0 0];
%X=[6;2.4;0;-4.2;1]; %3.8
X=[1.5;0.5;0;-0.5;1];
f=0.001;
lambda=10*[1 0 0 0 0;0 1 0 0 0;0 0 1 0 0;0 0 0 1 0;0 0 0 0 1];
D=diag(lambda);
yita=[0;0;0;0;0];
u=[0;0;0;0;0];
t=0:f:20;
for n=1:length(t)
    X(5,n+1)=sin((n-1)*f)*sin(2*X(5,n))*f+X(5,n); %Çóµ¼Êý
    if((n-1)*f>5 && (n-1)*f<7)
       Zf=50*abs(X(4,n)-X(3,n));
        X(4,5)=1/(1+Zf);
    end
    for i=1:length(A)-1
        u(i,n)=0;
        yita=[0;0;0;0;0];
        qwq=0;
        for j=1:length(A)
            yita(i)=yita(i)+A(i,j);
        end
        for j=1:length(A)
            qwq=qwq+A(i,j)*(X(i,n)-X(j,n));
            if n>=2
                u(i,n)=u(i,n)+(A(i,j)*((X(j,n)-X(j,n-1))/f)/yita(i));
            else
                u(i,n)=0;
            end
        end
        
        u(i,n)=u(i,n)+A(i,5)*(sin((n-1)*f)*sin(2*X(5,n)))/yita(i)-D(i)*tanh(qwq+A(i,5)*(X(i,n)-X(5,n)))/yita(i);
        if((n-1)*f>5 && (n-1)*f<7)
            u(3,n)=0; 
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
% figure(2);
% for i=1:row-1
%     plot(t,u(i,1:col-1));
%     hold on;
% end
% hold off


        
            