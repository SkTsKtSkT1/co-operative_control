%邻接矩阵A用于计算yita aii=0
A=[0 1 1 1 0;0 0 1 0 0;0 1 0 0 1;0 0 0 0 1;0 0 0 0 0];
X=[1.5;1;-1;-1.5;1];
t=0:0.001:20;
gamma=100; 
for n=1:length(t)
    X(5,n)=cos(2*pi*0.001*(n-1));
    for i=1:length(A)
        u(i,n)=0;
        yita(i)=0;
        for K=1:length(A)
            yita(i)=yita(i)+A(i,K);
        end %计算yita
        for j=1:length(A)
            if n>=2
            u(i,n)=u(i,n)+(A(i,j)*((X(j,n)-X(j,n-1))/0.001-gamma*(X(i,n)-X(j,n))))/yita(i);
            else
            u(i,n)=0;
            end
        end
        if n>=2
            u(i,n)=u(i,n)+A(i,5)*((X(5,n)-X(5,n-1))/0.001-gamma*(X(i,n)-X(5,n)))/yita(i);
        else 
            u(i,n)=0;
        end
       X(i,n+1)=X(i,n)+u(i,n)*0.001;
         
    end
        
end

[row,col]=size(X);
figure(1);
for i=1:row
    plot(t,X(i,1:col-1));
    hold on;
end
hold off
