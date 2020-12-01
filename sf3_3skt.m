%A=[-2 0 1 0 1;1 -2 1 0 0;0 1 -1 0 0;1 0 0 -1 0;0 0 0 0 0];  %(n+1)*(n+1)  n=4
A=[-2 0 1 0 1;1 -3 1 0 1;0 1 -2 0 1;1 0 0 -2 1;0 0 0 0 0];  %(n+1)*(n+1)  n=4
X=[-0.1;-0.3;0.1;0.3;0]; %初始状态
t=0:0.001:20;
k=100;  %% \alpha

for n=1:length(t)
    for i=1:size(A)
        u(i,n)=0;
        X(5,n)=0.1*cos(2*pi*0.001*(n-1));  %基准状态为余弦信号
        for j=1:size(A)
            u(i,n)=u(i,n)-A(i,j)*(X(i,n)-X(j,n))-A(i,5)*k*(X(i,n)-X(5,n))+A(i,5)*0.1*sin(2*pi*(n-1)*0.001);
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