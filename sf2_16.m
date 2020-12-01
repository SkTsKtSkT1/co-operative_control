BB(:,:,1)=[0 0 1 0 0;0 0 0 0 0;0 0 0 0 0;0 0 0 0 0;0 0 0 0 0];
BB(:,:,2)=[0 0 0 0 0;1 0 0 0 0;0 0 0 0 0;0 0 0 0 0;0 0 0 0 0];
BB(:,:,3)=[0 0 0 0 0;0 0 0 0 0;1 0 0 0 0;0 0 0 0 0;0 0 0 0 0];
BB(:,:,4)=[0 0 0 0 0;0 0 0 0 0;0 0 0 0 0;0 0 1 0 0;0 0 0 0 0];
BB(:,:,5)=[0 0 0 0 0;0 0 0 0 0;0 0 0 0 0;0 0 0 0 0;0 0 0 1 0];
X=[0.2;0.4;0.6;0.8;1];
t=0:0.01:50;
for n=1:length(t)
    m=round(rand(1,1)*4)+1;
    A=BB(:,:,m);
    for i=1:size(A)
        u(i,n)=0;
        for j=1:size(A)
            u(i,n)=u(i,n)-A(i,j)*(X(i,n)-X(j,n));
        end
    end
    for i=1:size(A)
        X(i,n+1)=X(i,n)+0.01*u(i,n);
    end
end
plot(t,X(1,1:length(t)),t,X(2,1:length(t)),t,X(3,1:length(t)),t,X(4,1:length(t)),t,X(5,1:length(t)));grid on
