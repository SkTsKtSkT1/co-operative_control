A=[0 0 1 0 1;1 0 1 0 0;0 1 0 0 0;1 0 0 0 0;0 0 0 0 0];
X=[-0.5;0;0.5;1.5;1];
t=0:0.01:20;
for n=1:length(t)
     for i=1:size(A)
        u(i,n)=0;
        for j=1:size(A)
            u(i,n)=u(i,n)-A(i,j)*(X(i,n)-X(j,n));
        end
     end
     for i=1:size(A)
         X(i,n+1)=X(i,n)+u(i,n)*0.01;
     end
end
plot(t,X(1,1:length(t)),t,X(2,1:length(t)),t,X(3,1:length(t)),t,X(4,1:length(t)),t,X(5,1:length(t)));grid on