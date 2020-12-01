A=5*[-1.5,1.5,0,0,0,0;2,-2,0,0,0,0;0.9,0,-2.8,0,1.9,0;0,1.2,0,-2.5,0,1.3;0,0,1.4,1.8,-3.2,0;0,0,0,0,0.7,-0.7]; %-Ln
X=[0.1;0.3;0.5;0.7;0.9;1.1];
tao=[0.2 ;0.4 ;0.6 ;0.8;1.0;1.2];
t=0:0.001:20;
for n=1:length(t)
     for i=1:size(A)
        u(i,n)=0;
        utao(i,n)=0;
        for j=1:size(A)
            u(i,n)=u(i,n)-A(i,j)*(X(i,n)-X(j,n));
            utao(i,n)= utao(i,n)-A(i,j)*(tao(i,n)-tao(j,n));
        end
        X(i,n+1)=X(i,n)+u(i,n)*0.001;
        tao(i,n+1)=tao(i,n)+utao(i,n)*0.001;
     end
end
[row,col]=size(X);
figure(1);
for i=1:row
    plot(t,X(i,1:col-1));
    hold on;
end
hold off
figure(2);
for i=1:row
    plot(t,tao(i,1:col-1));
    hold on;
end
hold off



