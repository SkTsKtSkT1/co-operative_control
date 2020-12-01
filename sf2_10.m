A=[0 1.5 0 0 0 0;2 0 0 0 0 0;0.9 0 0 0 1.9 0;0 1.2 0 0 0 1.3;0 0 1.4 1.8 0 0;0 0 0 0 0.7 0];
%A=5*a;
tao=[0.1;0.3;0.5;0.7;0.9;1.1];
t=0:0.01:20;
for n=1:length(t)
    for i=1:size(A)
        u(i,n)=0;
        u1(i,n)=0;
        for j=1:size(A)
            u1(i,n)=u1(i,n)-A(i,j)*(tao(i,n)-tao(j,n));
        end
        u(i,n)=u1(i,n)+0.2*abs(sin(2*pi*0.01*(n-1)));
    end
    for i=1:size(L)
        tao(i,n+1)=tao(i,n)+0.01*u(i,n);
    end
end
plot(t,tao(1,1:length(t)),t,tao(2,1:length(t)),t,tao(3,1:length(t)),t,tao(4,1:length(t)),t,tao(5,1:length(t)),t,tao(6,1:length(t)))