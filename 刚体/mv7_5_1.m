% 算法7.5
% 例7.3-P106
% 由有向图7.2-P106可知，vL={1},N1={4},N2={1},N3={2},N4={5},N5={2,6},N6={3}
% 故此，J1=N1并{r}={4,r};J2=N2={1};J3=N3={2};J4=N4={5};J5=N5={2,6};J6=N6={3};(注意这里写的Ji并不是刚体扭矩，只是证明中设置的一个中间参数)

clear 
clc
close all

n=6;
A=[0 0 0 1 0 0 1;      % 有向图7.2-P106且GnA与GnB是相同的
   1 0 0 0 0 0 0;
   0 1 0 0 0 0 0;
   0 0 0 0 1 0 0;
   0 1 0 0 0 1 0;
   0 0 1 0 0 0 0;
   0 0 0 0 0 0 0];
D=diag(sum(A, 2));           % 构建入度矩阵D
L=D-A;                       % L(n+1)*(n+1)         
M=L(1:n,1:n);                % L(n+1)*(n+1)去除第n+1行第n+1列之后才是M
b=L(1:n,n+1);                % L(n+1)*(n+1)第n+1列去掉最下方一个元素才是b

J1=[1 0.1 0.1;0.1 0.1 0.1;0.1 0.1 0.9];    % 6个刚体的惯性扭矩
J2=[1.5 0.2 0.3;0.2 0.9 0.4;0.3 0.4 2.0];
J3=[0.8 0.1 0.2;0.1 0.7 0.3;0.2 0.3 1.1];
J4=[1.2 0.3 0.7;0.3 0.9 0.2;0.7 0.2 1.4];
J5=[0.9 0.15 0.3;0.15 1.2 0.4;0.3 0.4 1.2];
J6=[1.1 0.35 0.45;0.35 1.0 0.5;0.45 0.5 1.3];

kq1=1;
kq2=1;
kq3=1;
kq4=1;
kq5=1;
kq6=1;

Kw1=2*eye(3);
Kw2=2*eye(3);
Kw3=2*eye(3);
Kw4=2*eye(3);
Kw5=2*eye(3);
Kw6=2*eye(3);

taor=[0;0;0];           
Jr=diag([1,2,1]);
qr_initial=[1;0;0;0];
wr_initial=[0.1;0.3;0.5];

q_initial=normalize(quaternion(rand(n,4)));   % q_initial是n*1维的四元数矢量
q_initial_quatParts=compact(q_initial);       % q_initial_quatParts(n*4)是n*1维的四元数矢量q_initial的系数矩阵
s_initial=q_initial_quatParts(:,1);           % s_initial也可以像w_guodu一样直接构建单位化的矩阵，而不是先通过构造单位四元数矢量然后再取出系数矩阵
x_initial=q_initial_quatParts(:,2);
y_initial=q_initial_quatParts(:,3);
z_initial=q_initial_quatParts(:,4);

w_guodu=rand(n,3);
for i=1:n
    w_guodu(i,:)=w_guodu(i,:)/norm(w_guodu(i,:));
end
% wbo_initial_quatParts=[zeros(n,1),w_guodu];   % wbo_initial_quatParts(n*4)是n*1维的四元数矢量wbo_initial的系数矩阵
% wbo_initial=quaternion(wbo_initial_quatParts);    % wbo_initial是n*1维的四元数矢量
wx_initial=w_guodu(:,1);
wy_initial=w_guodu(:,2);
wz_initial=w_guodu(:,3);

h0=[s_initial;x_initial;y_initial;z_initial;wx_initial;wy_initial;wz_initial;qr_initial;wr_initial];   

[t,h]=ode45(@f7_5_1,[0 40],h0,[],n,Jr,taor,J1,J2,J3,J4,J5,J6,kq1,kq2,kq3,kq4,kq5,kq6,Kw1,Kw2,Kw3,Kw4,Kw5,Kw6,M,b);       
s=h(:,1:n);   
x=h(:,(n+1):(2*n));
y=h(:,(2*n+1):(3*n));
z=h(:,(3*n+1):(4*n));
wx=h(:,(4*n+1):(5*n));
wy=h(:,(5*n+1):(6*n));
wz=h(:,(6*n+1):(7*n));
qrs=h(:,(7*n+1));
qrx=h(:,(7*n+2));
qry=h(:,(7*n+3));
qrz=h(:,(7*n+4));
wrx=h(:,(7*n+5));
wry=h(:,(7*n+6));
wrz=h(:,(7*n+7));

plot(t,s);  
hold on;
plot(t,qrs,'-.');    % 画出时变qr的第一维曲线
grid on;    
xlabel('t');
ylabel('s');
legend('刚体1','刚体2','刚体3','刚体4','刚体5','刚体6','qrs');
title('算法7.5-n个刚体的姿态四元数的标量部分');

figure;
plot(t,x);
hold on;
plot(t,qrx,'-.');    % 画出时变qr的第二维曲线
grid on;    
xlabel('t');
ylabel('x');
legend('刚体1','刚体2','刚体3','刚体4','刚体5','刚体6','qrx');
title('算法7.5-n个刚体的姿态四元数的向量部分的第一分量');

figure;
plot(t,y);
hold on;
plot(t,qry,'-.');    % 画出时变qr的第三维曲线
grid on;    
xlabel('t');
ylabel('y');
legend('刚体1','刚体2','刚体3','刚体4','刚体5','刚体6','qry');
title('算法7.5-n个刚体的姿态四元数的向量部分的第二分量');

figure;
plot(t,z);
hold on;
plot(t,qrz,'-.');    % 画出时变qr的第四维曲线
grid on;    
xlabel('t');
ylabel('z');
legend('刚体1','刚体2','刚体3','刚体4','刚体5','刚体6','qrz');
title('算法7.5-n个刚体的姿态四元数的向量部分的第三分量');

figure;
plot(t,wx);
hold on;
plot(t,wrx,'-.');    % 画出时变wr的第一维曲线
grid on;    
xlabel('t');
ylabel('wx');
legend('刚体1','刚体2','刚体3','刚体4','刚体5','刚体6','wrx');
title('算法7.5-n个刚体的角速度的第一分量');

figure;
plot(t,wy);
hold on;
plot(t,wry,'-.');    % 画出时变wr的第二维曲线
grid on;    
xlabel('t');
ylabel('wy');
legend('刚体1','刚体2','刚体3','刚体4','刚体5','刚体6','wry');
title('算法7.5-n个刚体的角速度的第二分量');

figure;
plot(t,wz);
hold on;
plot(t,wrz,'-.');    % 画出时变wr的第三维曲线
grid on;    
xlabel('t');
ylabel('wz');
legend('刚体1','刚体2','刚体3','刚体4','刚体5','刚体6','wrz');
title('算法7.5-n个刚体的角速度的第三分量');

