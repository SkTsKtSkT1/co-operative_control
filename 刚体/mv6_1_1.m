% 算法6.1
% 例6.5-P92-无向图6.1(情况a)
% 因为文件中的仿真图像上qi,wi初始值的选取基本上均在0~1之间，所以本.m文件中都是用的rand随机函数

clear 
clc
close all

n=6;    
A=[0 5 0 5 0 0            % 无向图6.1(情况a)
   5 0 5 0 5 0
   0 5 0 0 0 5
   5 0 0 0 5 0
   0 5 0 5 0 5
   0 0 5 0 5 0];          
B=[0 10 0 10 0 0
   10 0 10 0 10 0
   0 10 0 0 0 10
   10 0 0 0 10 0
   0 10 0 10 0 10
   0 0 10 0 10 0];                     % 假定无向图GnA与GnB是相同的
DA=diag(sum(A,2));   
DB=diag(sum(B,2));
LA=DA-A;               % 后面似乎没用到
LB=DB-B;
     
J1=[1 0.1 0.1;0.1 0.1 0.1;0.1 0.1 0.9];    % 6个刚体的惯性扭矩
J2=[1.5 0.2 0.3;0.2 0.9 0.4;0.3 0.4 2.0];
J3=[0.8 0.1 0.2;0.1 0.7 0.3;0.2 0.3 1.1];
J4=[1.2 0.3 0.7;0.3 0.9 0.2;0.7 0.2 1.4];
J5=[0.9 0.15 0.3;0.15 1.2 0.4;0.3 0.4 1.2];
J6=[1.1 0.35 0.45;0.35 1.0 0.5;0.45 0.5 1.3];

kG=40;             % 无向图6.1(情况a)
DG1=2*eye(3);
DG2=2*eye(3);
DG3=2*eye(3);
DG4=2*eye(3);
DG5=2*eye(3);
DG6=2*eye(3);

qr=quaternion(1,0,0,0);        % 给出时不变基准姿态qr

q_initial=normalize(quaternion(rand(n,4)));   % q_initial是n*1维的四元数矢量
q_initial_quatParts=compact(q_initial);       % q_initial_quatParts(n*4)是n*1维的四元数矢量q_initial的系数矩阵
s_initial=q_initial_quatParts(:,1);  %标量
x_initial=q_initial_quatParts(:,2);  
y_initial=q_initial_quatParts(:,3);
z_initial=q_initial_quatParts(:,4);

guodu=rand(n,3);
for i=1:n
    guodu(i,:)=guodu(i,:)/norm(guodu(i,:));
end
wbo_initial_quatParts=[zeros(n,1),guodu];   % wbo_initial_quatParts(n*4)是n*1维的四元数矢量wbo_initial的系数矩阵
wbo_initial=quaternion(wbo_initial_quatParts);    % wbo_initial是n*1维的四元数矢量
wx_initial=wbo_initial_quatParts(:,2);
wy_initial=wbo_initial_quatParts(:,3);
wz_initial=wbo_initial_quatParts(:,4);

h0=[s_initial;x_initial;y_initial;z_initial;wx_initial;wy_initial;wz_initial];   % 这是最终仿真的时候关注的刚体的7个信息的初始值,h0为7n*1维度

[t,h]=ode45(@f6_1_1,[0 30],h0,[],n,J1,J2,J3,J4,J5,J6,kG,DG1,DG2,DG3,DG4,DG5,DG6,qr,A,LB);       
s=h(:,1:n);   % 从理论推导上本来应该是x=z(1:n,:);但是ode45是以行向量来计算的，一个行向量对应一个时刻，所以即便自己初始的时候设置的是列向量，最后的结果也会掰成行向量
x=h(:,(n+1):(2*n));
y=h(:,(2*n+1):(3*n));
z=h(:,(3*n+1):(4*n));
wx=h(:,(4*n+1):(5*n));
wy=h(:,(5*n+1):(6*n));
wz=h(:,(6*n+1):(7*n));

plot(t,s);    
grid on;    
xlabel('t');
ylabel('s');
legend('刚体1','刚体2','刚体3','刚体4','刚体5','刚体6');
title('算法6.1-情况a-n个刚体的姿态四元数的标量部分');

figure;
plot(t,x);
grid on;    
xlabel('t');
ylabel('x');
legend('刚体1','刚体2','刚体3','刚体4','刚体5','刚体6');
title('算法6.1-情况a-n个刚体的姿态四元数的向量部分的第一分量');

figure;
plot(t,y);
grid on;    
xlabel('t');
ylabel('y');
legend('刚体1','刚体2','刚体3','刚体4','刚体5','刚体6');
title('算法6.1-情况a-n个刚体的姿态四元数的向量部分的第二分量');

figure;
plot(t,x);
grid on;    
xlabel('t');
ylabel('z');
legend('刚体1','刚体2','刚体3','刚体4','刚体5','刚体6');
title('算法6.1-情况a-n个刚体的姿态四元数的向量部分的第三分量');

figure;
plot(t,wx);
grid on;    
xlabel('t');
ylabel('wx');
legend('刚体1','刚体2','刚体3','刚体4','刚体5','刚体6');
title('算法6.1-情况a-n个刚体的角速度的第一分量');

figure;
plot(t,wy);
grid on;    
xlabel('t');
ylabel('wy');
legend('刚体1','刚体2','刚体3','刚体4','刚体5','刚体6');
title('算法6.1-情况a-n个刚体的角速度的第二分量');

figure;
plot(t,wz);
grid on;    
xlabel('t');
ylabel('wz');
legend('刚体1','刚体2','刚体3','刚体4','刚体5','刚体6');
title('算法6.1-情况a-n个刚体的角速度的第三分量');
