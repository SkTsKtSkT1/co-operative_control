% 算法7.1
% 刚体姿态达到相对姿态保持并且最终角速度为零
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

Adelta=[1 2 3 4;2 3 4 1;3 4 1 2;4 1 2 3;2 3 4 5 ;4 5 6 3];
     
J1=[1 0.1 0.1;0.1 0.1 0.1;0.1 0.1 0.9];    % 6个刚体的惯性扭矩
J2=[1.5 0.2 0.3;0.2 0.9 0.4;0.3 0.4 2.0];
J3=[0.8 0.1 0.2;0.1 0.7 0.3;0.2 0.3 1.1];
J4=[1.2 0.3 0.7;0.3 0.9 0.2;0.7 0.2 1.4];
J5=[0.9 0.15 0.3;0.15 1.2 0.4;0.3 0.4 1.2];
J6=[1.1 0.35 0.45;0.35 1.0 0.5;0.45 0.5 1.3];

kG=50;             % 无向图6.1(情况a)
DG1=5*eye(3);
DG2=5*eye(3);
DG3=5*eye(3);
DG4=5*eye(3);
DG5=5*eye(3);
DG6=5*eye(3);

qr=quaternion(1,0,0,0);        % 给出时不变基准姿态qr

q_initial=normalize(quaternion(rand(n,4)));   % q_initial是n*1维的单位四元数矢量
% q_initial=normalize(quaternion(Adelta));
q_initial_quatParts=compact(q_initial);    % q_initial_quatParts(n*4)是n*1维的四元数矢量q_initial的系数矩阵
s_initial=q_initial_quatParts(:,1);
x_initial=q_initial_quatParts(:,2);
y_initial=q_initial_quatParts(:,3);
z_initial=q_initial_quatParts(:,4);


% q_delta_1=normalize(quaternion(rand(n,4))); 
% q_delta_2=normalize(quaternion(rand(n,4))); % q_delta_1是n*1维的单位四元数矢量,表示i与期望姿态的偏差，随机生成
% q_delta_3=normalize(quaternion(rand(n,4))); 
% q_delta_4=normalize(quaternion(rand(n,4))); 
% q_delta_5=normalize(quaternion(rand(n,4))); 
% q_delta_6=normalize(quaternion(rand(n,4))); 
% q_delta_1_quatParts=compact(q_delta_1);  
% q_delta_2_quatParts=compact(q_delta_2);  
% q_delta_3_quatParts=compact(q_delta_3);  
% q_delta_4_quatParts=compact(q_delta_4);  
% q_delta_5_quatParts=compact(q_delta_5);  
% q_delta_6_quatParts=compact(q_delta_6);  
% delta0=[q_delta_1_quatParts;q_delta_2_quatParts;q_delta_3_quatParts;q_delta_4_quatParts;q_delta_5_quatParts;q_delta_6_quatParts];

q_delta=normalize(quaternion(rand(n,4)));   % q_initial是n*1维的单位四元数矢量
q_delta_quatParts=compact(q_delta);    % q_initial_quatParts(n*4)是n*1维的四元数矢量q_initial的系数矩阵
s_delta=q_delta_quatParts(:,1);
x_delta=q_delta_quatParts(:,2);
y_delta=q_delta_quatParts(:,3);
z_delta=q_delta_quatParts(:,4);


wbo_initial_quatParts=[zeros(n,1),rand(n,3)];   % wbo_initial_quatParts(n*4)是n*1维的四元数矢量wbo_initial的系数矩阵
wbo_initial=normalize(quaternion(wbo_initial_quatParts));    % wbo_initial是n*1维的四元数矢量
wx_initial=wbo_initial_quatParts(:,2);
wy_initial=wbo_initial_quatParts(:,3);
wz_initial=wbo_initial_quatParts(:,4);

h0=[s_initial;x_initial;y_initial;z_initial;wx_initial;wy_initial;wz_initial];   % 这是最终仿真的时候关注的刚体的7个信息的初始值,h0为7n*1维度

[t,h]=ode45(@f7_1,[0 50],h0,[],n,J1,J2,J3,J4,J5,J6,kG,DG1,DG2,DG3,DG4,DG5,DG6,s_delta,x_delta,y_delta,z_delta,qr,A,LB);       
s=h(:,1:n);   % 从理论推导上本来应该是x=z(1:n,:);但是ode45是以行向量来计算的，一个行向量对应一个时刻，所以即便自己初始的时候设置的是列向量，最后的结果也会掰成行向量
x=h(:,(n+1):(2*n));
y=h(:,(2*n+1):(3*n));
z=h(:,(3*n+1):(4*n));
wx=h(:,(4*n+1):(5*n));
wy=h(:,(5*n+1):(6*n));
wz=h(:,(6*n+1):(7*n));

plot(t,s);    % 暂不关注刚体姿态四元数的标量部分
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
title('刚体的姿态四元数的向量部分的第一分量');

figure;
plot(t,y);
grid on;    
xlabel('t');
ylabel('y');
legend('刚体1','刚体2','刚体3','刚体4','刚体5','刚体6');
title('刚体的姿态四元数的向量部分的第二分量');

figure;
plot(t,x);
grid on;    
xlabel('t');
ylabel('z');
legend('刚体1','刚体2','刚体3','刚体4','刚体5','刚体6');
title('刚体的姿态四元数的向量部分的第三分量');

figure;
plot(t,wx);
grid on;    
xlabel('t');
ylabel('wx');
legend('刚体1','刚体2','刚体3','刚体4','刚体5','刚体6');
title('刚体的角速度的第一分量');

figure;
plot(t,wy);
grid on;    
xlabel('t');
ylabel('wy');
legend('刚体1','刚体2','刚体3','刚体4','刚体5','刚体6');
title('刚体的角速度的第二分量');

figure;
plot(t,wz);
grid on;    
xlabel('t');
ylabel('wz');
legend('刚体1','刚体2','刚体3','刚体4','刚体5','刚体6');
title('刚体的角速度的第三分量');

% figure;
% quat1=1-x(1,1)^2+y(1,1)^2+z(1,1)^2;
% lzh1=quaternion(quat1,x(1,1),y(1,1),z(1,1));
% 
% quat2=1-x(1,2)^2+y(1,2)^2+z(1,2)^2;
% lzh2=quaternion(quat2,x(1,2),y(1,2),z(1,2));
% 
% quat3=1-x(1,3)^2+y(1,3)^2+z(1,3)^2;
% lzh3=quaternion(quat3,x(1,3),y(1,3),z(1,3));
% 
% quat4=1-x(1,4)^2+y(1,4)^2+z(1,4)^2;
% lzh4=quaternion(quat4,x(1,4),y(1,4),z(1,4));
% 
% quat5=1-x(1,5)^2+y(1,5)^2+z(1,5)^2;
% lzh5=quaternion(quat5,x(1,5),y(1,5),z(1,5));
% 
% quat6=1-x(1,6)^2+y(1,6)^2+z(1,6)^2;
% lzh6=quaternion(quat6,x(1,6),y(1,6),z(1,6));
% 
% normalize(lzh1);
% normalize(lzh2);
% normalize(lzh3);
% normalize(lzh4);
% normalize(lzh5);
% normalize(lzh6);
% rotlzh1=rotmat(lzh1,'frame');
% rotlzh2=rotmat(lzh2,'frame');
% rotlzh3=rotmat(lzh3,'frame');
% rotlzh4=rotmat(lzh4,'frame');
% rotlzh5=rotmat(lzh5,'frame');
% rotlzh6=rotmat(lzh6,'frame');
% triplot(rotlzh1,'color','r');
% hold on;
% triplot(rotlzh2,'color','b');tirplot(rotlzh3,'color','g');triplot(rotlzh4,'color','c');
% triplot(rotlzh5,'color','y');triplot(rotlzh6,'color','m');
% 
% quat1_final=1-x(1800,1)^2+y(1800,1)^2+z(1800,1)^2;
% lzh1_final=quaternion(quat1_final,x(1800,1),y(1800,1),z(1800,1));
% 
% quat2_final=1-x(1800,2)^2+y(1800,2)^2+z(1800,2)^2;
% lzh2_final=quaternion(quat2_final,x(1800,2),y(1800,2),z(1800,2));
% 
% quat3_final=1-x(1800,3)^2+y(1800,3)^2+z(1800,3)^2;
% lzh3_final=quaternion(quat3_final,x(1800,3),y(1800,3),z(1800,3));
% 
% quat4_final=1-x(1800,4)^2+y(1800,4)^2+z(1800,4)^2;
% lzh4_final=quaternion(quat4_final,x(1800,4),y(1800,4),z(1800,4));
% 
% 
% quat5_final=1-x(1800,5)^2+y(1800,5)^2+z(1800,5)^2;
% lzh5_final=quaternion(quat5_final,x(1800,5),y(1800,5),z(1800,5));
% 
% quat6_final=1-x(1800,6)^2+y(1800,6)^2+z(1800,6)^2;
% lzh6_final=quaternion(quat6_final,x(1800,6),y(1800,6),z(1800,6));
% 
% normalize(lzh1_final);
% normalize(lzh2_final);
% normalize(lzh3_final);
% normalize(lzh4_final);
% normalize(lzh5_final);
% normalize(lzh6_final);
% 
% rotlzh1_final=rotmat(lzh1_final,'frame');
% rotlzh2_final=rotmat(lzh2_final,'frame');
% rotlzh3_final=rotmat(lzh3_final,'frame');
% rotlzh4_final=rotmat(lzh4_final,'frame');
% rotlzh5_final=rotmat(lzh5_final,'frame');
% rotlzh6_final=rotmat(lzh6_final,'frame');
% figure;
% triplot(rotlzh1_final,'color','r');
% hold on;
% triplot(rotlzh2_final,'color','b');triplot(rotlzh3_final,'color','g');trplot(rotlzh4_final,'color','c');
% triplot(rotlzh5_final,'color','y');triplot(rotlzh6_final,'color','m');
