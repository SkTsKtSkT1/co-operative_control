% �㷨7.5
% ��7.3-P106
% ������ͼ7.2-P106��֪��vL={1},N1={4},N2={1},N3={2},N4={5},N5={2,6},N6={3}
% �ʴˣ�J1=N1��{r}={4,r};J2=N2={1};J3=N3={2};J4=N4={5};J5=N5={2,6};J6=N6={3};(ע������д��Ji�����Ǹ���Ť�أ�ֻ��֤�������õ�һ���м����)

clear 
clc
close all

n=6;
A=[0 0 0 1 0 0 1;      % ����ͼ7.2-P106��GnA��GnB����ͬ��
   1 0 0 0 0 0 0;
   0 1 0 0 0 0 0;
   0 0 0 0 1 0 0;
   0 1 0 0 0 1 0;
   0 0 1 0 0 0 0;
   0 0 0 0 0 0 0];
D=diag(sum(A, 2));           % ������Ⱦ���D
L=D-A;                       % L(n+1)*(n+1)         
M=L(1:n,1:n);                % L(n+1)*(n+1)ȥ����n+1�е�n+1��֮�����M
b=L(1:n,n+1);                % L(n+1)*(n+1)��n+1��ȥ�����·�һ��Ԫ�ز���b

J1=[1 0.1 0.1;0.1 0.1 0.1;0.1 0.1 0.9];    % 6������Ĺ���Ť��
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

q_initial=normalize(quaternion(rand(n,4)));   % q_initial��n*1ά����Ԫ��ʸ��
q_initial_quatParts=compact(q_initial);       % q_initial_quatParts(n*4)��n*1ά����Ԫ��ʸ��q_initial��ϵ������
s_initial=q_initial_quatParts(:,1);           % s_initialҲ������w_guoduһ��ֱ�ӹ�����λ���ľ��󣬶�������ͨ�����쵥λ��Ԫ��ʸ��Ȼ����ȡ��ϵ������
x_initial=q_initial_quatParts(:,2);
y_initial=q_initial_quatParts(:,3);
z_initial=q_initial_quatParts(:,4);

w_guodu=rand(n,3);
for i=1:n
    w_guodu(i,:)=w_guodu(i,:)/norm(w_guodu(i,:));
end
% wbo_initial_quatParts=[zeros(n,1),w_guodu];   % wbo_initial_quatParts(n*4)��n*1ά����Ԫ��ʸ��wbo_initial��ϵ������
% wbo_initial=quaternion(wbo_initial_quatParts);    % wbo_initial��n*1ά����Ԫ��ʸ��
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
plot(t,qrs,'-.');    % ����ʱ��qr�ĵ�һά����
grid on;    
xlabel('t');
ylabel('s');
legend('����1','����2','����3','����4','����5','����6','qrs');
title('�㷨7.5-n���������̬��Ԫ���ı�������');

figure;
plot(t,x);
hold on;
plot(t,qrx,'-.');    % ����ʱ��qr�ĵڶ�ά����
grid on;    
xlabel('t');
ylabel('x');
legend('����1','����2','����3','����4','����5','����6','qrx');
title('�㷨7.5-n���������̬��Ԫ�����������ֵĵ�һ����');

figure;
plot(t,y);
hold on;
plot(t,qry,'-.');    % ����ʱ��qr�ĵ���ά����
grid on;    
xlabel('t');
ylabel('y');
legend('����1','����2','����3','����4','����5','����6','qry');
title('�㷨7.5-n���������̬��Ԫ�����������ֵĵڶ�����');

figure;
plot(t,z);
hold on;
plot(t,qrz,'-.');    % ����ʱ��qr�ĵ���ά����
grid on;    
xlabel('t');
ylabel('z');
legend('����1','����2','����3','����4','����5','����6','qrz');
title('�㷨7.5-n���������̬��Ԫ�����������ֵĵ�������');

figure;
plot(t,wx);
hold on;
plot(t,wrx,'-.');    % ����ʱ��wr�ĵ�һά����
grid on;    
xlabel('t');
ylabel('wx');
legend('����1','����2','����3','����4','����5','����6','wrx');
title('�㷨7.5-n������Ľ��ٶȵĵ�һ����');

figure;
plot(t,wy);
hold on;
plot(t,wry,'-.');    % ����ʱ��wr�ĵڶ�ά����
grid on;    
xlabel('t');
ylabel('wy');
legend('����1','����2','����3','����4','����5','����6','wry');
title('�㷨7.5-n������Ľ��ٶȵĵڶ�����');

figure;
plot(t,wz);
hold on;
plot(t,wrz,'-.');    % ����ʱ��wr�ĵ���ά����
grid on;    
xlabel('t');
ylabel('wz');
legend('����1','����2','����3','����4','����5','����6','wrz');
title('�㷨7.5-n������Ľ��ٶȵĵ�������');

