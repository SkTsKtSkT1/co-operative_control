% �㷨6.1
% ��6.5-P92-����ͼ6.2(���b)
% ��ʵmv6_1_2��f6_1_2����û��Ҫ����д.m�ļ���ֻ��Ҫ��6_1_1�Ļ������޸ļ����������ɣ�����Ϊ�˸����ٵ�������������л���������д��.m�ļ�
% ��Ϊ�ļ��еķ���ͼ����qi,wi��ʼֵ��ѡȡ�����Ͼ���0~1֮�䣬���Ա�.m�ļ��ж����õ�rand�������

clear 
clc
close all

n=6;    
A=[0 5 0 5 0 0            % ����ͼ6.2(���b)
   5 0 5 0 0 0
   0 5 0 0 0 0
   5 0 0 0 5 0
   0 0 0 5 0 5
   0 0 0 0 5 0];          
B=[0 10 0 10 0 0
   10 0 10 0 0 0
   0 10 0 0 0 0
   10 0 0 0 10 0
   0 0 0 10 0 10
   0 0 0 0 10 0];      % �ٶ�����ͼGnA��GnB����ͬ��
DA=diag(sum(A,2));   
DB=diag(sum(B,2));
LA=DA-A;               % �����ƺ�û�õ�
LB=DB-B;
     
J1=[1 0.1 0.1;0.1 0.1 0.1;0.1 0.1 0.9];    % 6������Ĺ���Ť��
J2=[1.5 0.2 0.3;0.2 0.9 0.4;0.3 0.4 2.0];
J3=[0.8 0.1 0.2;0.1 0.7 0.3;0.2 0.3 1.1];
J4=[1.2 0.3 0.7;0.3 0.9 0.2;0.7 0.2 1.4];
J5=[0.9 0.15 0.3;0.15 1.2 0.4;0.3 0.4 1.2];
J6=[1.1 0.35 0.45;0.35 1.0 0.5;0.45 0.5 1.3];

kG=0;             % ����ͼ6.2(���b)
DG1=2*eye(3);
DG2=2*eye(3);
DG3=2*eye(3);
DG4=2*eye(3);
DG5=2*eye(3);
DG6=2*eye(3);

qr=quaternion(1,0,0,0);        % ����ʱ�����׼��̬qr

q_initial=normalize(quaternion(rand(n,4)));   % q_initial��n*1ά����Ԫ��ʸ��
q_initial_quatParts=compact(q_initial);    % q_initial_quatParts(n*4)��n*1ά����Ԫ��ʸ��q_initial��ϵ������
s_initial=q_initial_quatParts(:,1);
x_initial=q_initial_quatParts(:,2);
y_initial=q_initial_quatParts(:,3);
z_initial=q_initial_quatParts(:,4);

guodu=rand(n,3);
for i=1:n
    guodu(i,:)=guodu(i,:)/norm(guodu(i,:));
end
wbo_initial_quatParts=[zeros(n,1),guodu];   % wbo_initial_quatParts(n*4)��n*1ά����Ԫ��ʸ��wbo_initial��ϵ������
wbo_initial=quaternion(wbo_initial_quatParts);    % wbo_initial��n*1ά����Ԫ��ʸ��(ps:������ֱ�Ӱ�zeros(n,1),rand(n,3)д��quaternion(...)����bug)
wx_initial=wbo_initial_quatParts(:,2);
wy_initial=wbo_initial_quatParts(:,3);
wz_initial=wbo_initial_quatParts(:,4);

h0=[s_initial;x_initial;y_initial;z_initial;wx_initial;wy_initial;wz_initial];   % �������շ����ʱ���ע�ĸ����7����Ϣ�ĳ�ʼֵ,h0Ϊ7n*1ά��

[t,h]=ode45(@f6_1_2,[0 30],h0,[],n,J1,J2,J3,J4,J5,J6,kG,DG1,DG2,DG3,DG4,DG5,DG6,qr,A,LB);       
s=h(:,1:n);   % �������Ƶ��ϱ���Ӧ����x=z(1:n,:);����ode45����������������ģ�һ����������Ӧһ��ʱ�̣����Լ����Լ���ʼ��ʱ�����õ��������������Ľ��Ҳ������������
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
legend('����1','����2','����3','����4','����5','����6');
title('�㷨6.1-���b-n���������̬��Ԫ���ı�������');

figure;
plot(t,x);
grid on;    
xlabel('t');
ylabel('x');
legend('����1','����2','����3','����4','����5','����6');
title('�㷨6.1-���b-n���������̬��Ԫ�����������ֵĵ�һ����');

figure;
plot(t,y);
grid on;    
xlabel('t');
ylabel('y');
legend('����1','����2','����3','����4','����5','����6');
title('�㷨6.1-���b-n���������̬��Ԫ�����������ֵĵڶ�����');

figure;
plot(t,x);
grid on;    
xlabel('t');
ylabel('z');
legend('����1','����2','����3','����4','����5','����6');
title('�㷨6.1-���b-n���������̬��Ԫ�����������ֵĵ�������');

figure;
plot(t,wx);
grid on;    
xlabel('t');
ylabel('wx');
legend('����1','����2','����3','����4','����5','����6');
title('�㷨6.1-���b-n������Ľ��ٶȵĵ�һ����');

figure;
plot(t,wy);
grid on;    
xlabel('t');
ylabel('wy');
legend('����1','����2','����3','����4','����5','����6');
title('�㷨6.1-���b-n������Ľ��ٶȵĵڶ�����');

figure;
plot(t,wz);
grid on;    
xlabel('t');
ylabel('wz');
legend('����1','����2','����3','����4','����5','����6');
title('�㷨6.1-���b-n������Ľ��ٶȵĵ�������');
