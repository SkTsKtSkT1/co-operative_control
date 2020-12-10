% ������ͼ7.2-P106��֪��vL={1},N1={4},N2={1},N3={2},N4={5},N5={2,6},N6={3}
% �ʴˣ�J1=N1��{r}={4,r};J2=N2={1};J3=N3={2};J4=N4={5};J5=N5={2,6};J6=N6={3};(ע������д��Ji�����Ǹ���Ť�أ�ֻ��֤�������õ�һ���м����)

function dhdt=f7_5_1(t,h,n,Jr,taor,J1,J2,J3,J4,J5,J6,kq1,kq2,kq3,kq4,kq5,kq6,Kw1,Kw2,Kw3,Kw4,Kw5,Kw6,M,b)
dhdt=zeros((7*n+7),1);
s=h(1:n);                 % n���������̬�ĵ�һά
x=h((n+1):(2*n));         % n���������̬�ĵڶ�ά
y=h((2*n+1):(3*n));       % n���������̬�ĵ���ά
z=h((3*n+1):(4*n));       % n���������̬�ĵ���ά
wx=h((4*n+1):(5*n));      % n������Ľ��ٶȵĵ�һά
wy=h((5*n+1):(6*n));      % n������Ľ��ٶȵĵڶ�ά
wz=h((6*n+1):(7*n));      % n������Ľ��ٶȵĵ���ά
qr=h((7*n+1):(7*n+4));    % ʱ���׼��̬  
wr=h((7*n+5):(7*n+7));    % ʱ���׼���ٶ�

w1=[wx(1);wy(1);wz(1)];   % wi����3*1ά��������
w2=[wx(2);wy(2);wz(2)];
w3=[wx(3);wy(3);wz(3)];
w4=[wx(4);wy(4);wz(4)];
w5=[wx(5);wy(5);wz(5)];
w6=[wx(6);wy(6);wz(6)];
w=[w1;w2;w3;w4;w5;w6];            % w����3n*1ά��������

q=quaternion(s,x,y,z);                    % qΪn*1ά����Ԫ��ʸ��
wbo=quaternion(zeros(n,1),wx,wy,wz);      % wboΪn*1ά����Ԫ��ʸ��

spai1=(conj(q(4)).*q(1)).*(conj(quaternion(transpose(qr))).*q(1));   % spai1Ϊһ����Ԫ������ͬ
spai2=conj(q(1)).*q(2);
spai3=conj(q(2)).*q(3);
spai4=conj(q(5)).*q(4);
spai5=(conj(q(2)).*q(5)).*(conj(q(6)).*q(5));
spai6=conj(q(3)).*q(4);

spai1_quatParts=compact(spai1);   % spai1_quatPartsΪ1*4������������ͬ
spai2_quatParts=compact(spai2);
spai3_quatParts=compact(spai3);
spai4_quatParts=compact(spai4);
spai5_quatParts=compact(spai5);
spai6_quatParts=compact(spai6);

wsita1=(w1-w4)+(w1-wr);
wsita2=(w2-w1);
wsita3=(w3-w2);
wsita4=(w4-w5);
wsita5=(w5-w2)+(w5-w6);
wsita6=(w6-w3);

d1=-kq1*transpose(spai1_quatParts(2:4))-Kw1*wsita1;   % d1��3*1������������ͬ
d2=-kq2*transpose(spai2_quatParts(2:4))-Kw2*wsita2;   
d3=-kq3*transpose(spai3_quatParts(2:4))-Kw3*wsita3;  
d4=-kq4*transpose(spai4_quatParts(2:4))-Kw4*wsita4;   
d5=-kq5*transpose(spai5_quatParts(2:4))-Kw5*wsita5;   
d6=-kq6*transpose(spai6_quatParts(2:4))-Kw6*wsita6;  
d=[d1;d2;d3;d4;d5;d6];    % d��3n*1��������

%-------------------------------------------dwr��dqr΢�ַ��̵Ĺ���-----------------------------------------------------------------------------------
dwr=pinv(Jr)*[-cross(wr,Jr*wr)+taor];
dqr=0.5*quatmultiply(transpose(qr),transpose([0;wr]));     % �þ���ֱ�ӱ�ʾ��Ԫ���Ļ�����������������ʾ

%-------------------------------------------dq΢�ַ��̵Ĺ���------------------------------------------------------------------------------
dq=0.5*(q.*wbo);             % dq��n*1ά����Ԫ��ʸ��
dq_quatParts=compact(dq);    % dq_quatParts(n*4ά����)��dq��ϵ������

%--------------------------------------------dw΢�ַ��̵Ĺ���----------------------------------------------------------------------------
dw=pinv(kron(M,eye(3)))*(pinv(blkdiag(J1,J2,J3,J4,J5,J6))*d-kron(b,eye(3))*dwr);

%----------------------------------------------------------------------------------------------------------------------------------------
dhdt(1:n)=dq_quatParts(:,1);                % д��dhdt(:,1:n)����dhdt(1:n,:)�ᱨ��
dhdt((n+1):(2*n))=dq_quatParts(:,2);
dhdt((2*n+1):(3*n))=dq_quatParts(:,3);
dhdt((3*n+1):(4*n))=dq_quatParts(:,4);
dhdt((4*n+1):(5*n))=dw([1,4,7,10,13,16]);
dhdt((5*n+1):(6*n))=dw([2,5,8,11,14,17]);
dhdt((6*n+1):(7*n))=dw([3,6,9,12,15,18]);
dhdt((7*n+1):(7*n+4))=dqr;
dhdt((7*n+5):(7*n+7))=dwr;
