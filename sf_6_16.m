A=5*[0 1 0 1 0 0;1 0 1 0 1 0;0 1 0 0 0 1;1 0 0 0 1 0;0 1 0 1 0 1;0 0 1 0 1 0];  %5
B=10*[0 1 0 1 0 0;1 0 1 0 1 0;0 1 0 0 0 1;1 0 0 0 1 0;0 1 0 1 0 1;0 0 1 0 1 0];  %10
DA=diag(sum(A,2));   
DB=diag(sum(B,2));
LA=DA-A;               % �����ƺ�û�õ�
LB=DB-B;
num=6;

J(:,:,1)=[1 0.1 0.1;0.1 0.1 0.1;0.1 0.1 0.9];
J(:,:,2)=[1.5 0.2 0.3;0.2 0.9 0.4;0.3 0.4 2.0];
J(:,:,3)=[0.8 0.1 0.2;0.1 0.7 0.3;0.2 0.3 1.1];
J(:,:,4)=[1.2 0.3 0.7;0.3 0.9 0.2;0.7 0.2 1.4];
J(:,:,5)=[0.9 0.15 0.3;0.15 1.2 0.4;0.3 0.4 1.2];
J(:,:,6)=[1.1 0.35 0.45;0.35 1.0 0.5;0.45 0.5 1.3];
Dg=2*eye(3);
qr=quaternion(1 ,0, 0, 0);
Kg=1;
f=0.001;
t=0:f:30;

%��ʼ��q
q_init=normalize(quaternion(rand(num,4)));
q_init_quatParts=compact(q_init);
s_init=q_init_quatParts(:,1);  %����
x_init=q_init_quatParts(:,2);  
y_init=q_init_quatParts(:,3);
z_init=q_init_quatParts(:,4);
%��ʼ��w
omega=rand(num,3);
for i=1:num
    omega(i,:)=omega(i,:)/norm(omega(i,:));
end
wbo_initial_quatParts=[zeros(num,1),omega];   % wbo_initial_quatParts(n*4)��n*1ά����Ԫ��ʸ��wbo_initial��ϵ������
wbo_initial=quaternion(wbo_initial_quatParts);    % wbo_initial��n*1ά����Ԫ��ʸ��
wx_initial=wbo_initial_quatParts(:,2);
wy_initial=wbo_initial_quatParts(:,3);
wz_initial=wbo_initial_quatParts(:,4);

w1=[wx_initial(1);wy_initial(1);wz_initial(1)];           % wi����3*1ά��������
w2=[wx_initial(2);wy_initial(2);wz_initial(2)];
w3=[wx_initial(3);wy_initial(3);wz_initial(3)];
w4=[wx_initial(4);wy_initial(4);wz_initial(4)];
w5=[wx_initial(5);wy_initial(5);wz_initial(5)];
w6=[wx_initial(6);wy_initial(6);wz_initial(6)];
w=[w1;w2;w3;w4;w5;w6];            % w����3n*1ά��������

q=q_init;
qv=quaternion(zeros(num),zeros(num),zeros(num),zeros(num));   % ����qvΪn*nά����Ԫ������(��Ԫ�������ٽ���compact���ᷢ�������С��Ť��)
sumqv=zeros(3*num,1); %ÿ�ε����������
qv_vector=zeros(3*num,num); %Kg����������� %ÿ�ε����������
sumomega=zeros(3*num,1); %ÿ�ε����������
sumtotal=zeros(3*num,1); %ÿ�ε����������
tao=zeros(3*num,1); %ÿ�ε����������
wtemp=zeros(num,4); %ÿ�ε����������
dw=zeros(3*num,1);
for n=1:length(t)
    sumqv=zeros(3*num,1); %ÿ�ε����������
    sumomega=zeros(3*num,1); %ÿ�ε����������
    sumtotal=zeros(3*num,1); %ÿ�ε����������
    dw=zeros(3*num,1);
    for k=1:num
        qs=conj(qr)*q(k,n);
        qs_quatParts=compact(qs); %ÿ�ε����������
        for j=1:num %ÿ�ε����������
            qv(j,k)=conj(q(j,n))*q(k,n); %�㷨aij�������
            qv_n=compact(qv(j,k));
            qv_vector(3*j-2:3*j,k)=transpose(qv_n(2:4));
        end
        
        %����������һ����%ÿ�ε����������
        for j=1:num
            sumqv(3*k-2:3*k,1)=sumqv(3*k-2:3*k,1)+A(k,j)*qv_vector(3*j-2:3*j,k); %qingling
            sumomega(3*k-2:3*k,1)=sumomega(3*k-2:3*k,1)+B(k,j)*(w(3*k-2:3*k,n)-w(3*j-2:3*j,n));%qingling
        end
        %���������
        
        sumtotal=sumomega+sumqv;  %ÿ�ε����������
        test1=Dg*w(3*k-2:3*k,n);
     
        tao(3*k-2:3*k,1)=cross(w(3*k-2:3*k,n),J(:,:,k)*w(3*k-2:3*k,n))-J(:,:,k)*sumtotal(3*k-2:3*k,1);
        dw(3*k-2:3*k,1)=-inv(J(:,:,k))*(cross(w(3*k-2:3*k,n),J(:,:,k)*w(3*k-2:3*k,n))-tao(3*k-2:3*k,1));%qingling
        w(3*k-2:3*k,n+1)=dw(3*k-2:3*k,1)*f+w(3*k-2:3*k,n); %����w
        
        wtemp(k,1)=0;           
        wtemp(k,2:4)=w(3*k-2:3*k,n);

        w_temp_quaternion=quaternion(wtemp);
        
        q_temp_quaternion=1/2* q(k,n)*w_temp_quaternion(k);
        q_temp_vector1=compact(q(:,n));
        q_temp_vector=compact(q_temp_quaternion);
        q(k,n+1)=q_temp_quaternion*f+q(k,n);
        for i=1:num
            q_vector_plot(4*i-3:4*i,n)=transpose(q_temp_vector1(i,:));
        end
    end
end

%����w
[row,col]=size(w);
[row1,col1]=size(q_vector_plot);
figure(1)
for i=1:num
    plot(t,w(3*i,1:col-1))
    hold on;
end
title('w����1')
hold off;

figure(2)
for i=1:num
    plot(t,w(3*i-1,1:col-1))
    hold on;
end
title('w����2')
hold off;

figure(3)
for i=1:num
    plot(t,w(3*i-2,1:col-1))
    hold on;
end
title('w����3')
hold off;

figure(4)
for i=1:num
    plot(t,q_vector_plot(4*i-3,1:col1));
    hold on;
end
title('q����1-����')
hold off;
figure(5)
for i=1:num
    plot(t,q_vector_plot(4*i-2,1:col1));
    hold on;
end
title('q����2-����1')
hold off;
figure(6)
for i=1:num
    plot(t,q_vector_plot(4*i-1,1:col1));
    hold on;
end
title('q����3-����2')
hold off;
figure(7)
for i=1:num
    plot(t,q_vector_plot(4*i,1:col1));
    hold on;
end
title('q����4-����3')
hold off;

