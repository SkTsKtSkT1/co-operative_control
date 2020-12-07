

A=1*[0 1 0 1 0 0;1 0 1 0 1 0;0 1 0 0 0 1;1 0 0 0 1 0;0 1 0 1 0 1;0 0 1 0 1 0];  %5
B=1*[0 1 0 1 0 0;1 0 1 0 1 0;0 1 0 0 0 1;1 0 0 0 1 0;0 1 0 1 0 1;0 0 1 0 1 0];  %10
DA=diag(sum(A,2));   
DB=diag(sum(B,2));
LA=DA-A;               % 后面似乎没用到
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
f=0.01;
t=0:f:30;

%初始的q
q_init=normalize(quaternion(rand(num,4)));
q_init_quatParts=compact(q_init);
s_init=q_init_quatParts(:,1);  %标量
x_init=q_init_quatParts(:,2);  
y_init=q_init_quatParts(:,3);
z_init=q_init_quatParts(:,4);
%初始的w
omega=rand(num,3);
for i=1:num
    omega(i,:)=omega(i,:)/norm(omega(i,:));
end
wbo_initial_quatParts=[zeros(num,1),omega];   % wbo_initial_quatParts(n*4)是n*1维的四元数矢量wbo_initial的系数矩阵
wbo_initial=quaternion(wbo_initial_quatParts);    % wbo_initial是n*1维的四元数矢量
wx_initial=wbo_initial_quatParts(:,2);
wy_initial=wbo_initial_quatParts(:,3);
wz_initial=wbo_initial_quatParts(:,4);

w1=[wx_initial(1);wy_initial(1);wz_initial(1)];           % wi均是3*1维的列向量
w2=[wx_initial(2);wy_initial(2);wz_initial(2)];
w3=[wx_initial(3);wy_initial(3);wz_initial(3)];
w4=[wx_initial(4);wy_initial(4);wz_initial(4)];
w5=[wx_initial(5);wy_initial(5);wz_initial(5)];
w6=[wx_initial(6);wy_initial(6);wz_initial(6)];
w=[w1;w2;w3;w4;w5;w6];            % w均是3n*1维的列向量

q=q_init;
qv=quaternion(zeros(num),zeros(num),zeros(num),zeros(num));   % 定义qv为n*n维的四元数数组(四元数数组再进行compact，会发生矩阵大小的扭曲)
sumqv=zeros(3*num,1); %每次迭代都会更新
qv_vector=zeros(3*num,num); %Kg后的向量部分 %每次迭代都会更新
sumomega=zeros(3*num,1); %每次迭代都会更新
sumtotal=zeros(3*num,1); %每次迭代都会更新
tao=zeros(3*num,1); %每次迭代都会更新
wtemp=zeros(num,4); %每次迭代都会更新
dw=zeros(3*num,1);
for n=1:length(t)
    sumqv=zeros(3*num,1); %每次迭代都会更新
    sumomega=zeros(3*num,1); %每次迭代都会更新
    sumtotal=zeros(3*num,1); %每次迭代都会更新
    dw=zeros(3*num,1);
    for k=1:num
        qs=conj(qr)*q(k,n);
        qs_quatParts=compact(qs); %每次迭代都会更新
        for j=1:num %每次迭代都会更新
            qv(j,k)=conj(q(j,n))*q(k,n); %算法aij后的向量
            qv_n=compact(qv(j,k));
            qv_vector(3*j-2:3*j,k)=transpose(qv_n(2:4));
        end
        
        %计算求和项第一部分%每次迭代都会更新
        for j=1:num
            sumqv(3*k-2:3*k,1)=sumqv(3*k-2:3*k,1)+A(k,j)*qv_vector(3*j-2:3*j,k); %qingling
            sumomega(3*k-2:3*k,1)=sumomega(3*k-2:3*k,1)+B(k,j)*(w(3*k-2:3*k,n)-w(3*j-2:3*j,n));%qingling
        end
        %求和项完整
        %求和项和控制量都好大 怀疑是求和的问题(⊙o⊙)…  -》貌似是A和B太大了
        sumtotal=sumomega+sumqv;  %每次迭代都会更新
        test1=Dg*w(3*k-2:3*k,n);

       
        tao(3*k-2:3*k,1)=-Kg*[transpose(qs_quatParts(1,2:4))]-Dg*w(3*k-2:3*k,n)-sumtotal(3*k-2:3*k,1);
%         test2=(cross(w(3*k-2:3*k,n),J(:,:,k)*w(3*k-2:3*k,n))-tao(3*k-2:3*k,1));
%         test3=-inv(J(:,:,k));
%         test3*test2 %%
        dw(3*k-2:3*k,1)=-inv(J(:,:,k))*(cross(w(3*k-2:3*k,n),J(:,:,k)*w(3*k-2:3*k,n))-tao(3*k-2:3*k,1));%qingling
        w(3*k-2:3*k,n+1)=dw(3*k-2:3*k,1)*f+w(3*k-2:3*k,n); %更新w
        
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

%绘制w
[row,col]=size(w);
[row1,col1]=size(q_vector_plot);
figure(1)
for i=1:num
    plot(t,w(3*i,1:col-1))
    hold on;
end
title('w分量1')
hold off;

figure(2)
for i=1:num
    plot(t,w(3*i-1,1:col-1))
    hold on;
end
title('w分量2')
hold off;

figure(3)
for i=1:num
    plot(t,w(3*i-2,1:col-1))
    hold on;
end
title('w分量3')
hold off;

figure(4)
for i=1:num
    plot(t,q_vector_plot(4*i-3,1:col1));
    hold on;
end
title('q分量1-标量')
hold off;
figure(5)
for i=1:num
    plot(t,q_vector_plot(4*i-2,1:col1));
    hold on;
end
title('q分量2-向量1')
hold off;
figure(6)
for i=1:num
    plot(t,q_vector_plot(4*i-1,1:col1));
    hold on;
end
title('q分量3-向量2')
hold off;
figure(7)
for i=1:num
    plot(t,q_vector_plot(4*i,1:col1));
    hold on;
end
title('q分量4-向量3')
hold off;

