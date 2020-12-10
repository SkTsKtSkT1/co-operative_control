function dhdt=f7_2(t,h,n,A,LB)
dhdt=zeros(7*n,1);
s=h(1:n);              % s为h(7n*1维列向量)的第一个六段
x=h((n+1):(2*n));      % x为h(7n*1维列向量)的第二个六段，以下类似
y=h((2*n+1):(3*n));
z=h((3*n+1):(4*n));
wx=h((4*n+1):(5*n));
wy=h((5*n+1):(6*n));
wz=h((6*n+1):(7*n));

w1=[wx(1);wy(1);wz(1)];           % wi均是3*1维的列向量
w2=[wx(2);wy(2);wz(2)];
w3=[wx(3);wy(3);wz(3)];
w4=[wx(4);wy(4);wz(4)];
w5=[wx(5);wy(5);wz(5)];
w6=[wx(6);wy(6);wz(6)];
w=[w1;w2;w3;w4;w5;w6];            % w均是3n*1维的列向量
% w_delta1=[sin(t+2);sin(t+2);sin(t+2);sin(t+2);sin(t+2);sin(t+2)];
% w_delta2=[sin(2*t+1);sin(2*t+1);sin(2*t+1);sin(2*t+1);sin(2*t+1);sin(2*t+1)];
% w_delta3=[sin(4*t);sin(4*t);sin(4*t);sin(4*t);sin(4*t);sin(4*t)];
% 
% dw_delta1=[cos(t+2);cos(t+2);cos(t+2);cos(t+2);cos(t+2);cos(t+2)];
% dw_delta2=[cos(2*t+1);cos(2*t+1);cos(2*t+1);cos(2*t+1);cos(2*t+1);cos(2*t+1)];
% dw_delta3=[cos(4*t);cos(4*t);cos(4*t);cos(4*t);cos(4*t);cos(4*t)];

w_delta1=rand(n,1);
w_delta2=rand(n,1);
w_delta3=rand(n,1);


% dw1=[dw_delta1;dw_delta2;dw_delta3];
q=normalize(quaternion(s,x,y,z));                    % q为n*1维的四元数矢量
q_deltaall=normalize(quaternion(rand(n,4)));       %给出q_delta的初始值
wbo=normalize(quaternion(zeros(n,1),wx,wy,wz));      % wbo为n*1维的四元数矢量-----这一句是否需要放在后面某处？？？
w_delta=normalize(quaternion(zeros(n,1),w_delta1,w_delta2,w_delta3));
w_deltaall=[w_delta1;w_delta2;w_delta3];
%----------------------------------------------------------------------------------------------------------------------------------------
dq=0.5*(q.*wbo);             % dq是n*1维的四元数矢量
dq_quatParts=compact(dq);    % dq_quatParts(n*4维矩阵)是dq的系数矩阵
dq_deltaall=(q_deltaall.*w_delta);  
dq_deltaall_quatParts=compact(dq_deltaall);


%---------------------------------------------------------------------------------------------------------------------------------------

qv=quaternion(zeros(n),zeros(n),zeros(n),zeros(n));   % 定义qv为n*n维的四元数数组(四元数数组再进行compact，会发生矩阵大小的扭曲)
qv_vectorParts=zeros(3*n,n);         % qv_vectorParts是3n*n维的矩阵，该矩阵可以划分为多个3*1维列向量，它们的摆放规律为每行n个列向量，摆n行
for j=1:n
    for i=1:n
        qv(j,i)=conj(q(j)).*(dq_deltaall(j)).*conj(dq_deltaall(i)).*q(i);        % qv数组的每个元素是qjconj与qi的姿态误差，之所以写成qv(j,i)是为了与文件中的定义相一致
        qv_now_quatParts=compact(qv(j,i));   % qv_now_quatParts(1*4维列向量)只是一个过渡变量，每一次循环取值都不一样，所以没有对j,i进行检索设计
        qv_vectorParts((3*j-2):(3*j),i)=transpose(qv_now_quatParts(2:4));  % 这其实表示的是qv四元数数组中每个四元数元素的3*1维向量部分(取出之后要进行转置)
    end
end

complex=zeros(3*n,1);  % complex是3n*1维的列向量，装的是n个竖向拼接的3*1维列向量，其中的每个列向量表示的是i=1:n时，文件P86椭圆框分别代表的n个3*1维的列向量
for i=1:n
    for j=1:n
        complex((3*i-2):(3*i))=complex((3*i-2):(3*i))+A(i,j)*qv_vectorParts((3*j-2):(3*j),i);   % 这里不用加负号
    end
end
%----------------------------------------------------------------------------------------------------------------------------------------        
dw=-(+kron(LB,eye(3))*(w+w_deltaall)...
    +complex);
%----------------------------------------------------------------------------------------------------------------------------------------
dhdt(1:n)=dq_quatParts(:,1);                % 写成dhdt(:,1:n)或者dhdt(1:n,:)会报错
dhdt((n+1):(2*n))=dq_quatParts(:,2);
dhdt((2*n+1):(3*n))=dq_quatParts(:,3);
dhdt((3*n+1):(4*n))=dq_quatParts(:,4);
dhdt((4*n+1):(5*n))=dw([1,4,7,10,13,16]);
dhdt((5*n+1):(6*n))=dw([2,5,8,11,14,17]);
dhdt((6*n+1):(7*n))=dw([3,6,9,12,15,18]);



