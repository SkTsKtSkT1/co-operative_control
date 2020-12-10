function dhdt=f7_2(t,h,n,A,LB)
dhdt=zeros(7*n,1);
s=h(1:n);              % sΪh(7n*1ά������)�ĵ�һ������
x=h((n+1):(2*n));      % xΪh(7n*1ά������)�ĵڶ������Σ���������
y=h((2*n+1):(3*n));
z=h((3*n+1):(4*n));
wx=h((4*n+1):(5*n));
wy=h((5*n+1):(6*n));
wz=h((6*n+1):(7*n));

w1=[wx(1);wy(1);wz(1)];           % wi����3*1ά��������
w2=[wx(2);wy(2);wz(2)];
w3=[wx(3);wy(3);wz(3)];
w4=[wx(4);wy(4);wz(4)];
w5=[wx(5);wy(5);wz(5)];
w6=[wx(6);wy(6);wz(6)];
w=[w1;w2;w3;w4;w5;w6];            % w����3n*1ά��������
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
q=normalize(quaternion(s,x,y,z));                    % qΪn*1ά����Ԫ��ʸ��
q_deltaall=normalize(quaternion(rand(n,4)));       %����q_delta�ĳ�ʼֵ
wbo=normalize(quaternion(zeros(n,1),wx,wy,wz));      % wboΪn*1ά����Ԫ��ʸ��-----��һ���Ƿ���Ҫ���ں���ĳ��������
w_delta=normalize(quaternion(zeros(n,1),w_delta1,w_delta2,w_delta3));
w_deltaall=[w_delta1;w_delta2;w_delta3];
%----------------------------------------------------------------------------------------------------------------------------------------
dq=0.5*(q.*wbo);             % dq��n*1ά����Ԫ��ʸ��
dq_quatParts=compact(dq);    % dq_quatParts(n*4ά����)��dq��ϵ������
dq_deltaall=(q_deltaall.*w_delta);  
dq_deltaall_quatParts=compact(dq_deltaall);


%---------------------------------------------------------------------------------------------------------------------------------------

qv=quaternion(zeros(n),zeros(n),zeros(n),zeros(n));   % ����qvΪn*nά����Ԫ������(��Ԫ�������ٽ���compact���ᷢ�������С��Ť��)
qv_vectorParts=zeros(3*n,n);         % qv_vectorParts��3n*nά�ľ��󣬸þ�����Ի���Ϊ���3*1ά�����������ǵİڷŹ���Ϊÿ��n������������n��
for j=1:n
    for i=1:n
        qv(j,i)=conj(q(j)).*(dq_deltaall(j)).*conj(dq_deltaall(i)).*q(i);        % qv�����ÿ��Ԫ����qjconj��qi����̬��֮����д��qv(j,i)��Ϊ�����ļ��еĶ�����һ��
        qv_now_quatParts=compact(qv(j,i));   % qv_now_quatParts(1*4ά������)ֻ��һ�����ɱ�����ÿһ��ѭ��ȡֵ����һ��������û�ж�j,i���м������
        qv_vectorParts((3*j-2):(3*j),i)=transpose(qv_now_quatParts(2:4));  % ����ʵ��ʾ����qv��Ԫ��������ÿ����Ԫ��Ԫ�ص�3*1ά��������(ȡ��֮��Ҫ����ת��)
    end
end

complex=zeros(3*n,1);  % complex��3n*1ά����������װ����n������ƴ�ӵ�3*1ά�����������е�ÿ����������ʾ����i=1:nʱ���ļ�P86��Բ��ֱ�����n��3*1ά��������
for i=1:n
    for j=1:n
        complex((3*i-2):(3*i))=complex((3*i-2):(3*i))+A(i,j)*qv_vectorParts((3*j-2):(3*j),i);   % ���ﲻ�üӸ���
    end
end
%----------------------------------------------------------------------------------------------------------------------------------------        
dw=-(+kron(LB,eye(3))*(w+w_deltaall)...
    +complex);
%----------------------------------------------------------------------------------------------------------------------------------------
dhdt(1:n)=dq_quatParts(:,1);                % д��dhdt(:,1:n)����dhdt(1:n,:)�ᱨ��
dhdt((n+1):(2*n))=dq_quatParts(:,2);
dhdt((2*n+1):(3*n))=dq_quatParts(:,3);
dhdt((3*n+1):(4*n))=dq_quatParts(:,4);
dhdt((4*n+1):(5*n))=dw([1,4,7,10,13,16]);
dhdt((5*n+1):(6*n))=dw([2,5,8,11,14,17]);
dhdt((6*n+1):(7*n))=dw([3,6,9,12,15,18]);



