function idx = AP(S)
N = size(S,1);
%��ʼ��Ϊ0
A=zeros(N,N);%Availabilities
R=zeros(N,N); %Responsibility
 
lam=0.9; % Set damping factor
%ǰ�����ξ�������ͬ���ۼƴ�������ʼΪ-1
same_time = -1;
%������10000��ѭ���˳�
for iter=1:10000
    % Compute responsibilities
    Rold=R;
    AS=A+S;
    %�ҵ�ÿ��������ֵmax(a(i,k)+s(i,k)),Y��һ�����ֵ��I�����
    [Y,I]=max(AS,[],2);
    %ÿ�е����ֵ��Ϊ-1000������Ӱ��ڶ���ֵ��Ѱ��
    for i=1:N
        AS(i,I(i))=-1000;
    end
    %�ҵ�ÿ���еڶ����ֵY2���������ֵλ�õ�r(i,k)��
    [Y2,~]=max(AS,[],2);
    R=S-repmat(Y,[1,N]);
    for i=1:N
        R(i,I(i))=S(i,I(i))-Y2(i);
    end
    R=(1-lam)*R+lam*Rold; % Dampen responsibilities
    % Compute availabilities
    Aold=A;
    Rp=max(R,0);
    for k=1:N
        Rp(k,k)=R(k,k);
    end
    %���Ÿ���
    A=repmat(sum(Rp,1),[N,1])-Rp;
    dA=diag(A);
    A=min(A,0);
    for k=1:N
        A(k,k)=dA(k);
    end;
    A=(1-lam)*A+lam*Aold; % Dampen availabilities
 
    if(same_time == -1)%��һ��
        E=R+A;
        [~, idx_old] = max(E,[],2);
        same_time = 0;
    else
        E=R+A;
        [~, idx] = max(E,[],2);%����õ��µ�idx
 
        if(sum(abs(idx_old-idx)) == 0)%ǰ�����ε����õ��ķ�������ȫ��ͬ
            same_time = same_time + 1;%������һ
            if(same_time == 10)%�������ȶ��ڵ�ǰ���10��
                %���ѭ������
                disp(num2str(iter));
                break;
            end
        end
        %��δ��ȫ��ͬ����������
        idx_old = idx;
    end
end

E=R+A;
[~, idx] = max(E,[],2);


