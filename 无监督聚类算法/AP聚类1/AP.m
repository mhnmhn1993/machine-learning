function idx = AP(S)
N = size(S,1);
%初始化为0
A=zeros(N,N);%Availabilities
R=zeros(N,N); %Responsibility
 
lam=0.9; % Set damping factor
%前后两次聚类结果相同的累计次数，初始为-1
same_time = -1;
%最多进行10000次循环退出
for iter=1:10000
    % Compute responsibilities
    Rold=R;
    AS=A+S;
    %找到每行中最大的值max(a(i,k)+s(i,k)),Y是一列最大值，I是序号
    [Y,I]=max(AS,[],2);
    %每行的最大值设为-1000，不会影响第二大值的寻找
    for i=1:N
        AS(i,I(i))=-1000;
    end
    %找到每行中第二大的值Y2，计算最大值位置的r(i,k)用
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
    %竖着复制
    A=repmat(sum(Rp,1),[N,1])-Rp;
    dA=diag(A);
    A=min(A,0);
    for k=1:N
        A(k,k)=dA(k);
    end;
    A=(1-lam)*A+lam*Aold; % Dampen availabilities
 
    if(same_time == -1)%第一次
        E=R+A;
        [~, idx_old] = max(E,[],2);
        same_time = 0;
    else
        E=R+A;
        [~, idx] = max(E,[],2);%计算得到新的idx
 
        if(sum(abs(idx_old-idx)) == 0)%前后两次迭代得到的分类结果完全相同
            same_time = same_time + 1;%次数加一
            if(same_time == 10)%分类结果稳定在当前结果10次
                %输出循环次数
                disp(num2str(iter));
                break;
            end
        end
        %尚未完全相同，继续迭代
        idx_old = idx;
    end
end

E=R+A;
[~, idx] = max(E,[],2);


