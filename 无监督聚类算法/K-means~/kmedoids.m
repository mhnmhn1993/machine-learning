function label = kmedoids( data,k,start_data )
% kmedoids k���ĵ��㷨����
% data ����������ݼ�,ÿһ����һ���������ݵ�
% k �������
% start_data �����ʼ����ֵ,ÿһ��Ϊһ�����ĵ�,��cluster_n��
% class_idx ������,ÿ���������ǵ����
% ��ʼ������
n = length(data);
dist_temp1 = zeros(n,k);
dist_temp2 = zeros(n,k);
last = zeros(n,1);
a = 0;
b = 0;

if nargin==3
   centroid  = start_data;
else
   centroid  = data(randsample(n,k),:); 
end
for a = 1:k
    temp1 = ones(n,1)*centroid(a,:);
    dist_temp1(:,a) = sum((data-temp1).^2,2);    
end
[~,label] = min(dist_temp1,[],2);
while any(label~=last)
    for a = 1:k
        temp2 = ones(numel(data(label==a)),1);
        temp3 = data(label==a);    
        for b = 1:n        
            temp4 = temp2*data(b,:);
            temp5 = sum((temp3-temp4).^2,2);
            dist_temp2(b,a) = sum(temp5,1); 
        end
    end
    [~,centry_indx] = min(dist_temp2,[],1);
    last = label;    
    centroid  = data(centry_indx,:);
    for a = 1:k
        temp1 = ones(n,1)*centroid(a,:);
        dist_temp1(:,a) = sum((data-temp1).^2,2);    
    end
    [~,label] = min(dist_temp1,[],2);  
end

end