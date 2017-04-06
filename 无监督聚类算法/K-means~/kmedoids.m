function label = kmedoids( data,k,start_data )
% kmedoids k中心点算法函数
% data 待聚类的数据集,每一行是一个样本数据点
% k 聚类个数
% start_data 聚类初始中心值,每一行为一个中心点,有cluster_n行
% class_idx 聚类结果,每个样本点标记的类别
% 初始化变量
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