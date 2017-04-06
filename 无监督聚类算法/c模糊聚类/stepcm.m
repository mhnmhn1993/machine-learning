function [U_new, center, obj_fcn] = stepcm(data, U, cluster_n)
 % 通过迭代运算计算每一次迭代的聚类中心和距离矩阵（每一个像元离聚类中心的欧式距离）   

center = U*data./((ones(size(data, 2), 1)*sum(U'))'); % new center
dist = distcm(center, data);                          % fill the distance matrix
obj_fcn = sum(sum((dist.^2).*U));                     % objective function
%就是将距离矩阵转换成硬划分矩阵（就是每个像元开始是按距离隶属于每一个聚类中心）
%下面这样转换的目的是为每一次聚类求得聚类中心打基础
col_sum = min(dist);           % 计算每次迭代后的距离划分矩阵，然后求其每一列的最小值（说明离某个聚类中心最近）
m = cluster_n;
temp = repmat(col_sum, m, 1);  %repmat函数就是按指定的值（m*1）进行复制
U = dist;
U(U>temp)=0;                   %对应矩阵比较小于的就赋值为零
U(U==temp)=1;                  %对应矩阵比较大于的就赋值为1
U_new = U;