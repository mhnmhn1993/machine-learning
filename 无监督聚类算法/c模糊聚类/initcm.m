function U = initcm(cluster_n, data_n)
%  初始化一个随机U，将其转换成硬划分矩阵，为后面确定初始聚类中心做准备
U = rand(cluster_n, data_n);
col_sum = min(U); % 对产生的随机矩阵每列求最小值
m = cluster_n;
temp = repmat(col_sum, m, 1);  %  repmat函数就是按指定的值（m*1）进行复制
U(U>temp)=0;    %对应矩阵比较小于的就赋值为零
U(U==temp)=1;   %对应矩阵比较大于的就赋值为1
