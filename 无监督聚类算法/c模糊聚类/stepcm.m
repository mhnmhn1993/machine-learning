function [U_new, center, obj_fcn] = stepcm(data, U, cluster_n)
 % ͨ�������������ÿһ�ε����ľ������ĺ;������ÿһ����Ԫ��������ĵ�ŷʽ���룩   

center = U*data./((ones(size(data, 2), 1)*sum(U'))'); % new center
dist = distcm(center, data);                          % fill the distance matrix
obj_fcn = sum(sum((dist.^2).*U));                     % objective function
%���ǽ��������ת����Ӳ���־��󣨾���ÿ����Ԫ��ʼ�ǰ�����������ÿһ���������ģ�
%��������ת����Ŀ����Ϊÿһ�ξ�����þ������Ĵ����
col_sum = min(dist);           % ����ÿ�ε�����ľ��뻮�־���Ȼ������ÿһ�е���Сֵ��˵����ĳ���������������
m = cluster_n;
temp = repmat(col_sum, m, 1);  %repmat�������ǰ�ָ����ֵ��m*1�����и���
U = dist;
U(U>temp)=0;                   %��Ӧ����Ƚ�С�ڵľ͸�ֵΪ��
U(U==temp)=1;                  %��Ӧ����Ƚϴ��ڵľ͸�ֵΪ1
U_new = U;