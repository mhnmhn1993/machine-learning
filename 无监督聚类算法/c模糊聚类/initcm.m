function U = initcm(cluster_n, data_n)
%  ��ʼ��һ�����U������ת����Ӳ���־���Ϊ����ȷ����ʼ����������׼��
U = rand(cluster_n, data_n);
col_sum = min(U); % �Բ������������ÿ������Сֵ
m = cluster_n;
temp = repmat(col_sum, m, 1);  %  repmat�������ǰ�ָ����ֵ��m*1�����и���
U(U>temp)=0;    %��Ӧ����Ƚ�С�ڵľ͸�ֵΪ��
U(U==temp)=1;   %��Ӧ����Ƚϴ��ڵľ͸�ֵΪ1
