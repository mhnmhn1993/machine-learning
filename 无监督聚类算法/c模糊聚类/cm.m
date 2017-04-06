%% ����C-Means�����㷨:cm����function [center, U, obj_fcn] = cm(data, cluster_n)
% ���˼����Դ��FCM��̹���
% �������˵����data���������ݣ�������Ԫ����������ά���������������cluster_n����Ҫ���������
% �������˵����center���յ����ľ������ģ�U�����յ������Ӳ���־���obj_fcn�ǵ��������е�Ŀ�꺯��             
function [center, U, obj_fcn] = cm(data, cluster_n)
data_n = size(data, 1);
in_n = size(data, 2);

% ����������������������������ֵ����ʾÿ�ε�������
options = [100;	    % max. number of iteration
		   1e-5;	% min. amount of improvement
		   1];	    % info display during iteration 
       
max_iter = options(1);		% Max. iteration
min_impro = options(2);		% Min. improvement
display = options(3);		% Display info or not

obj_fcn = zeros(max_iter, 1);	% Array for objective function
U = initcm(cluster_n, data_n);	% Initial fuzzy partition
% Main loop
for i = 1:max_iter,
	[U, center, obj_fcn(i)] = stepcm(data, U, cluster_n);
	if display, 
		fprintf('Iteration count = %d, obj. fcn = %f\n', i, obj_fcn(i));
	end
	% check termination condition
	if i > 1,
		if abs(obj_fcn(i) - obj_fcn(i-1)) < min_impro, break; end,
	end
end

iter_n = i;	% Actual number of iterations 
obj_fcn(iter_n+1:max_iter) = [];