%% 基于C-Means聚类算法:cm函数function [center, U, obj_fcn] = cm(data, cluster_n)
% 编程思想来源于FCM编程过程
% 输入参量说明：data是输入数据，行是像元个数，列是维数即聚类类别数，cluster_n是需要聚类的类数
% 输出参量说明：center最终迭代的聚类中心，U是最终迭代后的硬划分矩阵，obj_fcn是迭代过程中的目标函数             
function [center, U, obj_fcn] = cm(data, cluster_n)
data_n = size(data, 1);
in_n = size(data, 2);

% 设置最大迭代次数，迭代结束的阈值，显示每次迭代过程
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