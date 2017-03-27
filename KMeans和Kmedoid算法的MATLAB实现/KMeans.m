function [ accuracy,MIhat ] = KMeans( K,mode )

% Artificial Intelligence & Data Mining - KMeans & K-Medoids Clustering
% Author: Sophia_qing @ ZJU
% CreateTime: 2012-11-18
% Function: Clustering
%  -K: number of clusters
%  -mode: 
%   1: use kmeans cluster algorithm in matlab
%   2: k_medroid algorithm: use data points as k centers
%   3: k_means algorithm: use average as k centers

global N_features;
global N_samples;
global fea;
global gnd;

switch (mode)
    case 1 %call system function KMeans
        label = kmeans(fea,K);
        [label,accuracy] = cal_accuracy(gnd,label);
        
    case 2%use kmedroid method
        for testcase = 1:10% do 10 times to get rid of the influence from Initial_center
            K_center = Initial_center(fea,K); %select initial points randomly
            changed_label = N_samples;
            label = zeros(1,N_samples);
            iteration_times = 0;
            while changed_label~=0
                cls_label = cell(1,K);
                for i = 1: N_samples
                    for j = 1 : K
                        D(j) = dis(fea(i,:),K_center(j,:));
                    end
                    [~,label(i)] = min(D);
                    cls_label{label(i)} = [cls_label{label(i)} i];
                end
                changed_label = 0;
                cls_center = zeros(K,N_features);
                for i = 1 : K
                    cls_center(i,:) = mean(fea(cls_label{i},:));
                    D1 = [];
                    for j = 1:size(cls_label{i},2)%number of samples clsutered in i-th class
                        D1(j) = dis(cls_center(i,:),fea(cls_label{i}(j),:));
                    end
                    [~,min_ind] = min(D1);
                    if ~isequal(K_center(i,:),fea(cls_label{i}(min_ind),:))
                        K_center(i,:) = fea(cls_label{i}(min_ind),:);
                        changed_label = changed_label+1;
                    end
                end
                iteration_times = iteration_times+1;
            end
            [label,acc(testcase)] = cal_accuracy(gnd,label);
        end
        accuracy = max(acc);
        
    case 3%use k-means method
        for testcase = 1:10% do 10 times to get rid of the influence from Initial_center
            K_center = Initial_center(fea,K); %select initial points randomly
            changed_label = N_samples;
            label = zeros(1,N_samples);
            label_new = zeros(1,N_samples);
            while changed_label~=0
                cls_label = cell(1,K);
                changed_label = 0;
                for i = 1: N_samples
                    for j = 1 : K
                        D(j) = dis(fea(i,:),K_center(j,:));
                    end
                    [~,label_new(i)] = min(D);
                    if(label_new(i)~=label(i))
                        changed_label = changed_label+1;
                    end;
                    cls_label{label_new(i)} = [cls_label{label_new(i)} i];
                end
                label = label_new;
                
                for i = 1 : K  %recalculate k centroid
                    K_center(i,:) = mean(fea(cls_label{i},:));
                end
            end
             [label,acc(testcase)] = cal_accuracy(gnd,label);
        end
        accuracy = max(acc);
end

MIhat = MutualInfo(gnd,label);


    function center = Initial_center(X,K)
        rnd_Idx = randperm(N_samples,K);
        center = X(rnd_Idx,:);
    end

    function res = dis(X1,X2)
        res = norm(X1-X2);
    end

    function [res,acc] = cal_accuracy(gnd,estimate_label)
        res = bestMap(gnd,estimate_label);
        acc = length(find(gnd == res))/length(gnd);
    end
end

