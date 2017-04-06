%% 基于C-Means聚类算法对进行分类
% 主程序中调用C-means聚类的函数完成：function [center,U,objFcn] = cm(data, cluster_n);  
% 主要分为四个基本步骤：
% 第一步是输入数据准备
% 第二步是初始化聚类中心
% 第三步是每一次迭代过程中计算新的聚类中心和距离划分矩阵
% 第四步是将最终聚类得到的距离划分矩阵转换成类别标签矩阵

clear all;  clc;
close all;
I = imread('tm2005mask.jpg');
I = double(I);
[ra,ca,layers_num]=size(I);
TM1 = I(:,:,1); %Blue
TM2 = I(:,:,2); %Green
TM3 = I(:,:,3); %Red
%TM4 = I(:,:,4); %NIR
%TM5 = I(:,:,5);
%TM6 = I(:,:,6);
%TM7 = I(:,:,7);
feature_layers=cat(3,TM3,TM1,TM2);
cmdata=[];
celldata=cell(3,1);
for i=1:3
    celldata{i,1}=feature_layers(:,:,i);
    cmdata(:,i)=celldata{i,1}(:);
end    
[center,U,objFcn] = cm(cmdata,2);  %进行模糊C均值聚类
[maxU,label]=max(U);
cla_label=reshape(label,ra,ca);    %分类标签矩阵
figure(3)
imshow(mat2gray(cla_label));
title('影像分类结果')
hold on


