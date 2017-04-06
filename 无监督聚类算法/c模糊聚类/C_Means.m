%% ����C-Means�����㷨�Խ��з���
% �������е���C-means����ĺ�����ɣ�function [center,U,objFcn] = cm(data, cluster_n);  
% ��Ҫ��Ϊ�ĸ��������裺
% ��һ������������׼��
% �ڶ����ǳ�ʼ����������
% ��������ÿһ�ε��������м����µľ������ĺ;��뻮�־���
% ���Ĳ��ǽ����վ���õ��ľ��뻮�־���ת��������ǩ����

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
[center,U,objFcn] = cm(cmdata,2);  %����ģ��C��ֵ����
[maxU,label]=max(U);
cla_label=reshape(label,ra,ca);    %�����ǩ����
figure(3)
imshow(mat2gray(cla_label));
title('Ӱ�������')
hold on


