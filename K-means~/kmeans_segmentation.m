    function kmeans_segmentation()
    clear;close all;clc;
    %% 读取测试图像
    im = imread('001.jpg');
    imshow(im), title('Imput image');  %%转换图像的颜色空间得到样本
    cform = makecform('srgb2lab');
    lab = applycform(im,cform);
    ab = double(lab(:,:,2:3));
    nrows = size(lab,1); ncols = size(lab,2);
    X = reshape(ab,nrows*ncols,2)';
          figure, scatter(X(1,:)',X(2,:)',3,'filled'),title('image 2');  box on; %显示颜色空间转换后的二维样本空间分布
    %% 对样本空间进行Kmeans聚类
    k = 5; % 聚类个数
    max_iter = 100; %最大迭代次数
    [centroids, labels] = run_kmeans(X, k, max_iter); 

    %% 显示聚类分割结果
    figure, scatter(X(1,:)',X(2,:)'3,labels,'filled'),title('image 3'); %显示二维样本空间聚类效果
    hold on; scatter(centroids(1,:),centroids(2,:), 60,'r','filled')
    hold on; scatter(centroids(1,:),centroids(2,:),30,'g','filled')
    box on; hold off;
    %print -dpdf 2D2.pdf

    pixel_labels = reshape(labels,nrows,ncols);
    rgb_labels = label2rgb(pixel_labels);
    figure, imshow(rgb_labels), title('Segmented Image');
    %print -dpdf Seg.pdf
    end