    function kmeans_segmentation()
    clear;close all;clc;
    %% ��ȡ����ͼ��
    im = imread('001.jpg');
    imshow(im), title('Imput image');  %%ת��ͼ�����ɫ�ռ�õ�����
    cform = makecform('srgb2lab');
    lab = applycform(im,cform);
    ab = double(lab(:,:,2:3));
    nrows = size(lab,1); ncols = size(lab,2);
    X = reshape(ab,nrows*ncols,2)';
          figure, scatter(X(1,:)',X(2,:)',3,'filled'),title('image 2');  box on; %��ʾ��ɫ�ռ�ת����Ķ�ά�����ռ�ֲ�
    %% �������ռ����Kmeans����
    k = 5; % �������
    max_iter = 100; %����������
    [centroids, labels] = run_kmeans(X, k, max_iter); 

    %% ��ʾ����ָ���
    figure, scatter(X(1,:)',X(2,:)'3,labels,'filled'),title('image 3'); %��ʾ��ά�����ռ����Ч��
    hold on; scatter(centroids(1,:),centroids(2,:), 60,'r','filled')
    hold on; scatter(centroids(1,:),centroids(2,:),30,'g','filled')
    box on; hold off;
    %print -dpdf 2D2.pdf

    pixel_labels = reshape(labels,nrows,ncols);
    rgb_labels = label2rgb(pixel_labels);
    figure, imshow(rgb_labels), title('Segmented Image');
    %print -dpdf Seg.pdf
    end