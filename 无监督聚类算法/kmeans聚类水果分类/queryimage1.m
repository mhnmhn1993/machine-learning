
clear all
close all
I_rgb=imread('tm2005mask.jpg');
figure,imshow(I_rgb);title('原始图像');
C = makecform('srgb2lab');      
I_lab = applycform(I_rgb, C);
ab = double(I_lab(:,:,2:3));  
nrows = size(ab,1);
ncols = size(ab,2);
ab = reshape(ab,nrows*ncols,2);
nColors = 3;      %number of segementation regions
%[cluster_idx cluster_center] = kmeans(ab,nColors,'distance','sqEuclidean','Replicates',3); 
[cluster_idx cluster_center] = kmeans(ab,nColors,'emptyaction','singleton'); 
pixel_labels = reshape(cluster_idx,nrows,ncols);
figure();imshow(pixel_labels,[]), title('聚类结果');
segmented_images = cell(1,3);
rgb_label = repmat(pixel_labels,[1 1 3]);
for k = 1:nColors
    n=strcat('分割结果——区域',num2str(k))
    color = I_rgb;
    color(rgb_label ~= k) = 0;
    segmented_images{k} = color;
    figure(),imshow(segmented_images{k});
    title(n);
end
image1=rgb2gray(segmented_images{2}); %interesting image
figure,imshow(image1);
tt=graythresh(image1);
image2=im2bw(image1,tt);
figure, imshow(image2) ;

SE = strel('disk',5);
i2= imclose(image2,SE);
i2= imfill(i2,'hole');
figure, imshow(i2);
SE = strel('disk',6);
i2= imopen(i2,SE);
figure, imshow(i2);
max=0;
I_lab2=imread('fruit2.jpg');
figure, imshow(I_lab2) ;hold on
[B,L]=bwboundaries(i2,'noholes');
STATS=regionprops(L,'all');
[L,num] = bwlabel(i2); 
for t=1:num       % find the region with maximum Area
    if STATS(t).Area>max
        max=STATS(t).Area;
        max_num=t;
    end
end
for t=1:num
    if t~=max_num
        rectangle('Position',[STATS(t).BoundingBox],'EdgeColor','g');
    end
end
