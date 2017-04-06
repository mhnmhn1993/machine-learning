clc;clear;clf;
img=imread('depth.png'); %读取深度图。
figure(1);imshow(img);title('原图');
Narrow = 0.25;
img_s=imresize(img,Narrow,'nearest'); %双线性插值 

figure(2); imshow(img_s);title('压缩后的图像');
[M,N]=size(img_s);            
[XX,YY]=meshgrid(1:N,1:M);
data=[XX(:) YY(:) double(img_s(:))*Narrow*0.5];%将深度图转换为3维数据集
data = double(data);
figure(3);plot3(data(:,1),data(:,2),data(:,3),'.');title('3维数据');
[n,~]=size(data);
%利用欧氏距离的相反数作为相似度矩阵
S = -squareform(pdist(data));
[row,col,v]=find(S);%找到不为0的元素
%v是一列n*(n-1)个元素的非0向量，排序之后取中位数
v=sortrows(v,1);
%med = v(floor(size(v,1)*0.005));
med = v(1)-95;
for i=1:n
    S(i,i)=med;
end
idx = AP(S);



figure(7);
hold on;
label = idx;
imgc = zeros(M,N,3);
for i=1:max(label)
  ii=find(label == i);
  if sum(data(ii,3))/size(ii,1)>240/8
      for j = 1:size(ii)
          k = ii(j);
        imgc(data(k,1),data(k,2),:) = [255,255,255];
      end
  else
      col=rand(1,3)*255; 
      for j = 1:size(ii)
          k = ii(j);
        imgc(data(k,1),data(k,2),:) = col;
      end
  end
end
imgc = uint8(imgc);
imgc = imrotate(imgc,90);%旋转
imgc = flipdim(imgc,1);
imgc=imresize(imgc,20,'nearest'); %双线性插值 
imshow(imgc);
title(['类的个数: ',num2str(size(unique(idx)))]);
%imshow(img);