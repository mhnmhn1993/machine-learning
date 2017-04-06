% close all;
[RGB,map]= imread ('tm2005mask.jpg'); %读入
figure,imshow(RGB);title(' 图一 原图像')
img=rgb2gray(RGB);
[m,n]=size(img);
figure
subplot(2,2,1),imshow(img);title(' 图二 原图像的灰度图像')
subplot(2,2,2),imhist(img);title(' 图三 聚类前的灰度图像直方图') 
img=double(img);
tic
for i=1:m*n
    c1(1)=25;
    c2(1)=125;
    c3(1)=200; 
    c4(4)=245;%选择k个初始聚类中心
    
    r=abs(img-c1(i));
    g=abs(img-c2(i)); 
    b=abs(img-c3(i));
    s=abs(img-c4(i));%计算各像素灰度与聚类中心的距离
  
    n_r=find(r-g<=0&r-b<=0&r-s<=0);%寻找聚类中心
    n_g=find(r-g>=0&g-b<=0&g-s<=0);
    n_b=find(g-b>0&r-b>0&b-s<=0);
    n_s=find(s-r<=0&s-g<=0&s-b<=0);
 
    i=i+1;%更新聚类中心
    
    c1(i)=sum(img(n_r))/length(n_r); %将所有低灰度求和取平均，作为下一个低灰度中心
    c2(i)=sum(img(n_g))/length(n_g); %将所有低灰度求和取平均，作为下一个中间灰度中心
    c3(i)=sum(img(n_b))/length(n_b); %将所有低灰度求和取平均，作为下一个中间灰度中心  
    c4(i)=sum(img(n_s))/length(n_s); %将所有低灰度求和取平均，作为下一个高灰度中心
    
    d1(i)=abs(c1(i)-c1(i-1));%聚类中心收敛准则
    d2(i)=abs(c2(i)-c2(i-1));
    d3(i)=abs(c3(i)-c3(i-1));
    d4(i)=abs(c4(i)-c4(i-1));

   
    if (d1(i)==0&&d2(i)==0&&d3(i)==0&&d4(i)==0)%聚类中心收敛准则
        V1=c1(i);%最终的聚类中心
        V2=c2(i);
        V3=c3(i);
        V4=c4(i);
        k=i; 
       break;
    end
    
end
toc
disp('算法迭代次数：')
(k-1)
disp('聚类中心：')
V1
V2
V3
V4

img(n_r)=V1;  
img(n_g)=V2;  
img(n_b)=V3;  
img(n_s)=V4;
img=uint8(img);
subplot(2,2,3),imshow(img);title(' 图四 聚类后的图像') 
%subplot(2,2,4),imhist(img);title(' 图五 聚类后的灰度图像直方图')
figure
imshow(img);title(' 图四 聚类后的图像') 
