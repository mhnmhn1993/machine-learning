clc;clear;clf;
img=imread('depth.png'); %��ȡ���ͼ��
figure(1);imshow(img);title('ԭͼ');
Narrow = 0.25;
img_s=imresize(img,Narrow,'nearest'); %˫���Բ�ֵ 

figure(2); imshow(img_s);title('ѹ�����ͼ��');
[M,N]=size(img_s);            
[XX,YY]=meshgrid(1:N,1:M);
data=[XX(:) YY(:) double(img_s(:))*Narrow*0.5];%�����ͼת��Ϊ3ά���ݼ�
data = double(data);
figure(3);plot3(data(:,1),data(:,2),data(:,3),'.');title('3ά����');
[n,~]=size(data);
%����ŷ�Ͼ�����෴����Ϊ���ƶȾ���
S = -squareform(pdist(data));
[row,col,v]=find(S);%�ҵ���Ϊ0��Ԫ��
%v��һ��n*(n-1)��Ԫ�صķ�0����������֮��ȡ��λ��
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
imgc = imrotate(imgc,90);%��ת
imgc = flipdim(imgc,1);
imgc=imresize(imgc,20,'nearest'); %˫���Բ�ֵ 
imshow(imgc);
title(['��ĸ���: ',num2str(size(unique(idx)))]);
%imshow(img);