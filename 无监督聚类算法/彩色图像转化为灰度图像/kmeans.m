%一种基于K均值的自然图像聚类方法
clc
clear 
pic = imread('ppp.png');  % RGB  三维矩阵
figure;
imshow(pic);
k = 10;  
pic = double(pic);       %double把任何类型数据转换成双精度数值
copy = pic;         
ima = pic(:);            % 变成一列
s = length(ima);
data1 = pic(:,:,1);
data2 = pic(:,:,2);
data3 = pic(:,:,3);

%xlswrite('D:\write2Excel.xls',data1,'sheet1');%数据存储至D盘根目录下
%xlswrite('D:\write2Excel.xls',data2,'sheet2')
%xlswrite('D:\write2Excel.xls',data3,'sheet3')

%建立灰度直方图
%灰度直方图是灰度级的函数，它表示图像中具有某种灰度级的像素的个数，反映了图像中某种灰度出现的频率。
m = max(ima);            % m = max(ima)+1;
h = zeros(1,m);          % 一共是1-m个灰度
hc = zeros(1,m);                                                                                                                           
for i = 1:s
  if(ima(i)>0) 
      h(ima(i )) = h(ima(i)) + 1;
  end
end
ind = find(h);           % 找到所有灰度。h矩阵是每个灰度的个数
hl = length(ind);        % 灰度长度

mu = (1:k)*m/(k);           %初始化聚类中心
while 1
  oldmu = mu;    
  for i = 1:hl      
      c = abs(ind(i)-mu);   %  距离
      cc = find(c==min(c)); %　最短的距离所在矩阵位置即是所属簇类
      hc(ind(i)) = cc(1);   %  hc矩阵是某一灰度所在的簇类 
  end  
  for i = 1:k               %  更新聚类中心
      a = find(hc==i);      %  在hc矩阵中找到属于第i类簇的序列     
      mu(i) = sum(a.*h(a))/sum(h(a));   %灰度定义   
  end  
  if(mu == oldmu)           %  聚类中心不再发生改变
      break;
  end
end

S = size(copy);
mask = zeros(S);
for i = 1:S(1)       %重新绘制图像
    for j = 1:S(2)     
        copy(i,j,1) - mu;
        c = abs(copy(i,j,1)-mu);
        a = find(c==min(c));  
        mask(i,j,1) = mu(a);
    end
end
mask  = uint8(mask); %8位的无符号整形数据,取值范围从0到255
figure;
imshow(mask(:,:,1));
figure;
imshow(mask(:,:,2));
figure;
imshow(mask(:,:,3));




