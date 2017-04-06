tic;
clear;
im=imread('hand.bmp');%读入一幅彩色图像
figure,imshow(im);
im=rgb2gray(im);
[height,width]=size(im);
f=freq(im,height,width);%主频率估计
f=round(1/f);
%f=6;
knoise=6;
%边缘检测
[M, m, or, featType, PC, EO] = phasecong2(im,4 ,6 , f, 2.1, 0.55, 1.2, knoise,0.4 ,10 );%相位一致性检测边缘
M=M+m;
nm = nonmaxsup(M, or, 1.5);  %非极大值抑制
[th,tl]=th_tl(nm);%自动高低阈值获取
edgim = hysthresh(nm, th, tl); %高低阈值进行二值化
edgeim = bwmorph(edgim,'skel',Inf);%骨架化
%figure,imshow(M);
%figure,imshow(edgeim);
%将短小边界去掉
LL=bwlabel(edgeim);
ma=max(LL(:));
for n=1:ma
x=find(LL==n);
len=length(x);
if len<30
for k=1:len
edgeim(x(k))=0;
end
end
end
[seedsh,seedsw]=find(edgeim==1);%寻找边缘为1的点并绘制直方图
len=length(seedsh);
for i=1:len
    hist(i)=im(seedsh(i),seedsw(i));
end
x=imhist(im);
num=histogsmooth(x,100);
%num=3;
[mu,mask]=kmeans1(hist,num);
%ed=zeros(height,width);
%for i=1:len
%ed(seedsh(i),seedsw(i))=mask(i);
%end
%figure,imagesc(ed);
me=zeros(1,num);
std=zeros(1,num);
for i=1:num
    Y=find(mask==i);
    lenY=length(Y);
    for j=1:lenY
       shuzu(j)=im(seedsh(Y(j)),seedsw(Y(j)));
    end
    me(i)=mean2(shuzu);
    std(i)=std2(shuzu);
end
im=double(im);
L=zeros(height,width);
seedh=zeros(1,height*width);%多种子点的区域生长
seedw=zeros(1,height*width);
for sp=1:len
seedh(1)=seedsh(sp);%设置种子点
seedw(1)=seedsw(sp);
stackpoint=1;
%求取阈值
T=1.5*std(mask(sp));

m=1;
count=0;
while stackpoint~=0
    count=count+1;
    h=seedh(stackpoint);
    w=seedw(stackpoint);
    stackpoint=stackpoint-1;
    L(h,w)=uint8(255*mask(sp)/num);edgeim(h,w)=1;
    %hist(m)=hm(h,w);
    %m=m+1;
%判断左边的点    
  if w>1
      d=abs(im(h,w-1)-me(mask(sp)));
     if (d<T)&(edgeim(h,w-1)==0)
      stackpoint=stackpoint+1;
      seedh(stackpoint)=h;
      seedw(stackpoint)=w-1;
  end
  d=0;
end
%判断上边的点
  if h>1
      d=abs(im(h-1,w)-me(mask(sp)));
     if (d<T)&(edgeim(h-1,w)==0)
      stackpoint=stackpoint+1;
      seedh(stackpoint)=h-1;
      seedw(stackpoint)=w;
  end
  d=0;
end
%判断右面的点
if w<width
    d=abs(im(h,w+1)-me(mask(sp)));
     if (d<T)&(edgeim(h,w+1)==0)
      stackpoint=stackpoint+1;
      seedh(stackpoint)=h;
      seedw(stackpoint)=w+1;
  end
  d=0;
end
%判断下边的点
if h<height
    d=abs(im(h+1,w)-me(mask(sp)));
     if (d<T)&(edgeim(h+1,w)==0)
      stackpoint=stackpoint+1;
      seedh(stackpoint)=h+1;
      seedw(stackpoint)=w;
  end
  d=0;
end
end
end
%figure,imagesc(L);
%figure,imshow(uint8(L));
se=strel('square',3);
LL=imdilate(L,se);
%figure,imagesc(LL);
%figure,imshow(uint8(LL));
for i=1:num
   [h,w]=find(LL==uint8(255*i/num));
   len=length(h);
   for j=1:len
       shuzu(j)=im(h(j),w(j));
   end
   me(i)=mean2(shuzu);
   clear shuzu;
end
X=find(LL==0);
len=length(X);
d=zeros(1,num);
for i=1:len
    for j=1:num
        d(j)=abs(me(j)-im(X(i)));
    end
    XX=find(d==min(d));
    biaoji=XX(1);
    LL(X(i))=uint8(biaoji*255/num);
    clear d,XX;
end    
LL=medfilt2(LL);
figure,imshow(uint8(LL));
toc;
t=toc;