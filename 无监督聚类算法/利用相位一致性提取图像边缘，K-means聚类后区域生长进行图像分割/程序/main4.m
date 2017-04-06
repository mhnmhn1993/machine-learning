tic;
clear;
im=imread('ppp.png');%����һ����ɫͼ��
figure,imshow(im);
im=rgb2hsv(im);
hm=im(:,:,1);figure,imshow(hm);sm=im(:,:,2);figure,imshow(sm);vm=im(:,:,3);figure,imshow(vm);
[height,width]=size(vm);
f=freq(vm,height,width);%��Ƶ�ʹ���
f=round(1/f);
%f=8;
knoise=6;
%��Ե���
[M1, m1, or1, featType, PC, EO] = phasecong2(hm,4 ,6 , f, 2.1, 0.55, 1.2, knoise,0.4 ,10 );%��λһ���Լ���Ե
M1=M1+m1;
nm1 = nonmaxsup(M1, or1, 1.5); 
[th,tl]=th_tl(nm1);
edgim1 = hysthresh(nm1, th, tl); 
edgeim1 = bwmorph(edgim1,'skel',Inf);
[M2, m2, or2, featType, PC, EO] = phasecong2(sm,4 ,6 , f, 2.1, 0.55, 1.2, knoise,0.4 ,10 );
M2=M2+m2;   
nm2= nonmaxsup(M2, or2, 1.5);  %�Ǽ���ֵ����
[th,tl]=th_tl(nm2);%�Զ��ߵ���ֵ��ȡ
edgim2 = hysthresh(nm2, th, tl); %�ߵ���ֵ���ж�ֵ��
edgeim2 = bwmorph(edgim2,'skel',Inf);%�Ǽܻ�
edgeim=edgeim1+edgeim2;
edgeim = bwmorph(edgeim,'spur');
figure,imshow(edgeim);
%����С�߽�ȥ��
LL=bwlabel(edgeim);
ma=max(LL(:));
for n=1:ma
x=find(LL==n);
len=length(x);
if len<50
for k=1:len
edgeim(x(k))=0;
end
end
end
figure,imshow(edgeim);
[seedsh,seedsw]=find(edgeim==1);%Ѱ�ұ�ԵΪ1�ĵ㲢����ֱ��ͼ
len=length(seedsh);
smax=max(max(vm));
smin=min(min(vm));
im=uint8(255/(smax-smin).*vm);
%im=histeq(im);
for i=1:len
    hist(i)=im(seedsh(i),seedsw(i));
end
figure,imhist(hist);
x=imhist(im);
num=histogsmooth(x,50);
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
seedh=zeros(1,height*width);%�����ӵ����������
seedw=zeros(1,height*width);
for sp=1:len
seedh(1)=seedsh(sp);%�������ӵ�
seedw(1)=seedsw(sp);
stackpoint=1;
%��ȡ��ֵ
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
%�ж���ߵĵ�    
  if w>1
      d=abs(im(h,w-1)-me(mask(sp)));
     if (d<T)&(edgeim(h,w-1)==0)
      stackpoint=stackpoint+1;
      seedh(stackpoint)=h;
      seedw(stackpoint)=w-1;
  end
  d=0;
end
%�ж��ϱߵĵ�
  if h>1
      d=abs(im(h-1,w)-me(mask(sp)));
     if (d<T)&(edgeim(h-1,w)==0)
      stackpoint=stackpoint+1;
      seedh(stackpoint)=h-1;
      seedw(stackpoint)=w;
  end
  d=0;
end
%�ж�����ĵ�
if w<width
    d=abs(im(h,w+1)-me(mask(sp)));
     if (d<T)&(edgeim(h,w+1)==0)
      stackpoint=stackpoint+1;
      seedh(stackpoint)=h;
      seedw(stackpoint)=w+1;
  end
  d=0;
end
%�ж��±ߵĵ�
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

se=strel('square',1);
LL=imopen(L,se);
%LL=imdilate(L,se);
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
LL=uint8(medfilt2(LL));
%LL=uint8(medfilt2(LL));
figure,imshow(LL);
%se=strel('disk',4);
%LL=imopen(LL,se);
se=strel('disk',6);
LL=imclose(LL,se);
figure,imshow(LL);
toc;
t=toc;