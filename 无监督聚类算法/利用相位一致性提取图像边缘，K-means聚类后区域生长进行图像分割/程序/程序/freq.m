%clear;
function [f0]=freq(im,height,width)
%im=imread('2.bmp');
%figure,imshow(im);
%im=rgb2gray(im);
%[height,width]=size(im);
m=min(height,width);
n=m;
f=fft2(im,m,n);
h=fftshift(abs(f));
if mod(m,2)==0   %获取频谱图像的中心
    centerh=m/2+1;
    centerw=m/2+1;
else
    centerh=(m+1)/2;
    centerw=(m+1)/2;
end
h(centerh-5:centerh+5,centerw-5:centerw+5)=0;%将零频率对应的点的周围区域置零
%h(76-3:76+3,75-3:75+3)=0;
%h(79-3:79+3,80-3:80+3)=0;
%h(71-3:71+3,72-3:72+3)=0;
[m,n]=size(h);
%将k*k大小的图像分割成4个大小（k-1)*(k-1)的图像
iformal1=1,iformaln=m;
jformal1=1,jformaln=n;

rformal1=sum(h(1,jformal1+1:jformaln-1));
rformaln=sum(h(m,jformal1+1:jformaln-1));
cformal1=sum(h(iformal1+1:iformaln-1,1));
cformaln=sum(h(iformal1+1:iformaln-1,n));
%分别求取四个子图像中所包含点的模值的和
s1=rformal1+cformal1+h(iformal1,jformal1);
s2=rformal1+cformaln+h(iformal1,jformaln);
s3=rformaln+cformal1+h(iformaln,jformal1);
s4=rformaln+cformaln+h(iformaln,jformaln);
s=[s1,s2,s3,s4];
%比较四个和的大小，以和最小的作为下一轮检测的输入图像
smax=find(s==max(s));
flag=1;                  
while flag
switch (smax)
    case 1
        icurrent1=iformal1; icurrentn=iformaln-1;
        jcurrent1=jformal1; jcurrentn=jformaln-1;
        rcurrent1=rformal1-h(icurrent1,jcurrentn);
        rcurrentn=sum(h(icurrentn,jcurrent1+1:jcurrentn-1));
        ccurrent1=cformal1-h(icurrentn,jcurrent1);
        ccurrentn=sum(h(icurrent1+1:icurrentn-1,jcurrentn));
    case 2
         icurrent1=iformal1; icurrentn=iformaln-1;
         jcurrent1=jformal1+1; jcurrentn=jformaln;
         rcurrent1=rformal1-h(icurrent1,jcurrent1);
         rcurrentn=sum(h(icurrentn,jcurrent1+1:jcurrentn-1));
         ccurrent1=sum(h(icurrent1+1:icurrentn-1,jcurrent1));
         ccurrentn=cformaln-h(icurrentn,jcurrentn);
     case 3
         icurrent1=iformal1+1; icurrentn=iformaln;
         jcurrent1=jformal1; jcurrentn=jformaln-1;
         rcurrent1=sum(h(icurrent1,jcurrent1+1:jcurrentn-1));
         rcurrentn=rformaln-h(icurrentn,jcurrentn);
         ccurrent1=cformal1-h(icurrent1,jcurrent1);
         ccurrentn=sum(h(icurrent1+1:icurrentn-1,jcurrentn));
     case 4
         icurrent1=iformal1+1; icurrentn=iformaln;
         jcurrent1=jformal1+1; jcurrentn=jformaln;
         rcurrent1=sum(h(icurrent1,jcurrent1+1:jcurrentn-1));
         rcurrentn=rformaln-h(icurrentn,jcurrent1);
         ccurrent1=sum(h(icurrent1+1:icurrentn-1,jcurrent1));
         ccurrentn=cformaln-h(icurrent1,jcurrentn);
 end
 
  s1=rcurrent1+ccurrent1+h(icurrent1,jcurrent1);
  s2=rcurrent1+ccurrentn+h(icurrent1,jcurrentn);
  s3=rcurrentn+ccurrent1+h(icurrentn,jcurrent1);
  s4=rcurrentn+ccurrentn+h(icurrentn,jcurrentn);
  s=[s1,s2,s3,s4];
  smaxcurrent=find(s==max(s));
  %if size(smaxcurrent)>1
      smax=smaxcurrent(1);
      % end
 a=h(icurrent1:icurrentn,jcurrent1:jcurrentn);
  [x,y]=size(a);
 if x*y==4
      flag=0;
 else
     flag=1;
      iformal1=icurrent1; iformaln=icurrentn;
      jformal1=jcurrent1; jformaln=jcurrentn;
      rformal1=rcurrent1; rformaln=rcurrentn;
      cformal1=ccurrent1; cformaln=ccurrentn;
  end
end
[fx fy]=find(a==max(a(:)));
b=[icurrent1 icurrentn];
c=[jcurrent1 jcurrentn];
f1=b(fx);
f2=c(fy);
f=sqrt((f1-centerh)^2+(f2-centerw)^2);
u=abs(f1-centerh);
v=abs(f2-centerw);
f0=sqrt(u^2+v^2);
if mod(m,2)==0
    f0=f0/(m/2);
else
    f0=f0/((m+1)/2);
end
