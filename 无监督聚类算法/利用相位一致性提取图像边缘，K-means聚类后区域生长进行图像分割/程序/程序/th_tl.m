%im=imread('ziran.bmp');
%figure,imshow(im);
%im=rgb2gray(im);
function [th,tl]=th_tl(im)
mi=min(im(:));
ma=max(im(:));
[height,width]=size(im);
GDT0=abs(double(mi)-double(ma))/4;
GD=zeros(1,9);
ncount=0;
TotalGDH=zeros(1,height*width);
TotalGDL=zeros(1,height*width);
for h=2:height-1
 for w=2:width-1
   k=1;
   for i=h-1:h+1
     for j=w-1:w+1
         GD(k)=abs(double(im(h,w))-double(im(i,j)));
         k=k+1;
     end
    end
 maxGD=max(GD);
 AGD=mean(GD);
 if maxGD>GDT0
     ncount=ncount+1;
     TotalGDH(ncount)=maxGD;
     TotalGDL(ncount)=AGD;
 end
maxGD=0;
AGD=0;
end
end
th=sum(TotalGDH)/ncount;
tl=sum(TotalGDL)/ncount;

