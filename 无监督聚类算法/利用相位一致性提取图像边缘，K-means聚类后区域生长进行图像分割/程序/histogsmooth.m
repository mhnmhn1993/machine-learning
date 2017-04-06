%clear;
%i=imread('xibao.bmp');
%figure,imshow(i);
%i=rgb2gray(i);
%x=imhist(i);
%figure,plot(x);
function [num]=histogsmooth(x,n);
%n=50;
for j=1:n
   k=x(1);
   for i=2:255
       x(i-1)=(x(i-1)+x(i)+x(i+1))/3;
   end
   for i=255:-1:2
       x(i)=x(i-1);
   end
   x(1)=k;
end
figure,plot(x);
num=0;
for i=3:253
    if  (x(i)>x(i-1))&x(i-1)>x(i-2)&(x(i)>x(i+1)&x(i+1)>x(i+2))
        num=num+1;
    end
end
%num