clc;
I=imread('ppp.png');
% I=imnoise(I,'salt & pepper', 0.05);
I=rgb2gray(I);
% figure;
 figure(1);
 subplot(1,2,1),imshow(I);
title('ԭͼ')
M=I;
num=1;
I=double(I);
M=double(M);
flag11=1;
H=256;
L=256;
for i=1:H
    for j=1:L
        flag(i,j)=1;
    end
end
%i������
%j������
for i=1:H%��ѭ��
    for j=1:L%��ѭ��
        omiga=2;
        %%%%%%%%ȷ������
       while flag(i,j)==1
        zuo=i-omiga;
        xia=j-omiga;
        you=i+omiga;
        shang=j+omiga;
    if zuo<1
        zuo=1;
    end
    if xia<1
        xia=1;
    end
    if you>L
        you=L;
    end
    if shang>H
        shang=H;
    end
    %%%%%%%����ȷ������
    %%%%%%%%%%%ȷ�������Сֵ
    smin=I(i,j);
    smax=I(i,j);
    total=(you-zuo+1)*(shang-xia+1);
    vect1=zeros(1,total-1);
    kn=1;
    for in=zuo:you
        for jn=xia:shang
            if ((in==i&jn==j)==0)
        vect1(1,kn)=I(in,jn);
        kn=kn+1;
    end
    end
end
    smin=nanmin(vect1);
    smax=nanmax(vect1);
    smed=nanmedian(vect1);
    %%%%%%%%%%%ȷ�������Сֵ����
    if (smed-smin)>0&(smax-smed)>0
    if smin<M(i,j)&M(i,j)<smax
        flag(i,j)=0;
    else 
        M(i,j)=smed;
        I(i,j)=smed;
        flag(i,j)=0;
    end
    
else
    omiga=omiga+2;
    if omiga>=5
    flag11=0;
    end
    if omiga>=17
  flag11=0;
    M(i,j)=smed;
    flag(i,j)=0;
    end
   end 
    
end%while
    end%��ѭ��
end%��ѭ��
I=uint8(M);
imwrite(I,'lvbo.bmp','bmp');%����ͼ��
figure(1);
subplot(1,2,2),imshow(I);
title('��ֵ�˲���ͼ��')
%����Ϊͼ����ǿ����
M=imread('lvbo.bmp');
%M=rgb2gray(M);
[m,n]=size(M);
GP=zeros(1,256);
for k=0:255
    GP(k+1)=length(find(M==k))/(m*n);
end
figure(2);
subplot(2,2,1),imshow(M);
title('��ֵ�˲���ͼ��');
subplot(2,2,2),bar(0:255,GP,'g')
title('�˲�ͼ��ֱ��ͼ')
xlabel('�Ҷ�ֵ')
ylabel('���ָ���')
S1=zeros(1,256);
for i=1:256
    for j=1:i
        S1(i)=GP(j)+S1(i);
    end
end
S2=round((S1*256)+0.5);
for i=1:256
    GPeq(i)=sum(GP(find(S2==i)));
end
figure(2);
 subplot(2,2,3),bar(0:255,GPeq,'b')
title('���⻯���ֱ��ͼ')
xlabel('�Ҷ�ֵ')
ylabel('���ָ���')
I1=I;
for i=0:255
    I1(find(I==i))=S2(i+1);
end
imwrite(I1,'junzhi.jpg','jpg');%����ͼ��
figure(2);
subplot(2,2,4),imshow(I1);
title('���⻯��ͼ��')
% imwrite(I1,'PicEqual,bmp');
%����Ϊͼ��ָ��
FCM(I1);

