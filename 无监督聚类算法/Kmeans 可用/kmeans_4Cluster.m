% close all;
[RGB,map]= imread ('tm2005mask.jpg'); %����
figure,imshow(RGB);title(' ͼһ ԭͼ��')
img=rgb2gray(RGB);
[m,n]=size(img);
figure
subplot(2,2,1),imshow(img);title(' ͼ�� ԭͼ��ĻҶ�ͼ��')
subplot(2,2,2),imhist(img);title(' ͼ�� ����ǰ�ĻҶ�ͼ��ֱ��ͼ') 
img=double(img);
tic
for i=1:m*n
    c1(1)=25;
    c2(1)=125;
    c3(1)=200; 
    c4(4)=245;%ѡ��k����ʼ��������
    
    r=abs(img-c1(i));
    g=abs(img-c2(i)); 
    b=abs(img-c3(i));
    s=abs(img-c4(i));%��������ػҶ���������ĵľ���
  
    n_r=find(r-g<=0&r-b<=0&r-s<=0);%Ѱ�Ҿ�������
    n_g=find(r-g>=0&g-b<=0&g-s<=0);
    n_b=find(g-b>0&r-b>0&b-s<=0);
    n_s=find(s-r<=0&s-g<=0&s-b<=0);
 
    i=i+1;%���¾�������
    
    c1(i)=sum(img(n_r))/length(n_r); %�����еͻҶ����ȡƽ������Ϊ��һ���ͻҶ�����
    c2(i)=sum(img(n_g))/length(n_g); %�����еͻҶ����ȡƽ������Ϊ��һ���м�Ҷ�����
    c3(i)=sum(img(n_b))/length(n_b); %�����еͻҶ����ȡƽ������Ϊ��һ���м�Ҷ�����  
    c4(i)=sum(img(n_s))/length(n_s); %�����еͻҶ����ȡƽ������Ϊ��һ���߻Ҷ�����
    
    d1(i)=abs(c1(i)-c1(i-1));%������������׼��
    d2(i)=abs(c2(i)-c2(i-1));
    d3(i)=abs(c3(i)-c3(i-1));
    d4(i)=abs(c4(i)-c4(i-1));

   
    if (d1(i)==0&&d2(i)==0&&d3(i)==0&&d4(i)==0)%������������׼��
        V1=c1(i);%���յľ�������
        V2=c2(i);
        V3=c3(i);
        V4=c4(i);
        k=i; 
       break;
    end
    
end
toc
disp('�㷨����������')
(k-1)
disp('�������ģ�')
V1
V2
V3
V4

img(n_r)=V1;  
img(n_g)=V2;  
img(n_b)=V3;  
img(n_s)=V4;
img=uint8(img);
subplot(2,2,3),imshow(img);title(' ͼ�� ������ͼ��') 
%subplot(2,2,4),imhist(img);title(' ͼ�� �����ĻҶ�ͼ��ֱ��ͼ')
figure
imshow(img);title(' ͼ�� ������ͼ��') 
