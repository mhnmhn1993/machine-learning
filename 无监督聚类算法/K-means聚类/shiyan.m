% K-mean��������
clear all;
%%%%%%%%%%%%����7��ң��ͼ������%%%%%%%%%%%%%%%%%%
b(1,:,:)=imread('L1.bmp');
b(2,:,:)=imread('L2.bmp');
b(3,:,:)=imread('L3.bmp');
b(4,:,:)=imread('L4.bmp');
b(5,:,:)=imread('L5.bmp');
b(6,:,:)=imread('L6.bmp');
b(7,:,:)=imread('L7.bmp');
b=double(b);          %ת��Ϊ˫��������%
[temp,M,N]=size(b);
N1=[10 20 20 10 10 15 10]';N2=[20 30 20 100 120 60 50]';N3=[80 80 80 200 180 80 110]';    %��ͼ��ѡȡ��ʼ��������%
number1=0;number2=0;number3=0;             %�������������ص������־��%
%�趨������ֹ���������ӽ���������ֹ���㣩%
T=zeros(M,N);                      %TΪ�жϷ����־��%
e_new=10;e_old=20;times=0;         %�ж������̶ȵı���%
format long;
while(abs(e_old-e_new)>0.00001)
    times=times+1
    e_new
    e_old=e_new;
    e_new=0;
    for i=1:M
        for j=1:N
    %%%%%%%%%%%%%%%%%%%%%%%%%�ж��������%%%%%%%%%%%%%%%%%%%%%%%%%%
            d1=sum((b(:,i,j)-N1).^2);                 %������������ʼ��������֮��ľ���%
            d2=sum((b(:,i,j)-N2).^2);
            d3=sum((b(:,i,j)-N3).^2);   
            min_d=min([d1,d2,d3]);              %��K-mean�㷨�ж��������������%
            if (d1==min_d) 
                T(i,j)=1;    
                number1=number1+1; 
            elseif (d2==min_d)
                T(i,j)=2;   
                number2=number2+1;
            elseif (d3==min_d)
                T(i,j)=3;   
                number3=number3+1;
            end       
        e_new=e_new+sqrt(min_d);                       %���������ֹ����%
        end
    end
   N1=N1*0;N2=N2*0;N3=N3*0;                             %����������%
    e_new=e_new/M/N;
    %%%%%%%%%%%%%%%%%%%%%%���������������%%%%%%%%%%%%%%%%%%%%%%%%%
    for i=1:M
        for j=1:N
            if (T(i,j)==1) N1=N1+b(:,i,j);
            elseif (T(i,j)==2) N2=N2+b(:,i,j);
            elseif (T(i,j)==3) N3=N3+b(:,i,j);
            end
        end
    end
    N1=N1/number1; number1=0;    
    N2=N2/number2; number2=0;
    N3=N3/number3; number3=0;    
end

result=zeros(M,N,3); % �����ϲ�ɫ��RGB��ͼ
for  i=1:M
    for j=1:N                 %     R����            G����            B����
        if (T(i,j)==1)     result(i,j,1)=135;result(i,j,2)=206; result(i,j,3)=250;        %LightSkyBlue%
        elseif(T(i,j)==2)  result(i,j,1)=255;result(i,j,2)=236; result(i,j,3)=139;        %LightGoldenrod1%
        elseif(T(i,j)==3)  result(i,j,1)=244;result(i,j,2)=164;result(i,j,3)=96;          %SandyBrown%
        end
    end
end
result=uint8(result);
figure;
imshow(result);                            %��ʾ�����Ժ��ϲ�ɫ�Ľ��ͼ%
imwrite(result,'result.bmp');


    
    
    
    
    
