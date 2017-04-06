% K-mean方法聚类
clear all;
%%%%%%%%%%%%读入7幅遥感图像数据%%%%%%%%%%%%%%%%%%
b(1,:,:)=imread('L1.bmp');
b(2,:,:)=imread('L2.bmp');
b(3,:,:)=imread('L3.bmp');
b(4,:,:)=imread('L4.bmp');
b(5,:,:)=imread('L5.bmp');
b(6,:,:)=imread('L6.bmp');
b(7,:,:)=imread('L7.bmp');
b=double(b);          %转换为双精度类型%
[temp,M,N]=size(b);
N1=[10 20 20 10 10 15 10]';N2=[20 30 20 100 120 60 50]';N3=[80 80 80 200 180 80 110]';    %由图像选取初始分类中心%
number1=0;number2=0;number3=0;             %计量各类中象素点个数标志符%
%设定迭代终止条件（若接近收敛则终止计算）%
T=zeros(M,N);                      %T为判断分类标志符%
e_new=10;e_old=20;times=0;         %判断收敛程度的变量%
format long;
while(abs(e_old-e_new)>0.00001)
    times=times+1
    e_new
    e_old=e_new;
    e_new=0;
    for i=1:M
        for j=1:N
    %%%%%%%%%%%%%%%%%%%%%%%%%判定各点归类%%%%%%%%%%%%%%%%%%%%%%%%%%
            d1=sum((b(:,i,j)-N1).^2);                 %计算各点与各初始分类中心之间的距离%
            d2=sum((b(:,i,j)-N2).^2);
            d3=sum((b(:,i,j)-N3).^2);   
            min_d=min([d1,d2,d3]);              %由K-mean算法判定各点所属的类别%
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
        e_new=e_new+sqrt(min_d);                       %计算迭代终止条件%
        end
    end
   N1=N1*0;N2=N2*0;N3=N3*0;                             %调聚类中心%
    e_new=e_new/M/N;
    %%%%%%%%%%%%%%%%%%%%%%调整各聚类的中心%%%%%%%%%%%%%%%%%%%%%%%%%
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

result=zeros(M,N,3); % 进行上彩色，RGB彩图
for  i=1:M
    for j=1:N                 %     R分量            G分量            B分量
        if (T(i,j)==1)     result(i,j,1)=135;result(i,j,2)=206; result(i,j,3)=250;        %LightSkyBlue%
        elseif(T(i,j)==2)  result(i,j,1)=255;result(i,j,2)=236; result(i,j,3)=139;        %LightGoldenrod1%
        elseif(T(i,j)==3)  result(i,j,1)=244;result(i,j,2)=164;result(i,j,3)=96;          %SandyBrown%
        end
    end
end
result=uint8(result);
figure;
imshow(result);                            %显示处理以后上彩色的结果图%
imwrite(result,'result.bmp');


    
    
    
    
    
