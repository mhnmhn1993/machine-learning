clear all;

%<----------------------------!基于小波变换和K-Mean算法的遥感图像分类-------------------------
%读入样本11,即遥感图像的背景
I=imread('11.jpg');
I=rgb2gray(I);
%灰度值归一化
I=im2double(I);
%二维小波分解
[cA,cH,cV,cD]=dwt2(I,'db1');
[M N]=size(I);
M=M/2;N=N/2;
A11=0;H11=0;V11=0;D11=0;
for i=1:M
    for j=1:N
        A11=A11+cA(i,j)/(M*N);
        H11=H11+cH(i,j)/(M*N);
        V11=V11+cV(i,j)/(M*N);
        D11=D11+cD(i,j)/(M*N);
    end
end
%得出特征向量T11
T11=[A11;H11;V11;D11]';

%读入样本图像1
I=imread('1.jpg');
I=rgb2gray(I);
%灰度值归一化
I=im2double(I);
%二维小波分解
[cA,cH,cV,cD]=dwt2(I,'db1');
[M N]=size(I);
M=M/2;N=N/2;
A1=0;H1=0;V1=0;D1=0;
for i=1:M
    for j=1:N
        A1=A1+cA(i,j)/(M*N);
        H1=H1+cH(i,j)/(M*N);
        V1=V1+cV(i,j)/(M*N);
        D1=D1+cD(i,j)/(M*N);
    end
end
T1=[A1;H1;V1;D1]';

%读入样本图像2
I=imread('2.jpg');
I=rgb2gray(I);
%灰度值归一化
I=im2double(I);
%二维小波分解
[cA,cH,cV,cD]=dwt2(I,'db1');
[M N]=size(I);
M=M/2;N=N/2;
A2=0;H2=0;V2=0;D2=0;
for i=1:M
    for j=1:N
        A2=A2+cA(i,j)/(M*N);
        H2=H2+cH(i,j)/(M*N);
        V2=V2+cV(i,j)/(M*N);
        D2=D2+cD(i,j)/(M*N);
    end
end
T2=[A2;H2;V2;D2]';

%读入样本图像3
I=imread('3.jpg');
I=rgb2gray(I);
%灰度值归一化
I=im2double(I);
%二维小波分解
[cA,cH,cV,cD]=dwt2(I,'db1');
[M N]=size(I);
M=M/2;N=N/2;
A3=0;H3=0;V3=0;D3=0;
for i=1:M
    for j=1:N
        A3=A3+cA(i,j)/(M*N);
        H3=H3+cH(i,j)/(M*N);
        V3=V3+cV(i,j)/(M*N);
        D3=D3+cD(i,j)/(M*N);
    end
end
T3=[A3;H3;V3;D3]';

%读入样本图像4
I=imread('4.jpg');
I=rgb2gray(I);
%灰度值归一化
I=im2double(I);
%二维小波分解
[cA,cH,cV,cD]=dwt2(I,'db1');
[M N]=size(I);
M=M/2;N=N/2;
A4=0;H4=0;V4=0;D4=0;
for i=1:M
    for j=1:N
        A4=A4+cA(i,j)/(M*N);
        H4=H4+cH(i,j)/(M*N);
        V4=V4+cV(i,j)/(M*N);
        D4=D4+cD(i,j)/(M*N);
    end
end
T4=[A4;H4;V4;D4]';

%读入样本图像5
I=imread('5.jpg');
I=rgb2gray(I);
%灰度值归一化
I=im2double(I);
%二维小波分解
[cA,cH,cV,cD]=dwt2(I,'db1');
[M N]=size(I);
M=M/2;N=N/2;
A5=0;H5=0;V5=0;D5=0;
for i=1:M
    for j=1:N
        A5=A5+cA(i,j)/(M*N);
        H5=H5+cH(i,j)/(M*N);
        V5=V5+cV(i,j)/(M*N);
        D5=D5+cD(i,j)/(M*N);
    end
end
T5=[A5;H5;V5;D5]';

%读入样本图像6
I=imread('6.jpg');
I=rgb2gray(I);
%灰度值归一化
I=im2double(I);
%二维小波分解
[cA,cH,cV,cD]=dwt2(I,'db1');
[M N]=size(I);
M=M/2;N=N/2;
A6=0;H6=0;V6=0;D6=0;
for i=1:M
    for j=1:N
        A6=A6+cA(i,j)/(M*N);
        H6=H6+cH(i,j)/(M*N);
        V6=V6+cV(i,j)/(M*N);
        D6=D6+cD(i,j)/(M*N);
    end
end
T6=[A6;H6;V6;D6]';

%读入样本图像7
I=imread('7.jpg');
I=rgb2gray(I);
%灰度值归一化
I=im2double(I);
%二维小波分解
[cA,cH,cV,cD]=dwt2(I,'db1');
[M N]=size(I);
M=M/2;N=N/2;
A7=0;H7=0;V7=0;D7=0;
for i=1:M
    for j=1:N
        A7=A7+cA(i,j)/(M*N);
        H7=H7+cH(i,j)/(M*N);
        V7=V7+cV(i,j)/(M*N);
        D7=D7+cD(i,j)/(M*N);
    end
end
T7=[A7;H7;V7;D7]';

%读入样本图像8
I=imread('8.jpg');
I=rgb2gray(I);
%灰度值归一化
I=im2double(I);
%二维小波分解
[cA,cH,cV,cD]=dwt2(I,'db1');
[M N]=size(I);
M=M/2;N=N/2;
A8=0;H8=0;V8=0;D8=0;
for i=1:M
    for j=1:N
        A8=A8+cA(i,j)/(M*N);
        H8=H8+cH(i,j)/(M*N);
        V8=V8+cV(i,j)/(M*N);
        D8=D8+cD(i,j)/(M*N);
    end
end
T8=[A8;H8;V8;D8]';


%读入样本图像9
I=imread('9.jpg');
I=rgb2gray(I);
%灰度值归一化
I=im2double(I);
%二维小波分解
[cA,cH,cV,cD]=dwt2(I,'db1');
[M N]=size(I);
M=M/2;N=N/2;
A9=0;H9=0;V9=0;D9=0;
for i=1:M
    for j=1:N
        A9=A9+cA(i,j)/(M*N);
        H9=H9+cH(i,j)/(M*N);
        V9=V9+cV(i,j)/(M*N);
        D9=D9+cD(i,j)/(M*N);
    end
end
T9=[A9;H9;V9;D9]';

%读入样本图像10
I=imread('10.jpg');
I=rgb2gray(I);
%灰度值归一化
I=im2double(I);
%二维小波分解
[cA,cH,cV,cD]=dwt2(I,'db1');
[M N]=size(I);
M=M/2;N=N/2;
A10=0;H10=0;V10=0;D10=0;
for i=1:M
    for j=1:N
        A10=A10+cA(i,j)/(M*N);
        H10=H10+cH(i,j)/(M*N);
        V10=V10+cV(i,j)/(M*N);
        D10=D10+cD(i,j)/(M*N);
    end
end
T10=[A10;H10;V10;D10]';

%读入待分类遥感图像
I=imread('tm2000mask.jpg');
I=rgb2gray(I);
%灰度值归一化
I=im2double(I);
%二维小波分解
[cA,cH,cV,cD]=dwt2(I,'db1');
[M N]=size(I);
M=M/2;N=N/2;
P=[cA;cH;cV;cD];

%K-Mean均值分类
K=zeros(M,N);                      %K为判断分类标志符%
e_new=10;e_old=20;times=0;         %判断收敛程度的变量%
number1=0;number2=0;number3=0;number4=0;number5=0;number6=0;number7=0;number8=0;number9=0;number10=0;number11=0;
format long;
while(abs(e_old-e_new)>0.00001)
    times=times+1;
    e_new;
    e_old=e_new;
    e_new=0;
    for i=1:M
        for j=1:N
     %%%%%%%%%%%%%%%%%%%%%%%%%判定各点归类%%%%%%%%%%%%%%%%%%%%%%%%%%
            d(1)=((cA(i,j)-T1(1))^2+(cH(i,j)-T1(2))^2+(cV(i,j)-T1(3))^2+(cD(i,j)-T1(4))^2);   %计算各点与各初始分类中心之间的距离
            d(2)=((cA(i,j)-T2(1))^2+(cH(i,j)-T2(2))^2+(cV(i,j)-T2(3))^2+(cD(i,j)-T2(4))^2);
            d(3)=((cA(i,j)-T3(1))^2+(cH(i,j)-T3(2))^2+(cV(i,j)-T3(3))^2+(cD(i,j)-T3(4))^2);
            d(4)=((cA(i,j)-T4(1))^2+(cH(i,j)-T4(2))^2+(cV(i,j)-T4(3))^2+(cD(i,j)-T4(4))^2);
            d(5)=((cA(i,j)-T5(1))^2+(cH(i,j)-T5(2))^2+(cV(i,j)-T5(3))^2+(cD(i,j)-T5(4))^2);
            d(6)=((cA(i,j)-T6(1))^2+(cH(i,j)-T6(2))^2+(cV(i,j)-T6(3))^2+(cD(i,j)-T6(4))^2);
            d(7)=((cA(i,j)-T7(1))^2+(cH(i,j)-T7(2))^2+(cV(i,j)-T7(3))^2+(cD(i,j)-T7(4))^2);
            d(8)=((cA(i,j)-T8(1))^2+(cH(i,j)-T8(2))^2+(cV(i,j)-T8(3))^2+(cD(i,j)-T8(4))^2);
            d(9)=((cA(i,j)-T9(1))^2+(cH(i,j)-T9(2))^2+(cV(i,j)-T9(3))^2+(cD(i,j)-T9(4))^2);
            d(10)=((cA(i,j)-T10(1))^2+(cH(i,j)-T10(2))^2+(cV(i,j)-T10(3))^2+(cD(i,j)-T10(4))^2);
            d(11)=((cA(i,j)-T11(1))^2+(cH(i,j)-T11(2))^2+(cV(i,j)-T11(3))^2+(cD(i,j)-T11(4))^2);  
            min_d=min(d);              %由K-mean算法判定各点所属的类别
            if (d(1)==min_d) 
                K(i,j)=1;    
                number1=number1+1; 
            elseif (d(2)==min_d)
                K(i,j)=2;   
                number2=number2+1;
            elseif (d(3)==min_d)
                K(i,j)=3;   
                number3=number3+1;
            elseif (d(4)==min_d)
                K(i,j)=4;
                number4=number4+1;
            elseif (d(5)==min_d)
                K(i,j)=5;
                number5=number5+1;
            elseif (d(6)==min_d)
                K(i,j)=6;
                number6=number6+1;
            elseif (d(7)==min_d)
                K(i,j)=7;
                number7=number7+1;
            elseif (d(8)==min_d)
                K(i,j)=8;
                number8=number8+1;
            elseif (d(9)==min_d)
                K(i,j)=9;
                number9=number9+1;
            elseif (d(10)==min_d)
                K(i,j)=10;
                number10=number10+1;
            elseif (d(11)==min_d)
                K(i,j)=11;
                number11=number11+1;
            end       
        e_new=e_new+sqrt(min_d);                       %计算迭代终止条件
        end
    end
   T1=T1*0;T2=T2*0;T3=T3*0;T4=T4*0;T5=T5*0;T6=T6*0;T7=T7*0;T8=T8*0;T9=T9*0;T10=T10*0;T11=T11*0;    %调整聚类中心
    e_new=e_new/M/N;
    %%%%%%%%%%%%%%%%%%%%%%调整各聚类的中心%%%%%%%%%%%%%%%%%%%%%%%%%
    for i=1:M
        for j=1:N
            if (K(i,j)==1) T1(1)=T1(1)+cA(i,j);T1(2)=T1(2)+cH(i,j);T1(3)=T1(3)+cV(i,j);T1(4)=T1(4)+cD(i,j);
            elseif (K(i,j)==2) T2(1)=T2(1)+cA(i,j);T2(2)=T2(2)+cH(i,j);T2(3)=T2(3)+cV(i,j);T2(4)=T2(4)+cD(i,j);
            elseif (K(i,j)==3) T3(1)=T3(1)+cA(i,j);T3(2)=T3(2)+cH(i,j);T3(3)=T3(3)+cV(i,j);T3(4)=T3(4)+cD(i,j);
            elseif (K(i,j)==4) T4(1)=T4(1)+cA(i,j);T4(2)=T4(2)+cH(i,j);T4(3)=T4(3)+cV(i,j);T4(4)=T4(4)+cD(i,j);
            elseif (K(i,j)==5) T5(1)=T5(1)+cA(i,j);T5(2)=T5(2)+cH(i,j);T5(3)=T5(3)+cV(i,j);T5(4)=T5(4)+cD(i,j);
            elseif (K(i,j)==6) T6(1)=T6(1)+cA(i,j);T6(2)=T6(2)+cH(i,j);T6(3)=T6(3)+cV(i,j);T6(4)=T6(4)+cD(i,j);
            elseif (K(i,j)==7) T7(1)=T7(1)+cA(i,j);T7(2)=T7(2)+cH(i,j);T7(3)=T7(3)+cV(i,j);T7(4)=T7(4)+cD(i,j);
            elseif (K(i,j)==8) T8(1)=T8(1)+cA(i,j);T8(2)=T8(2)+cH(i,j);T8(3)=T8(3)+cV(i,j);T8(4)=T8(4)+cD(i,j);
            elseif (K(i,j)==9) T9(1)=T9(1)+cA(i,j);T9(2)=T9(2)+cH(i,j);T9(3)=T9(3)+cV(i,j);T9(4)=T9(4)+cD(i,j);
            elseif (K(i,j)==10) T10(1)=T10(1)+cA(i,j);T10(2)=T10(2)+cH(i,j);T10(3)=T10(3)+cV(i,j);T10(4)=T10(4)+cD(i,j);
            elseif (K(i,j)==11) T11(1)=T11(1)+cA(i,j);T11(2)=T11(2)+cH(i,j);T11(3)=T11(3)+cV(i,j);T11(4)=T11(4)+cD(i,j);
            end
        end
    end
    
    T1=T1/number1; number1=0;    
    T2=T2/number2; number2=0;
    T3=T3/number3; number3=0;
    T4=T4/number4; number4=0;
    T5=T5/number5; number5=0;
    T6=T6/number6; number6=0;
    T7=T7/number7; number7=0;
    T8=T8/number8; number8=0;
    T9=T9/number9; number9=0;
    T10=T10/number10; number10=0;
    T11=T11/number11; number11=0;
    end
    
    
a1=0;a2=0;a3=0;a4=0;a5=0;a6=0;a7=0;a8=0;a9=0;a10=0;a11=0;
result=zeros(M,N,3); % 进行上彩色，RGB彩图
for  i=1:M
    for j=1:N                 %     R分量              G分量                B分量
        if (K(i,j)==1)     result(i,j,1)=0;      result(i,j,2)=0;       result(i,j,3)=0;    a1=a1+1;    
        elseif(K(i,j)==2)  result(i,j,1)=25;     result(i,j,2)=76;      result(i,j,3)=127;  a2=a2+1;     
        elseif(K(i,j)==3)  result(i,j,1)=51;     result(i,j,2)=102;     result(i,j,3)=153;  a3=a3+1;
        elseif(K(i,j)==4)  result(i,j,1)=76;     result(i,j,2)=51;      result(i,j,3)=178;  a4=a4+1;
        elseif(K(i,j)==5)  result(i,j,1)=102;    result(i,j,2)=102;     result(i,j,3)=178;  a5=a5+1;
        elseif(K(i,j)==6)  result(i,j,1)=127;    result(i,j,2)=127;     result(i,j,3)=102;  a6=a6+1;
        elseif(K(i,j)==7)  result(i,j,1)=153;    result(i,j,2)=229;     result(i,j,3)=51;   a7=a7+1;
        elseif(K(i,j)==8)  result(i,j,1)=178;    result(i,j,2)=51;      result(i,j,3)=204;  a8=a8+1;
        elseif(K(i,j)==9)  result(i,j,1)=204;    result(i,j,2)=76;      result(i,j,3)=127;  a9=a9+1;
        elseif(K(i,j)==10) result(i,j,1)=229;    result(i,j,2)=51;      result(i,j,3)=25;   a10=a10+1;
        elseif(K(i,j)==11) result(i,j,1)=255;    result(i,j,2)=127;     result(i,j,3)=153;  a11=a11+1;
        end
    end
end
result=uint8(result);
figure;
imshow(result); 
 title('基于小波变换的图像分类效果');
 
 %绘制地物面积饼状图
  a=[a1,a2,a3,a4,a5,a6,a7,a8,a9,a10];
  figure;
  pie(a);
  legend('海水','农地','绿林地','房屋','养殖场','芦苇','互花米草','海三棱藨草','光滩','未利用地',-1);            %标注图例                                 