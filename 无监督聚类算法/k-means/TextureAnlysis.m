clc
clear all;

Image=imread('3B.bmp');
imsize=size(Image,1);
winsize=4;
winnum=imsize/winsize;

ENTmean=zeros(winnum,winnum);    % 熵的均值
% ENTvar=zeros(winnum,winnum);   % 熵的方差，无明显效果
% CONmean=zeros(winnum,winnum);  % 对比度度量的均值，无明显效果
% K1=2;K2=2;
stepsize=4;

for i=1:winsize:imsize
    for j=1:winsize:imsize
         Qmatrix=zeros(winsize,winsize);
         Pmatrix=zeros(stepsize,stepsize,4);
         
         block=Image(i:i+winsize-1,j:j+winsize-1);  % 图像子块
         
         maxpxl=max(max(block));
         minpxl=min(min(block));
         step=(maxpxl-minpxl)/stepsize;
         
         for k=1:stepsize  % 量化
             index=find(block>=minpxl+(k-1)*step & block<minpxl+k*step);
             Qmatrix(index)=k;
         end
         index=find(block>=minpxl+(k-1)*step & block<=maxpxl);
         Qmatrix(index)=k;  
         
         for m=1:stepsize   % 四个方向的共发矩阵（0，45，90，135度）
             for n=1:stepsize
                 
                 if n==1
                     Pmatrix(Qmatrix(m,n),Qmatrix(m,n+1),1)=Pmatrix(Qmatrix(m,n),Qmatrix(m,n+1),1)+1;                     
                     if m~=1
                         Pmatrix(Qmatrix(m,n),Qmatrix(m-1,n+1),2)=Pmatrix(Qmatrix(m,n),Qmatrix(m-1,n+1),2)+1;
                     elseif m~=stepsize
                         Pmatrix(Qmatrix(m,n),Qmatrix(m+1,n+1),4)=Pmatrix(Qmatrix(m,n),Qmatrix(m+1,n+1),4)+1;
                     end                     
                 elseif n==stepsize
                     Pmatrix(Qmatrix(m,n),Qmatrix(m,n-1),1)=Pmatrix(Qmatrix(m,n),Qmatrix(m,n-1),1)+1;                     
                     if m~=1
                         Pmatrix(Qmatrix(m,n),Qmatrix(m-1,n-1),4)=Pmatrix(Qmatrix(m,n),Qmatrix(m-1,n-1),4)+1;
                     elseif m~=stepsize
                         Pmatrix(Qmatrix(m,n),Qmatrix(m+1,n-1),2)=Pmatrix(Qmatrix(m,n),Qmatrix(m+1,n-1),2)+1;
                     end                     
                 else
                     Pmatrix(Qmatrix(m,n),Qmatrix(m,n+1),1)=Pmatrix(Qmatrix(m,n),Qmatrix(m,n+1),1)+1;
                     Pmatrix(Qmatrix(m,n),Qmatrix(m,n-1),1)=Pmatrix(Qmatrix(m,n),Qmatrix(m,n-1),1)+1;                     
                     if m==1
                         Pmatrix(Qmatrix(m,n),Qmatrix(m+1,n-1),2)=Pmatrix(Qmatrix(m,n),Qmatrix(m+1,n-1),2)+1;
                         Pmatrix(Qmatrix(m,n),Qmatrix(m+1,n+1),4)=Pmatrix(Qmatrix(m,n),Qmatrix(m+1,n+1),4)+1;
                     elseif m==stepsize
                         Pmatrix(Qmatrix(m,n),Qmatrix(m-1,n+1),2)=Pmatrix(Qmatrix(m,n),Qmatrix(m-1,n+1),2)+1;
                         Pmatrix(Qmatrix(m,n),Qmatrix(m-1,n-1),4)=Pmatrix(Qmatrix(m,n),Qmatrix(m-1,n-1),4)+1;
                     else
                         Pmatrix(Qmatrix(m,n),Qmatrix(m-1,n+1),2)=Pmatrix(Qmatrix(m,n),Qmatrix(m-1,n+1),2)+1;
                         Pmatrix(Qmatrix(m,n),Qmatrix(m+1,n-1),2)=Pmatrix(Qmatrix(m,n),Qmatrix(m+1,n-1),2)+1;
                         Pmatrix(Qmatrix(m,n),Qmatrix(m+1,n+1),4)=Pmatrix(Qmatrix(m,n),Qmatrix(m+1,n+1),4)+1;
                         Pmatrix(Qmatrix(m,n),Qmatrix(m-1,n-1),4)=Pmatrix(Qmatrix(m,n),Qmatrix(m-1,n-1),4)+1;
                     end                     
                 end  % end if
                 
                 if m==1
                     Pmatrix(Qmatrix(m,n),Qmatrix(m+1,n),3)=Pmatrix(Qmatrix(m,n),Qmatrix(m+1,n),3)+1;
                 elseif m==stepsize
                     Pmatrix(Qmatrix(m,n),Qmatrix(m-1,n),3)=Pmatrix(Qmatrix(m,n),Qmatrix(m-1,n),3)+1;
                 else
                     Pmatrix(Qmatrix(m,n),Qmatrix(m+1,n),3)=Pmatrix(Qmatrix(m,n),Qmatrix(m+1,n),3)+1;
                     Pmatrix(Qmatrix(m,n),Qmatrix(m-1,n),3)=Pmatrix(Qmatrix(m,n),Qmatrix(m-1,n),3)+1;
                 end  % end if               
             end
         end  % end for
        
%          CON=zeros(1,4);
        for k=1:4    
            P=Pmatrix(:,:,k);
            pos=find(P>0);
            pxl=P(pos);
            M(k)=(-1)*sum(sum(pxl.*log10(pxl))); % 四个灰度共发矩阵的熵
                        
%             for m=1:stepsize           % 四个灰度共发矩阵的对比度度量
%                 for n=1:stepsize
%                     pxlabs=abs(m-n);
%                     CON(1,k)=CON(1,k)+(pxlabs^K1)*(Pmatrix(m,n,k)^K2);
%                 end
%             end
        end
        
        ENTmean(floor(i/stepsize)+1,floor(j/stepsize)+1)=mean(M);               % 熵的均值        
%         ENTvar(floor(i/stepsize)+1,floor(j/stepsize)+1)=mean((M-mean(M)).^2); % 熵的方差 
%         CONmean(floor(i/stepsize)+1,floor(j/stepsize)+1)=mean(CON);    % 对比度度量的均值           
    end
end  % end for

Texture=zeros(2,imsize,imsize);  % 特征矩阵

for i=1:winnum
    for j=1:winnum
        Texture(1,[(i-1)*winsize+1:i*winsize],[(i-1)*winsize+1:i*winsize])=ENTmean(i,j);
%         Texture(1,[(i-1)*winsize+1:i*winsize],[(i-1)*winsize+1:i*winsize])=ENTvar(i,j);  % 无明显效果
%         Texture(2,[(i-1)*winsize+1:i*winsize],[(i-1)*winsize+1:i*winsize])=CONmean(i,j); % 无明显效果
    end
end

Imext=wextend('2D','zpd',Image,[4,4]); % zpd 补零  sp0  平滑
Direction=[1,0,0,0,1,0,0,0,1; ...   % 9×9的矩阵，八个方向为1
    0,1,0,0,1,0,0,1,0; ...
    0,0,1,0,1,0,1,0,0; ...
    0,0,0,1,1,1,0,0,0; ...
    1,1,1,1,0,1,1,1,1; ...  % 中心点取0 1无明显差别
    0,0,0,1,1,1,0,0,0; ...
    0,0,1,0,1,0,1,0,0; ...
    0,1,0,0,1,0,0,1,0; ...
    1,0,0,0,1,0,0,0,1];
Direction=double(Direction);
for i=1:imsize
    for j=1:imsize
%         pxlmean=mean(mean(Imext([i:i+2],[j:j+2])));
%         Imsmooth(i,j)=sum(sum((Imext([i:i+2],[j:j+2])-pxlmean).^2)); % 单位面积灰度变化总量，效果不好
%         Imsmooth(i,j)=mean(mean(Imext([i:i+4],[j:j+4])));       % 平均灰度,效果好于单位面积灰度变化总量
        Direction2=double(Imext([i:i+8],[j:j+8]));
        Imsmooth(i,j)=mean(mean(Direction2.*Direction));  % 在邻域内八个方向取均值（应该有比均值更好的）
    end
end
Texture(2,:,:)=Imsmooth;           % 亦可用均值平方，增加差距，但效果各有所长

C1=Texture(:,30,30);               % 初始聚类中心
C2=Texture(:,150,200);
C3=Texture(:,250,400);

T=zeros(imsize,imsize);            % 分类标识矩阵
number1=0; number2=0; number3=0;   % 累计各类中象素个数
e_new=10;e_old=20;

while abs(e_new-e_old)>1
    e_old=e_new;
    e_new=0;
   
    for m=1:imsize
        for n=1:imsize
            d1=sum((Texture(:,m,n)-C1).^2);
            d2=sum((Texture(:,m,n)-C2).^2);  
            d3=sum((Texture(:,m,n)-C3).^2);
            dmin=min([d1,d2,d3]); % 最小距离 的平方
            
            if d1==dmin           % 按照最小距离聚类
                T(m,n)=1;
                number1=number1+1;
            elseif d2==dmin
                T(m,n)=2;
                number2=number2+1;
            elseif d3==dmin
                T(m,n)=3;
                number3=number3+1;
            end
            
            e_new=e_new+sqrt(dmin);            
        end
    end % end for
    
    C1=C1.*0; C2=C2.*0; C3=C3.*0; % 中心置零，稍后更新
    e_new=e_new/(imsize^2);       % 更新误差
    
    for m=1:imsize                % 更新聚类中心
        for n=1:imsize
            if T(m,n)==1
                C1=C1+Texture(:,m,n);
            elseif T(m,n)==2
                C2=C2+Texture(:,m,n);
            elseif T(m,n)==3
                C3=C3+Texture(:,m,n);
            end
        end
    end
    C1=C1./number1; number1=0;    % 新的聚类中心
    C2=C2./number2; number2=0;
    C3=C3./number3; number3=0;       
end % end while

RGBresult=zeros(imsize,imsize,3); % 进行上彩色，RGB彩图
for  i=1:imsize
    for j=1:imsize           
        if (T(i,j)==1)     % Red
            RGBresult(i,j,1)=255;  % R分量
            RGBresult(i,j,2)=0;    % G分量
            RGBresult(i,j,3)=0;    % B分量   
        elseif(T(i,j)==2)  % Green
            RGBresult(i,j,1)=0; 
            RGBresult(i,j,2)=255; 
            RGBresult(i,j,3)=0;        
        elseif(T(i,j)==3)  % Blue
            RGBresult(i,j,1)=0; 
            RGBresult(i,j,2)=0; 
            RGBresult(i,j,3)=255;        
        end
    end
end
RGBresult=uint8(RGBresult);
imwrite(RGBresult,'3Bresult.bmp');

figure(1),imshow(Image);
figure,imshow(RGBresult);