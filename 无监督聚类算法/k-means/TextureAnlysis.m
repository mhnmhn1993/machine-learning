clc
clear all;

Image=imread('3B.bmp');
imsize=size(Image,1);
winsize=4;
winnum=imsize/winsize;

ENTmean=zeros(winnum,winnum);    % �صľ�ֵ
% ENTvar=zeros(winnum,winnum);   % �صķ��������Ч��
% CONmean=zeros(winnum,winnum);  % �Աȶȶ����ľ�ֵ��������Ч��
% K1=2;K2=2;
stepsize=4;

for i=1:winsize:imsize
    for j=1:winsize:imsize
         Qmatrix=zeros(winsize,winsize);
         Pmatrix=zeros(stepsize,stepsize,4);
         
         block=Image(i:i+winsize-1,j:j+winsize-1);  % ͼ���ӿ�
         
         maxpxl=max(max(block));
         minpxl=min(min(block));
         step=(maxpxl-minpxl)/stepsize;
         
         for k=1:stepsize  % ����
             index=find(block>=minpxl+(k-1)*step & block<minpxl+k*step);
             Qmatrix(index)=k;
         end
         index=find(block>=minpxl+(k-1)*step & block<=maxpxl);
         Qmatrix(index)=k;  
         
         for m=1:stepsize   % �ĸ�����Ĺ�������0��45��90��135�ȣ�
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
            M(k)=(-1)*sum(sum(pxl.*log10(pxl))); % �ĸ��Ҷȹ����������
                        
%             for m=1:stepsize           % �ĸ��Ҷȹ�������ĶԱȶȶ���
%                 for n=1:stepsize
%                     pxlabs=abs(m-n);
%                     CON(1,k)=CON(1,k)+(pxlabs^K1)*(Pmatrix(m,n,k)^K2);
%                 end
%             end
        end
        
        ENTmean(floor(i/stepsize)+1,floor(j/stepsize)+1)=mean(M);               % �صľ�ֵ        
%         ENTvar(floor(i/stepsize)+1,floor(j/stepsize)+1)=mean((M-mean(M)).^2); % �صķ��� 
%         CONmean(floor(i/stepsize)+1,floor(j/stepsize)+1)=mean(CON);    % �Աȶȶ����ľ�ֵ           
    end
end  % end for

Texture=zeros(2,imsize,imsize);  % ��������

for i=1:winnum
    for j=1:winnum
        Texture(1,[(i-1)*winsize+1:i*winsize],[(i-1)*winsize+1:i*winsize])=ENTmean(i,j);
%         Texture(1,[(i-1)*winsize+1:i*winsize],[(i-1)*winsize+1:i*winsize])=ENTvar(i,j);  % ������Ч��
%         Texture(2,[(i-1)*winsize+1:i*winsize],[(i-1)*winsize+1:i*winsize])=CONmean(i,j); % ������Ч��
    end
end

Imext=wextend('2D','zpd',Image,[4,4]); % zpd ����  sp0  ƽ��
Direction=[1,0,0,0,1,0,0,0,1; ...   % 9��9�ľ��󣬰˸�����Ϊ1
    0,1,0,0,1,0,0,1,0; ...
    0,0,1,0,1,0,1,0,0; ...
    0,0,0,1,1,1,0,0,0; ...
    1,1,1,1,0,1,1,1,1; ...  % ���ĵ�ȡ0 1�����Բ��
    0,0,0,1,1,1,0,0,0; ...
    0,0,1,0,1,0,1,0,0; ...
    0,1,0,0,1,0,0,1,0; ...
    1,0,0,0,1,0,0,0,1];
Direction=double(Direction);
for i=1:imsize
    for j=1:imsize
%         pxlmean=mean(mean(Imext([i:i+2],[j:j+2])));
%         Imsmooth(i,j)=sum(sum((Imext([i:i+2],[j:j+2])-pxlmean).^2)); % ��λ����Ҷȱ仯������Ч������
%         Imsmooth(i,j)=mean(mean(Imext([i:i+4],[j:j+4])));       % ƽ���Ҷ�,Ч�����ڵ�λ����Ҷȱ仯����
        Direction2=double(Imext([i:i+8],[j:j+8]));
        Imsmooth(i,j)=mean(mean(Direction2.*Direction));  % �������ڰ˸�����ȡ��ֵ��Ӧ���бȾ�ֵ���õģ�
    end
end
Texture(2,:,:)=Imsmooth;           % ����þ�ֵƽ�������Ӳ�࣬��Ч����������

C1=Texture(:,30,30);               % ��ʼ��������
C2=Texture(:,150,200);
C3=Texture(:,250,400);

T=zeros(imsize,imsize);            % �����ʶ����
number1=0; number2=0; number3=0;   % �ۼƸ��������ظ���
e_new=10;e_old=20;

while abs(e_new-e_old)>1
    e_old=e_new;
    e_new=0;
   
    for m=1:imsize
        for n=1:imsize
            d1=sum((Texture(:,m,n)-C1).^2);
            d2=sum((Texture(:,m,n)-C2).^2);  
            d3=sum((Texture(:,m,n)-C3).^2);
            dmin=min([d1,d2,d3]); % ��С���� ��ƽ��
            
            if d1==dmin           % ������С�������
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
    
    C1=C1.*0; C2=C2.*0; C3=C3.*0; % �������㣬�Ժ����
    e_new=e_new/(imsize^2);       % �������
    
    for m=1:imsize                % ���¾�������
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
    C1=C1./number1; number1=0;    % �µľ�������
    C2=C2./number2; number2=0;
    C3=C3./number3; number3=0;       
end % end while

RGBresult=zeros(imsize,imsize,3); % �����ϲ�ɫ��RGB��ͼ
for  i=1:imsize
    for j=1:imsize           
        if (T(i,j)==1)     % Red
            RGBresult(i,j,1)=255;  % R����
            RGBresult(i,j,2)=0;    % G����
            RGBresult(i,j,3)=0;    % B����   
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