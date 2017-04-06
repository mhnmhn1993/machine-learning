function [ output_args ] = FCM(IM)
%tmp=imread('junzhi.jpg');
tmp=imread('ppp.png');
tmp=tmp(:,:,1);
IM=double(tmp);
figure(3);
subplot(2,2,1);
imshow(uint8(IM));
title('图像增强处理后图像');
[maxX,maxY]=size(IM);
 
IMM=cat(4,IM,IM,IM,IM);

cc1=8;
cc2=80;
cc3=150;
cc4=220;
 
ttFcm=0;
 
while(ttFcm<31)
    ttFcm=ttFcm+1;
    
    c1=repmat(cc1,maxX,maxY);
    c2=repmat(cc2,maxX,maxY);
    c3=repmat(cc3,maxX,maxY);
    c4=repmat(cc4,maxX,maxY);
    c=cat(4,c1,c2,c3,c4);
    
    ree=repmat(0.000001,maxX,maxY);
    ree1=cat(4,ree,ree,ree,ree);
    
    distance=IMM-c;
    distance=distance.*distance+ree1;
    
    daoShu=1./distance;
    
    daoShu2=daoShu(:,:,:,1)+daoShu(:,:,:,2)+daoShu(:,:,:,3)+daoShu(:,:,:,4);
  
    distance1=distance(:,:,:,1).*daoShu2;
    u1=1./distance1;
    distance2=distance(:,:,:,2).*daoShu2;
    u2=1./distance2;
    distance3=distance(:,:,:,3).*daoShu2;
    u3=1./distance3;
    distance4=distance(:,:,:,4).*daoShu2;
    u4=1./distance4;
   
    ccc1=sum(sum(u1.*u1.*IM))/sum(sum(u1.*u1));
    ccc2=sum(sum(u2.*u2.*IM))/sum(sum(u2.*u2));
    ccc3=sum(sum(u3.*u3.*IM))/sum(sum(u3.*u3));
    ccc4=sum(sum(u4.*u4.*IM))/sum(sum(u4.*u4));
    
    tmpMatrix=[abs(cc1-ccc1)/cc1,abs(cc2-ccc2)/cc2,abs(cc3-ccc3)/cc3,abs(cc4-ccc4)/cc4];
    
    pp=cat(4,u1,u2,u3,u4);
    
    for i=1:maxX
        for j=1:maxY
            if max(pp(i,j,:))==u1(i,j)
                IX2(i,j)=1;
            elseif max(pp(i,j,:))==u2(i,j)
                IX2(i,j)=2;
            elseif max(pp(i,j,:))==u3(i,j)
                IX2(i,j)=3;
            else
                IX2(i,j)=4;
            end
        end
    end
   
    if max(tmpMatrix)<0.0001
        break;
    else
        cc1=ccc1;
        cc2=ccc2;
        cc3=ccc3;
        cc4=ccc4;
    end
    
    for i=1:maxX
        for j=1:maxY
            if IX2(i,j)==4
                IMMM(i,j)=255;
            elseif IX2(i,j)==3;
                IMMM(i,j)=180;
            elseif IX2(i,j)==2;
                IMMM(i,j)=40;
            else
                IMMM(i,j)=0;
            end
        end
    end
  figure(3);
    subplot(2,2,2);
    imshow(uint8(IMMM));
    title('处理后的灰度图');
end
 
for i=1:maxX
    for j=1:maxY
        if IX2(i,j)==4
            IMMM(i,j)=255;
        elseif IX2(i,j)==3
            IMMM(i,j)=180;
        elseif IX2(i,j)==2
            IMMM(i,j)=40;
        else
            IMMM(i,j)=0;
        end
    end
end

IMMM=uint8(IMMM);
colorbar;
figure(3);
subplot(2,2,3.5);
imshow(IMMM);
title('图像分割后'); 

figure
imshow(IMMM)
title('图像分割结果')
end



