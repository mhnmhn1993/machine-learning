function Test_AP_ImgSegment()
clc;close all;clear all;
%说明：用AP聚类的方法进行图像分割

%1 读出待分割的图像 359146_b680c62513.jpg 1122267216_cbef8498a7.jpg
%1176827894_441080cda4.jpg  3011639697_6f98ec3d06.jpg
ImgPath='D:\福特安第二轮测试系统\Ldx图像分类与标注系统2.0\测试图像集\海边风景\';


Img=imread([ImgPath '359146_b680c62513.jpg']);  %444918_ba43b1675b.jpg  8251_917c7dd5ec
%   Img=imread([ImgPath '118310_3deb29c197.jpg']);
  Img=imread([ImgPath '184535991_978afa7e4e.jpg']);


%  ImgPath='D:\福特安第二轮测试系统\Ldx图像分类与标注系统2.0\测试图像集\园林风景\';
%  Img=imread([ImgPath '235898503_48d1223f34.jpg']);  %291067_a34b1e9d20  235898503_48d1223f34.jpg


[N M dim]=size(Img);
mrows=floor(N/2);
ncols=floor(M/2);

Img = imresize(Img,[mrows ncols],'bicubic');
figure,imshow(Img),title('原始图像')

%2 提取图像的特征
BoxSize=10;
[N M dim]=size(Img);
H=floor(N/BoxSize);W=floor(M/BoxSize);
FeaSet=[];Num=1;
FeaBuff=[];
for n=1:BoxSize:N
    for m=1:BoxSize:M
        y=n;x=m;
        if (y+BoxSize-1>N)
            y=N-BoxSize+1;
        end
        if (x+BoxSize-1>M)
            x=M-BoxSize+1;
        end
        Box=Img(y:y+BoxSize-1,x:x+BoxSize-1,:);
        
        %(1)颜色特征
        hsvBox=256*rgb2hsv(Box);
        clear mu ta
        for i=1:3
            temp=hsvBox(:,:,i);
            temp=temp(:);
            mu(i)=mean(temp); %均值
            ta(i)=sqrt(sum((temp-mu(i)).^2)/N); %方差
            S=sum((temp-mu(i)).^3)/N;
            if S<0
               is(i)=-((-S).^(1/3));  %二阶方差
            else
               is(i)=(S).^(1/3);
            end
        end
        FeaColor=[mu ta is];
        %Fea=Fea/sum(Fea);
        
        %(2)小波纹理
        origSize = size(Box);
        if length(origSize)==3 Im=rgb2gray(Box);end
        Im=double(Im);
        nbcol = 256;%size(colormap,1);
        [cA1,cH1,cV1,cD1]=dwt2(Im,'db1');
        tmpF1=mean(cA1(:));  tmpF2=std(cA1(:));
        tmpF3=mean(cH1(:));  tmpF4=std(cH1(:));
        tmpF5=mean(cV1(:));  tmpF6=std(cV1(:));
        tmpF7=mean(cD1(:));  tmpF8=std(cD1(:));
        WaveLetFea=[tmpF3 tmpF4 tmpF5 tmpF6 tmpF7 tmpF8];
          
        %(3)LBP纹理
        mapping=getmapping(8,'u2'); 
        LBPFea=LBP(rgb2gray(Box),1,8,mapping,'h'); %归一化的直方图LBP histogram in (8,1) neighborhood

        Fea=[FeaColor LBPFea WaveLetFea];
        Fea=[Fea];
        %Fea=[Fea y/N x/M];
        %记下块信息
        FeaSet(Num).boxsize=BoxSize;
        FeaSet(Num).y=y;
        FeaSet(Num).x=x;
        %FeaSet(Num).fea=Fea;
        FeaBuff=[FeaBuff;Fea];

        Num=Num+1;
    end %m loop
end
%特征归一化
% Tmax=max(FeaBuff);Tmax=repmat(Tmax,length(FeaSet),1);
% Tmin=min(FeaBuff);Tmin=repmat(Tmin,length(FeaSet),1);
% FeaBuff=(FeaBuff-Tmin)./(Tmax-Tmin);

%3 AP聚类
%3.1 计算相似性
N=length(FeaSet);
DisMat=zeros(N,N);
TSum=0;TNum=0;
for n=1:N
    for m=n+1:N
        %DisMat(n,m)=10/sqrt(sum((InstSet(n,:)-InstSet(m,:)).^2));
        D=sum( (FeaBuff(n,:)-FeaBuff(m,:) ).^2) ;
        DD=( FeaSet(n).y-FeaSet(m).y ).^2+( FeaSet(n).x-FeaSet(m).x ).^2;
        DD=sqrt(DD);
        D=sqrt(D);
        DisMat(n,m)=-(D+DD); %exp(-D/5);
        DisMat(m,n)=DisMat(n,m);
        TSum=TSum+DisMat(m,n)*2;TNum=TNum+2;
    end
end
% TAvg=TSum/TNum;
% for n=1:N
%    p(n)=TAvg;
% end
%开始聚类
%p=median(DisMat); %中位数
p=mean(DisMat); %中位数

%计算p最优范围
% [pmin,pmax]=preferenceRange(DisMat,'exact');
% for n=1:N
%    p(n)=(pmin+pmax)/2;
% end
% P=pmin
%方法一：由p控制聚成自适应类
%[idx,netsim,dpsim,expref]=apcluster(DisMat,p,'plot',0);
[idx,netsim,dpsim,expref] = apclustermex(DisMat,p,'plot',0);

%方法二：聚成指定的K类
%K=20
%[idx,netsim,dpsim,expref,pref]=apclusterK(DisMat,K);

%3 将同一聚类中的标上一个号
clear DisMat
C=unique(idx);
[H W dim]=size(Img);
Label=zeros(H,W);
Center=[];
for n=1:length(C)
    [index v]=find(idx==C(n));
    F=[];
    for m=1:length(index)
        P=index(m);
        Label( FeaSet(P).y: FeaSet(P).y+ FeaSet(P).boxsize-1,FeaSet(P).x: FeaSet(P).x+ FeaSet(P).boxsize-1)=n;
        F=[F;FeaBuff(P,:)];
    end
    F=mean(F,1);
    Center=[Center;F];
end
DrawBox(Img,Label),title('最初分割结果')
clear FeaBuff

%4 用区域生长法，对标号进行重新分配
[ImgSegLabel,FeaBuff]=RegionGraw(Center,Label);
DrawBox(Img,ImgSegLabel),title('区域生长后分割结果')

%5 相似性合并
%相似阈值
TSum=0;TNum=0;
for n=1:size(Center,1)
    for m=n+1:size(Center,1)
        D=sum( (Center(n,:)-Center(m,:) ).^2) ;
        D=sqrt(D);
        TSum=TSum+D;TNum=TNum+1;
    end
end
Th=TSum/TNum;
%合并
[TempImgSegLabel,Fea]=MergeBaseFeaSim(ImgSegLabel,FeaBuff,Th*0.6);
DrawBox(Img,TempImgSegLabel),title('相似合并后分割结果')

%6 将小区域合成进大区域
%[TempImgSegLabel]=MergeBaseReginNew(TempImgSegLabel,Fea); %将相似的小区域合并成一大区域
[TempImgSegLabel]=MergeBaseRegin(TempImgSegLabel,Fea);     %将小区域合并到与它相似且最相似的区域中

DrawBox(Img,TempImgSegLabel),title('小区域合并后分割结果')
max(TempImgSegLabel(:))

disp('Cluster is over...')

[N M]=size(TempImgSegLabel);
for n=1:max(TempImgSegLabel(:))
    [index v]=find(TempImgSegLabel==n);
    figure(1)
    Im=Img;
    for r=1:N
        for c=1:M
            if not(TempImgSegLabel(r,c)==n)
                Im(r,c,:)=0;Im(r,c,1)=250;
            end
        end
    end
    imshow(Im),title(['分割的第' num2str(n) '个区域'])
    pause(3)
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%说明6: 如果分出来的区域太小了，则合并 
%参数: ImgSegLabel--图像的结构标号 JoinTable--区域邻接表 
%      everyLabelPixNum - 每类的像素数目
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [TempImgSegLabel]=MergeBaseReginNew(ImgSegLabel,Fea)
    [N M]=size(ImgSegLabel);
    MaxLabel=max(ImgSegLabel(:));
    
    JoinTable=zeros(MaxLabel,MaxLabel);
    everyLabelPixNum=zeros(1,MaxLabel);
    
    %建立邻接表与统计每个区的像素个数
    for r=1:N
        for c=1:M
           t=ImgSegLabel(r,c);
           everyLabelPixNum(t)=everyLabelPixNum(t)+1;
           for rr=-1:1
               for cc=-1:1
                   if not(rr*cc==0)
                      continue; 
                   end
                   y=r+rr;x=c+cc;
                   if (y>0 && x>0 && y<=N && x<=M) %下标不要出界
                       t1=ImgSegLabel(y,x);
                       if (not(t1==t))
                           JoinTable(t,t1)=JoinTable(t,t1)+1; %有了一次邻接关系 
                           JoinTable(t1,t)=JoinTable(t1,t)+1; 
                       end
                   end
               end
           end
        end %c loop
    end %r loop
    
    [v idx]=min(everyLabelPixNum);
    %h=Fea(idx,:);
    while (v<1000)
        %找相邻且最相似的
        tSim=inf;
        ptr=0;
        for n=1:length(everyLabelPixNum) %找最相似的
            if ((JoinTable(idx,n)>0 || JoinTable(n,idx)>0 )  && not(idx==n))
                 %h1=Fea(n,:);
                 tSim1=everyLabelPixNum(n); %sum((h-h1).^2);
                 if (tSim1<tSim)
                     tSim=tSim1;ptr=n;
                 end
            end
        end
        if (ptr==0)
           a=1 
        end
        
        %将 idx 类合并到 ptr 类
        idx1=find(ImgSegLabel==idx);
        ImgSegLabel(idx1)=ptr;

        JoinTable(ptr,:)=JoinTable(ptr,:)+JoinTable(idx,:);
        JoinTable(:,ptr)=JoinTable(:,ptr)+JoinTable(:,idx);
        JoinTable(idx,:)=0;
        JoinTable(:,idx)=0;

        everyLabelPixNum(ptr)=everyLabelPixNum(ptr)+everyLabelPixNum(idx);
        everyLabelPixNum(idx)=inf;
        [v idx]=min(everyLabelPixNum);
    end %while
    
    %将不存的标号去掉
    TempImgSegLabel=zeros(N,M);
    tLebal=1;
    for r=1:MaxLabel
        if not(everyLabelPixNum(r)==inf)
            idx=find(ImgSegLabel==r);
            TempImgSegLabel(idx)=tLebal;
            %everyLabelCHist(tLebal,:)=everyLabelCHist(r,:);
            everyLabelPixNum(tLebal)=everyLabelPixNum(r);
            tLebal=tLebal+1;
        end
    end

    disp('基于大小合并完成')
    return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%说明6: 如果分出来的区域太小了，则合并 
%参数: ImgSegLabel--图像的结构标号 JoinTable--区域邻接表 
%      everyLabelPixNum - 每类的像素数目
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [TempImgSegLabel]=MergeBaseRegin(ImgSegLabel,Fea)
    [N M]=size(ImgSegLabel);
    MaxLabel=max(ImgSegLabel(:));
    
    JoinTable=zeros(MaxLabel,MaxLabel);
    everyLabelPixNum=zeros(1,MaxLabel);
    
    %建立邻接表与统计每个区的像素个数
    for r=1:N
        for c=1:M
           t=ImgSegLabel(r,c);
           everyLabelPixNum(t)=everyLabelPixNum(t)+1;
           for rr=-1:1
               for cc=-1:1
                   if not(rr*cc==0)
                      continue; 
                   end
                   y=r+rr;x=c+cc;
                   if (y>0 && x>0 && y<=N && x<=M) %下标不要出界
                       t1=ImgSegLabel(y,x);
                       if (not(t1==t))
                           JoinTable(t,t1)=JoinTable(t,t1)+1; %有了一次邻接关系 
                           JoinTable(t1,t)=JoinTable(t1,t)+1; 
                       end
                   end
               end
           end
        end %c loop
    end %r loop    
    
    [v idx]=min(everyLabelPixNum);
    h=Fea(idx,:);
    while (v<1000)
        %找相邻且最相似的
        tSim=inf;
        ptr=0;
        for n=1:length(everyLabelPixNum) %找最相似的
            if ((JoinTable(idx,n)>0 || JoinTable(n,idx)>0 )  && not(idx==n))
                 h1=Fea(n,:);
                 tSim1=sum((h-h1).^2);
                 if (tSim1<tSim)
                     tSim=tSim1;ptr=n;
                 end
            end
        end
        if (ptr==0)
           disp('合并时出现错误') 
        end
        
        %将 idx 类合并到 ptr 类
        idx1=find(ImgSegLabel==idx);
        ImgSegLabel(idx1)=ptr;

        JoinTable(ptr,:)=JoinTable(ptr,:)+JoinTable(idx,:);
        JoinTable(:,ptr)=JoinTable(:,ptr)+JoinTable(:,idx);
        JoinTable(idx,:)=0;
        JoinTable(:,idx)=0;

        everyLabelPixNum(ptr)=everyLabelPixNum(ptr)+everyLabelPixNum(idx);
        everyLabelPixNum(idx)=inf;
        [v idx]=min(everyLabelPixNum);
    end %while
    
    %将不存的标号去掉
    TempImgSegLabel=zeros(N,M);
    tLebal=1;
    for r=1:MaxLabel
        if not(everyLabelPixNum(r)==inf)
            idx=find(ImgSegLabel==r);
            TempImgSegLabel(idx)=tLebal;
            %everyLabelCHist(tLebal,:)=everyLabelCHist(r,:);
            everyLabelPixNum(tLebal)=everyLabelPixNum(r);
            tLebal=tLebal+1;
        end
    end
    disp('利用相似对小区域进行合并完成')



function DrawBox(Img,ImgSegLabel)
    %在分界处画上线　以观看
    [N M tp]=size(Img);
    tBW=zeros(N,M);
    II=Img;
    for r=2:N-1
        for c=2:M-1
            if (not(ImgSegLabel(r,c)==ImgSegLabel(r,c-1)) || not(ImgSegLabel(r,c)==ImgSegLabel(r-1,c)))
                tBW(r,c)=1;
                II(r-1:r+1,c-1:c+1,:)=0; II(r-1:r+1,c-1:c+1,1)=255;
            end
        end
    end
    %figure,imshow(tBW),title('平滑与非平滑区分割后边界图')
    figure,imshow(II),title('平滑与非平滑区分割图')
  



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%说明5: 在分出来各区上基于颜色的相似性合并 
%参数: ImgSegLabel--图像的结构标号 everyLabelCHist--图像的颜色直方图 everyLabelPixNum
%      everyLabelPixNum - 每类的像素数目
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [TempImgSegLabel,tFea]=MergeBaseFeaSim(ImgSegLabel,Fea,cTh)
    [N M]=size(ImgSegLabel);
    MaxLabel=max(ImgSegLabel(:));
    JoinTable=zeros(MaxLabel,MaxLabel);    %邻接关系表
    everyLabelPixNum=zeros(1,MaxLabel);
    everyLabelCHist=zeros(MaxLabel,9);
    
    %建立邻接表与统计每个区的像素个数
    for r=1:N
        for c=1:M
           t=ImgSegLabel(r,c);
           everyLabelPixNum(t)=everyLabelPixNum(t)+1;
           if (c+1<=M)
              t1=ImgSegLabel(r,c+1);
              if (not(t1==t))
                  JoinTable(t,t1)=JoinTable(t,t1)+1; 
                  JoinTable(t1,t)=JoinTable(t1,t)+1;%有了一次邻接关系 
              end
           end
           if (r+1<=N)
              t1=ImgSegLabel(r+1,c);
              if (not(t1==t))
                  JoinTable(t,t1)=JoinTable(t,t1)+1; 
                  JoinTable(t1,t)=JoinTable(t1,t)+1; %有了一次邻接关系 
              end
           end
%            if (r+1<=N && c+1<=M)
%               t1=ImgSegLabel(r+1,c+1);
%               if (not(t1==t))
%                   JoinTable(t,t1)=JoinTable(t,t1)+1; JoinTable(t1,t)=JoinTable(t1,t)+1; %有了一次邻接关系 
%               end
%            end
           
        end %c loop
    end %r loop

    %合并%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %cTh=0.2;          %相似阈值
    FeaBuff=[];
    for r=1:MaxLabel
        if (everyLabelPixNum(r)<=0)
            continue;        
        end
        OK_flg=1;
        F=Fea(r,:);
        while (OK_flg)  %将所颜色相似全合并
            tSim=inf;ptr=0;h=Fea(r,:);
            for c=1:MaxLabel %找最相似的
                if ((JoinTable(r,c)>0 || JoinTable(c,r)>0) && not(c==r)) % && everyLabelPixNum(r,:)>0)
                    h1=Fea(c,:);
                    tSim1=sum( (h-h1).^2 );
                    tSim1=sqrt(tSim1);
                    if (tSim1<tSim)
                        tSim=tSim1;ptr=c;
                    end
                end
            end

            if (tSim<cTh) % && not(ptr==r))  %颜色相似
                %将 ptr 类合并到 pr 类
                idx=find(ImgSegLabel==ptr);
                ImgSegLabel(idx)=r;

                JoinTable(r,:)=JoinTable(r,:)+JoinTable(ptr,:);
                JoinTable(:,r)=JoinTable(:,r)+JoinTable(:,ptr);
                JoinTable(ptr,:)=0;
                JoinTable(:,ptr)=0;

                everyLabelPixNum(r)=everyLabelPixNum(r)+everyLabelPixNum(ptr);
                everyLabelPixNum(ptr)=0;
                F=[F;Fea(ptr,:)];
            else
                OK_flg=0;
            end
        end %while
        FeaBuff(r,:)=mean(F,1);
    end % r loop

    %将不存的标号去掉
    TempImgSegLabel=zeros(N,M);
    tLebal=1;
    tFea=[];
    for r=1:MaxLabel
        if (everyLabelPixNum(r)>0)
            idx=find(ImgSegLabel==r);
            TempImgSegLabel(idx)=tLebal;
            tFea(tLebal,:)=FeaBuff(r,:);
            tLebal=tLebal+1;
        end
    end
    clear FeaBuff
    disp('基于相似合并完成')
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%说明4: 利用区域生长法，对聚类的标号进行重新分配
%参数: ImgStructLabel--图像的结构标号 ImgColorHist--图像的颜色直方图  ImgSegLabel- 分割标号
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ImgSegLabel]=Segmentation(Fea,N,M,Th)
    MaxLabel=1;
    ImgSegLabel=zeros(N,M);
    Visted=zeros(N,M);  %访问标志
    
    %3.1 区域生长
    for r=2:N-1
        for c=2:M-1
           if (Visted(r,c)==0) %没有访问过，且是光滑的
               %它周围的4个像素相似，用周围的标它
               if (Visted(r-1,c)==1 && ImgSegLabel(r-1,c)==ImgSegLabel(r,c-1) && ImgSegLabel(r-1,c)==ImgSegLabel(r,c+1))
                   t=ImgSegLabel(r-1,c);
                   ImgSegLabel(r,c)=t;Visted(r,c)=1;
                   continue;
               else
                   ImgSegLabel(r,c)=MaxLabel;Visted(r,c)=1;
                   Q=[r c]; %入队
               end

               while (length(Q)>0)
                   y0=Q(1);x0=Q(2); %取队头
                   cHist0=Fea(y0,x0,:); %everyLabelCHist(MaxLabel,:)/everyLabelPixNum(MaxLabel,:); %区内的平均颜色
                   for rr=-1:1
                       for cc=-1:1
                           y=y0+rr;x=x0+cc;
                           if (y>0 && x>0 && y<N && x<M) %下标不要出界
                           if (Visted(y,x)==0)
                               %比较(y,x)(y0,x0)点的颜色相似性
                               cHist1=Fea(y,x,:); % ImgColorHist(y,x).cHist; 
                               D=sum((cHist1-cHist0).^2);
                               if (D<=Th)
                                   ImgSegLabel(y,x)=MaxLabel;Visted(y,x)=1;
                                   Q=[Q y x]; %入队
                               else
                                   %虽然它颜色不相似，但它周围的4个像素相似，用周围的标它
                                   if (y-1>0 && x-1>0 && x+1<M && Visted(y-1,x)==1 && ImgSegLabel(y-1,x)==ImgSegLabel(y,x-1) && ImgSegLabel(y-1,x)==ImgSegLabel(y,x+1))
                                        t=ImgSegLabel(y-1,x);
                                        ImgSegLabel(y,x)=t;Visted(y,x)=1;
                                   end
                               end %if
                           end % if
                           end %下标不要出界
                       end %cc loop
                   end %rr loop
                   tL=length(Q);
                   if (tL>2)
                       Q=Q(3:tL);
                   else
                       Q=[];
                   end
               end %while loop
               MaxLabel=MaxLabel+1;
           end %if
        end % c loop
    end % r loop
figure,
imagesc(ImgSegLabel); figure(gcf),title('区域生长后的结果')

    %边界处理
    for r=2:M-1
        if (Visted(1,r)==0)
            Visted(1,r)=1;ImgSegLabel(1,r)=ImgSegLabel(2,r);t=ImgSegLabel(2,r);
            %everyLabelCHist(t,:)=everyLabelCHist(t,:)+ImgColorHist(1,r).cHist;everyLabelPixNum(t,:)=everyLabelPixNum(t,:)+1;
        end
        if (Visted(N,r)==0)
            Visted(N,r)=1;ImgSegLabel(N,r)=ImgSegLabel(N-1,r);t=ImgSegLabel(N-1,r);
            %everyLabelCHist(t,:)=everyLabelCHist(t,:)+ImgColorHist(N-1,r).cHist;everyLabelPixNum(t,:)=everyLabelPixNum(t,:)+1;
        end
    end
    
    for r=1:N
        if (Visted(r,1)==0)
            Visted(r,1)=1;ImgSegLabel(r,1)=ImgSegLabel(r,2);t=ImgSegLabel(r,2);
            %everyLabelCHist(t,:)=everyLabelCHist(t,:)+ImgColorHist(r,1).cHist;everyLabelPixNum(t,:)=everyLabelPixNum(t,:)+1;
        end
        if (Visted(r,M)==0)
            Visted(r,M)=1;ImgSegLabel(r,M)=ImgSegLabel(r,M-1);t=ImgSegLabel(r,M-1);
            %everyLabelCHist(t,:)=everyLabelCHist(t,:)+ImgColorHist(r,M-1).cHist;everyLabelPixNum(t,:)=everyLabelPixNum(t,:)+1;
        end
    end

   disp('第一次分割完成')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%说明4: 利用区域生长法，对聚类的标号进行重新分配
%参数: Label--图像聚类得到的标号 Fea--图像的特征  ImgSegLabel,FeaBuff- 重新分配的标号与特征
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ImgSegLabel,FeaBuff]=RegionGraw(Fea,Label)
    MaxLabel=1;
    [N M]=size(Label);
    ImgSegLabel=zeros(N,M);
    Visted=zeros(N,M);  %访问标志
    FeaBuff=[];
    %3.1 对平滑区进行分割
    for r=1:N
        for c=1:M
           if (Visted(r,c)==0) %没有访问过
               ImgSegLabel(r,c)=MaxLabel;
               Visted(r,c)=1;
               Q=[r c]; %入队
               y0=Q(1);x0=Q(2); %取队头
               tt0=Label(y0,x0);
               
               while (length(Q)>0)
                   y0=Q(1);x0=Q(2); %取队头
                   t0=Label(y0,x0);
                   %它的4领域
                   for rr=-1:1
                       for cc=-1:1
                           if not(rr*cc==0)  %四领域
                              continue;; 
                           end
                           y=y0+rr;x=x0+cc;
                           if (y>0 && x>0 && y<=N && x<=M) %下标不要出界
                               if (Visted(y,x)==0)
                                   t1=Label(y,x);
                                   if (t1==t0)
                                       ImgSegLabel(y,x)=MaxLabel;
                                       Visted(y,x)=1;
                                       Q=[Q y x]; %入队
                                   end %if
                               end % if
                           end %下标不要出界
                       end %cc loop
                   end %rr loop
                   tL=length(Q);
                   if (tL>2)
                       Q=Q(3:tL);
                   else
                       Q=[];
                   end
               end %while loop
               FeaBuff=[FeaBuff;Fea(tt0,:)];
               MaxLabel=MaxLabel+1;
           end %if
        end % c loop
    end % r loop
    % figure,
    % imagesc(ImgSegLabel); figure(gcf),title('区域生长后的结果')
    disp('区域生长完成...')   

%APCLUSTER Affinity Propagation Clustering (Frey/Dueck, Science 2007) 
% [idx,netsim,dpsim,expref]=APCLUSTER(s,p) clusters data, using a set  
% of real-valued pairwise data point similarities as input. Clusters  
% are each represented by a cluster center data point (the "exemplar").  
% The method is iterative and searches for clusters so as to maximize  
% an objective function, called net similarity. 
%  
% For N data points, there are potentially N^2-N pairwise similarities;  
% this can be input as an N-by-N matrix 's', where s(i,k) is the  
% similarity of point i to point k (s(i,k) needn?t equal s(k,i)).  In  
% fact, only a smaller number of relevant similarities are needed; if  
% only M similarity values are known (M < N^2-N) they can be input as  
% an M-by-3 matrix with each row being an (i,j,s(i,j)) triple. 
%  
% APCLUSTER automatically determines the number of clusters based on  
% the input preference 'p', a real-valued N-vector. p(i) indicates the  
% preference that data point i be chosen as an exemplar. Often a good  
% choice is to set all preferences to median(s); the number of clusters  
% identified can be adjusted by changing this value accordingly. If 'p'  
% is a scalar, APCLUSTER assumes all preferences are that shared value. 
%  
% The clustering solution is returned in idx. idx(j) is the index of  
% the exemplar for data point j; idx(j)==j indicates data point j  
% is itself an exemplar. The sum of the similarities of the data points to  
% their exemplars is returned as dpsim, the sum of the preferences of  
% the identified exemplars is returned in expref and the net similarity  
% objective function returned is their sum, i.e. netsim=dpsim+expref. 
%  
% 	[ ... ]=apcluster(s,p,'NAME',VALUE,...) allows you to specify  
% 	  optional parameter name/value pairs as follows: 
%  
%   'maxits'     maximum number of iterations (default: 1000) 
%   'convits'    if the estimated exemplars stay fixed for convits  
%          iterations, APCLUSTER terminates early (default: 100) 
%   'dampfact'   update equation damping level in [0.5, 1).  Higher  
%        values correspond to heavy damping, which may be needed  
%        if oscillations occur. (default: 0.9) 
%   'plot'       (no value needed) Plots netsim after each iteration 
%   'details'    (no value needed) Outputs iteration-by-iteration  
%      details (greater memory requirements) 
%   'nonoise'    (no value needed) APCLUSTER adds a small amount of  
%      noise to 's' to prevent degenerate cases; this disables that. 
%  
% Copyright (c) B.J. Frey & D. Dueck (2006). This software may be  
% freely used and distributed for non-commercial purposes. 
%          (RUN APCLUSTER WITHOUT ARGUMENTS FOR DEMO CODE) 
function [idx,netsim,dpsim,expref]=apcluster(s,p,varargin); 
if nargin==0, % display demo 
	fprintf('Affinity Propagation (APCLUSTER) sample/demo code\n\n'); 
	fprintf('N=100; x=rand(N,2); % Create N, 2-D data points\n'); 
	fprintf('M=N*N-N; s=zeros(M,3); % Make ALL N^2-N similarities\n'); 
	fprintf('j=1;\n'); 
	fprintf('for i=1:N\n'); 
	fprintf('  for k=[1:i-1,i+1:N]\n'); 
	fprintf('    s(j,1)=i; s(j,2)=k; s(j,3)=-sum((x(i,:)-x(k,:)).^2);\n'); 
	fprintf('    j=j+1;\n'); 
	fprintf('  end;\n'); 
	fprintf('end;\n'); 
	fprintf('p=median(s(:,3)); % Set preference to median similarity\n'); 
	fprintf('[idx,netsim,dpsim,expref]=apcluster(s,p,''plot'');\n'); 
	fprintf('fprintf(''Number of clusters: %%d\\n'',length(unique(idx)));\n'); 
	fprintf('fprintf(''Fitness (net similarity): %%g\\n'',netsim);\n'); 
	fprintf('figure; % Make a figures showing the data and the clusters\n'); 
	fprintf('for i=unique(idx)''\n'); 
	fprintf('  ii=find(idx==i); h=plot(x(ii,1),x(ii,2),''o''); hold on;\n'); 
	fprintf('  col=rand(1,3); set(h,''Color'',col,''MarkerFaceColor'',col);\n'); 
	fprintf('  xi1=x(i,1)*ones(size(ii)); xi2=x(i,2)*ones(size(ii)); \n'); 
	fprintf('  line([x(ii,1),xi1]'',[x(ii,2),xi2]'',''Color'',col);\n'); 
	fprintf('end;\n'); 
	fprintf('axis equal tight;\n\n'); 
	return; 
end; 
start = clock; 
% Handle arguments to function 
if nargin<2 error('Too few input arguments'); 
else 
    maxits=1000; convits=100; lam=0.8; plt=0; details=0; nonoise=0; 
    %maxits=1000; convits=300; lam=0.6; plt=0; details=0; nonoise=0; 
    
    i=1; 
    while i<=length(varargin) 
        if strcmp(varargin{i},'plot') 
            plt=1; i=i+1; 
        elseif strcmp(varargin{i},'details') 
            details=1; i=i+1; 
		elseif strcmp(varargin{i},'sparse') 
% 			[idx,netsim,dpsim,expref]=apcluster_sparse(s,p,varargin{:}); 
			fprintf('''sparse'' argument no longer supported; see website for additional software\n\n'); 
			return; 
        elseif strcmp(varargin{i},'nonoise') 
            nonoise=1; i=i+1; 
        elseif strcmp(varargin{i},'maxits') 
            maxits=varargin{i+1}; 
            i=i+2; 
            if maxits<=0 error('maxits must be a positive integer'); end; 
        elseif strcmp(varargin{i},'convits') 
            convits=varargin{i+1}; 
            i=i+2; 
            if convits<=0 error('convits must be a positive integer'); end; 
        elseif strcmp(varargin{i},'dampfact') 
            lam=varargin{i+1}; 
            i=i+2; 
            if (lam<0.5)||(lam>=1) 
                error('dampfact must be >= 0.5 and < 1'); 
            end; 
        else i=i+1; 
        end; 
    end; 
end; 
if lam>0.9 
    fprintf('\n*** Warning: Large damping factor in use. Turn on plotting\n'); 
    fprintf('    to monitor the net similarity. The algorithm will\n'); 
    fprintf('    change decisions slowly, so consider using a larger value\n'); 
    fprintf('    of convits.\n\n'); 
end; 
 
% Check that standard arguments are consistent in size 
if length(size(s))~=2 error('s should be a 2D matrix'); 
elseif length(size(p))>2 error('p should be a vector or a scalar'); 
elseif size(s,2)==3 
    tmp=max(max(s(:,1)),max(s(:,2))); 
    if length(p)==1 N=tmp; else N=length(p); end; 
    if tmp>N 
        error('data point index exceeds number of data points'); 
    elseif min(min(s(:,1)),min(s(:,2)))<=0 
        error('data point indices must be >= 1'); 
    end; 
elseif size(s,1)==size(s,2) 
    N=size(s,1); 
    if (length(p)~=N)&&(length(p)~=1) 
        error('p should be scalar or a vector of size N'); 
    end; 
else error('s must have 3 columns or be square'); end; 
 
% Construct similarity matrix 
if N>3000 
    fprintf('\n*** Warning: Large memory request. Consider activating\n'); 
    fprintf('    the sparse version of APCLUSTER.\n\n'); 
end; 
if size(s,2)==3 && size(s,1)~=3, 
    S=-Inf*ones(N,N,class(s));  
    for j=1:size(s,1), S(s(j,1),s(j,2))=s(j,3); end; 
else S=s; 
end; 
 
if S==S', symmetric=true; else symmetric=false; end; 
realmin_=realmin(class(s)); realmax_=realmax(class(s)); 
 
% In case user did not remove degeneracies from the input similarities, 
% avoid degenerate solutions by adding a small amount of noise to the 
% input similarities 
if ~nonoise 
    rns=randn('state'); randn('state',0); 
    S=S+(eps*S+realmin_*100).*rand(N,N); 
    randn('state',rns); 
end; 
 
% Place preferences on the diagonal of S 
if length(p)==1 for i=1:N S(i,i)=p; end; 
else for i=1:N S(i,i)=p(i); end; 
end; 
 
% Numerical stability -- replace -INF with -realmax 
n=find(S<-realmax_); if ~isempty(n), warning('-INF similarities detected; changing to -REALMAX to ensure numerical stability'); S(n)=-realmax_; end; clear('n'); 
if ~isempty(find(S>realmax_,1)), error('+INF similarities detected; change to a large positive value (but smaller than +REALMAX)'); end; 
 
 
% Allocate space for messages, etc 
dS=diag(S); A=zeros(N,N,class(s)); R=zeros(N,N,class(s)); t=1; 
if plt, netsim=zeros(1,maxits+1); end; 
if details 
    idx=zeros(N,maxits+1); 
    netsim=zeros(1,maxits+1);  
    dpsim=zeros(1,maxits+1);  
    expref=zeros(1,maxits+1);  
end; 
 
% Execute parallel affinity propagation updates 
e=zeros(N,convits); dn=0; i=0; 
if symmetric, ST=S; else ST=S'; end; % saves memory if it's symmetric 
while ~dn 
    i=i+1;  
 
    % Compute responsibilities 
	%%A=A'; R=R';
    T=A+R;
    A=A';  
    R=R'; 
	for ii=1:N, 
		old = R(:,ii); 
		AS = A(:,ii) + ST(:,ii); [Y,I]=max(AS); AS(I)=-Inf; 
		[Y2,I2]=max(AS); 
		R(:,ii)=ST(:,ii)-Y; 
		R(I,ii)=ST(I,ii)-Y2; 
		R(:,ii)=(1-lam)*R(:,ii)+lam*old; % Damping 
        R(R(:,ii)>realmax_,ii)=realmax_; 
	end; 
	A=A'; R=R'; 
 
    % Compute availabilities 
	for jj=1:N, 
		old = A(:,jj); 
		Rp = max(R(:,jj),0); Rp(jj)=R(jj,jj); 
		A(:,jj) = sum(Rp)-Rp; 
		dA = A(jj,jj); A(:,jj) = min(A(:,jj),0); A(jj,jj) = dA; 
		A(:,jj) = (1-lam)*A(:,jj) + lam*old; % Damping 
	end; 
	 
    % Check for convergence 
    E=((diag(A)+diag(R))>0); e(:,mod(i-1,convits)+1)=E; K=sum(E); 
    if i>=convits || i>=maxits, 
        se=sum(e,2); 
        unconverged=(sum((se==convits)+(se==0))~=N); 
        if ( ~unconverged && (K>0) )||(i==maxits) dn=1; end; 
    end; 
 
    % Handle plotting and storage of details, if requested 
    if plt||details 
        if K==0 
            tmpnetsim=nan; tmpdpsim=nan; tmpexpref=nan; tmpidx=nan; 
        else 
            I=find(E); notI=find(~E); [tmp c]=max(S(:,I),[],2); c(I)=1:K; tmpidx=I(c); 
            tmpdpsim=sum(S(sub2ind([N N],notI,tmpidx(notI)))); 
            tmpexpref=sum(dS(I)); 
            tmpnetsim=tmpdpsim+tmpexpref; 
        end; 
    end; 
    if details 
        netsim(i)=tmpnetsim; dpsim(i)=tmpdpsim; expref(i)=tmpexpref; 
        idx(:,i)=tmpidx; 
    end; 
    if plt, 
        netsim(i)=tmpnetsim; 
		figure(234); 
        plot(((netsim(1:i)/10)*100)/10,'r-'); xlim([0 i]); % plot barely-finite stuff as infinite 
        xlabel('# Iterations'); 
        ylabel('Fitness (net similarity) of quantized intermediate solution'); 
%         drawnow;  
    end; 
end; % iterations 
I=find((diag(A)+diag(R))>0); K=length(I); % Identify exemplars 
if K>0 
    [tmp c]=max(S(:,I),[],2); c(I)=1:K; % Identify clusters 
    % Refine the final set of exemplars and clusters and return results 
    for k=1:K ii=find(c==k); [y j]=max(sum(S(ii,ii),1)); I(k)=ii(j(1)); end; notI=reshape(setdiff(1:N,I),[],1); 
    [tmp c]=max(S(:,I),[],2); c(I)=1:K; tmpidx=I(c); 
	tmpdpsim=sum(S(sub2ind([N N],notI,tmpidx(notI)))); 
	tmpexpref=sum(dS(I)); 
	tmpnetsim=tmpdpsim+tmpexpref; 
else 
    tmpidx=nan*ones(N,1); tmpnetsim=nan; tmpexpref=nan; 
end; 
if details 
    netsim(i+1)=tmpnetsim; netsim=netsim(1:i+1); 
    dpsim(i+1)=tmpdpsim; dpsim=dpsim(1:i+1); 
    expref(i+1)=tmpexpref; expref=expref(1:i+1); 
    idx(:,i+1)=tmpidx; idx=idx(:,1:i+1); 
else 
    netsim=tmpnetsim; dpsim=tmpdpsim; expref=tmpexpref; idx=tmpidx; 
end; 
if plt||details 
    fprintf('\nNumber of exemplars identified: %d  (for %d data points)\n',K,N); 
    fprintf('Net similarity: %g\n',tmpnetsim); 
    fprintf('  Similarities of data points to exemplars: %g\n',dpsim(end)); 
    fprintf('  Preferences of selected exemplars: %g\n',tmpexpref); 
    fprintf('Number of iterations: %d\n\n',i); 
	fprintf('Elapsed time: %g sec\n',etime(clock,start)); 
end; 
if unconverged 
	fprintf('\n*** Warning: Algorithm did not converge. Activate plotting\n'); 
	fprintf('    so that you can monitor the net similarity. Consider\n'); 
	fprintf('    increasing maxits and convits, and, if oscillations occur\n'); 
	fprintf('    also increasing dampfact.\n\n'); 
end;
