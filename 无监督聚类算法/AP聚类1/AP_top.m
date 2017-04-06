clear
load('aggregation');
[n,~]=size(data);
%利用欧氏距离的相反数作为相似度矩阵
S = -squareform(pdist(data));
[row,col,v]=find(S);%找到不为0的元素
%v是一列n*(n-1)个元素的非0向量，排序之后取中位数
v=sortrows(v',1)';
med = median(v)-170;
for i=1:n
    S(i,i)=med;
end
idx = AP(S);
figure
clf
hold on
id_num = 0;
%共id_num个聚类
for i=unique(idx)'
    %找到第i类中的点，序号记为ii
    ii=find(idx==i);
    h=plot(data(ii,1),data(ii,2),'o');
    hold on;
    col=rand(1,3);
    set(h,'Color',col,'MarkerFaceColor',col,'MarkerSize',2);
    %以聚类中心为起点，其余点为终点，画出聚类内的连线
     xi1=data(i,1)*ones(size(ii)); xi2=data(i,2)*ones(size(ii));
     line([data(ii,1),xi1],[data(ii,2),xi2],'Color',col,'LineWidth',2);
    
    h1=plot(data(i,1),data(i,2),'o');
    set(h1,'Color',col,'MarkerEdgeColor','k','MarkerSize',6);
    id_num = id_num + 1;
end;
title({'AP聚类';...
    ['共',num2str(id_num),'个聚类']});