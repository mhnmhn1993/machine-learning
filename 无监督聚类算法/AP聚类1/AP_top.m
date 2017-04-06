clear
load('aggregation');
[n,~]=size(data);
%����ŷ�Ͼ�����෴����Ϊ���ƶȾ���
S = -squareform(pdist(data));
[row,col,v]=find(S);%�ҵ���Ϊ0��Ԫ��
%v��һ��n*(n-1)��Ԫ�صķ�0����������֮��ȡ��λ��
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
%��id_num������
for i=unique(idx)'
    %�ҵ���i���еĵ㣬��ż�Ϊii
    ii=find(idx==i);
    h=plot(data(ii,1),data(ii,2),'o');
    hold on;
    col=rand(1,3);
    set(h,'Color',col,'MarkerFaceColor',col,'MarkerSize',2);
    %�Ծ�������Ϊ��㣬�����Ϊ�յ㣬���������ڵ�����
     xi1=data(i,1)*ones(size(ii)); xi2=data(i,2)*ones(size(ii));
     line([data(ii,1),xi1],[data(ii,2),xi2],'Color',col,'LineWidth',2);
    
    h1=plot(data(i,1),data(i,2),'o');
    set(h1,'Color',col,'MarkerEdgeColor','k','MarkerSize',6);
    id_num = id_num + 1;
end;
title({'AP����';...
    ['��',num2str(id_num),'������']});