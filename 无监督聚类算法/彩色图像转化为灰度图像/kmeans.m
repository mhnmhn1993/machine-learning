%һ�ֻ���K��ֵ����Ȼͼ����෽��
clc
clear 
pic = imread('ppp.png');  % RGB  ��ά����
figure;
imshow(pic);
k = 10;  
pic = double(pic);       %double���κ���������ת����˫������ֵ
copy = pic;         
ima = pic(:);            % ���һ��
s = length(ima);
data1 = pic(:,:,1);
data2 = pic(:,:,2);
data3 = pic(:,:,3);

%xlswrite('D:\write2Excel.xls',data1,'sheet1');%���ݴ洢��D�̸�Ŀ¼��
%xlswrite('D:\write2Excel.xls',data2,'sheet2')
%xlswrite('D:\write2Excel.xls',data3,'sheet3')

%�����Ҷ�ֱ��ͼ
%�Ҷ�ֱ��ͼ�ǻҶȼ��ĺ���������ʾͼ���о���ĳ�ֻҶȼ������صĸ�������ӳ��ͼ����ĳ�ֻҶȳ��ֵ�Ƶ�ʡ�
m = max(ima);            % m = max(ima)+1;
h = zeros(1,m);          % һ����1-m���Ҷ�
hc = zeros(1,m);                                                                                                                           
for i = 1:s
  if(ima(i)>0) 
      h(ima(i )) = h(ima(i)) + 1;
  end
end
ind = find(h);           % �ҵ����лҶȡ�h������ÿ���Ҷȵĸ���
hl = length(ind);        % �Ҷȳ���

mu = (1:k)*m/(k);           %��ʼ����������
while 1
  oldmu = mu;    
  for i = 1:hl      
      c = abs(ind(i)-mu);   %  ����
      cc = find(c==min(c)); %����̵ľ������ھ���λ�ü�����������
      hc(ind(i)) = cc(1);   %  hc������ĳһ�Ҷ����ڵĴ��� 
  end  
  for i = 1:k               %  ���¾�������
      a = find(hc==i);      %  ��hc�������ҵ����ڵ�i��ص�����     
      mu(i) = sum(a.*h(a))/sum(h(a));   %�Ҷȶ���   
  end  
  if(mu == oldmu)           %  �������Ĳ��ٷ����ı�
      break;
  end
end

S = size(copy);
mask = zeros(S);
for i = 1:S(1)       %���»���ͼ��
    for j = 1:S(2)     
        copy(i,j,1) - mu;
        c = abs(copy(i,j,1)-mu);
        a = find(c==min(c));  
        mask(i,j,1) = mu(a);
    end
end
mask  = uint8(mask); %8λ���޷�����������,ȡֵ��Χ��0��255
figure;
imshow(mask(:,:,1));
figure;
imshow(mask(:,:,2));
figure;
imshow(mask(:,:,3));




