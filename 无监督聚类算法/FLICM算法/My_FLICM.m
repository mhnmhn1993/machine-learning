% ���룺
% imgIn������Ҷ�ͼ��
% cNum����������Ĭ��Ϊ2
% m:ģ��ϵ����Ĭ��Ϊ2
% winSize������Ĵ��ڴ�С��������Ĭ��Ϊ3
% maxIter��������������Ĭ��Ϊ100
% thrE�����ȣ�Ĭ��Ϊ1.0e-2

% �����
% imOut�� ���ͼ��
% C����������
% U�������Ⱦ���
% iter����������



function [imOut,C,U,iter] = My_FLICM( imIn, cNum, m, winSize, maxIter, thrE )
if nargin < 6
    thrE = 1.0e-2;
end

if nargin < 5
    maxIter = 100;
end

if nargin < 4
    winSize = 3;
end

if nargin < 3
    m = 2;
end

if nargin < 2
    cNum = 2;
end


if( size(imIn,3)>1 )
    imIn = rgb2gray( imIn );
end
image = double( imIn );

% �����Ⱦ���ĳ�ʼ��
iter = 0;
[H,W,p] = size( image );
N = H*W;
U = rand( H, W, cNum-1 )*(1/cNum);
U(:,:,cNum) = 1 - sum(U,3);

X = reshape(image, H*W, p);
P = reshape(U, H*W, cNum);

hh = round((winSize-1)/2);
ww = round((winSize-1)/2);
Uplus = ones(H+2*hh, W+2*ww, cNum);
DDplus = zeros(H+2*hh, W+2*ww, cNum);
J_prev = inf; 
J = [];
while true
   iter = iter + 1; 
   if iter>maxIter
       break;
   end
   
   t = P.^m;
   C = (t'*X)./(sum(t)'*ones(1, p));
   dist = sum(X.*X, 2)*ones(1, cNum) + (sum(C.*C, 2)*ones(1, N))'-2*X*C';
      
   DDplus(1+hh:H+hh, 1+ww:W+ww, :) = reshape(dist, H, W, cNum);
   Uplus(1+hh:H+hh, 1+ww:W+ww, :) = U;
   
   GG = zeros(H, W, cNum);
   for ii = -hh:hh
       for jj = -ww:ww
           if ii~=0&&jj~=0
               GG = GG + 1/(1+ii*ii+jj*jj)*(1-Uplus(1+hh+ii:H+hh+ii, 1+ww+jj:W+ww+jj, :)).^m.*DDplus(1+hh+ii:H+hh+ii, 1+ww+jj:W+ww+jj, :);
           end
       end
   end
   
   G = reshape(GG, H*W, cNum);
   rd = dist + G;
   t2 = (1./rd).^(1/(m-1));
   P = t2./(sum(t2, 2)*ones(1, cNum));
   U = reshape(P, H, W, cNum);
   J_cur = sum(sum((P.^m).*dist+G))/N;
   J = [J J_cur];
   if norm(J_cur-J_prev, 'fro') < thrE
       break;
   end
   J_prev = J_cur;
end

[~, label] = max(P, [], 2);

imOut = reshape(C(label, :), H, W, p);











