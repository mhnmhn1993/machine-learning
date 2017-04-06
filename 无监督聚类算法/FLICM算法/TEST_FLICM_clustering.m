clear;
clc;
imIn = imread('ppp.png');
[imOut,C,U,iter] = My_FLICM( imIn,5);
figure
imshow(imIn, [0 255]), title('‘≠Õº');
figure
imshow(uint8(imOut),[0 255]), title('æ€¿‡Õº');