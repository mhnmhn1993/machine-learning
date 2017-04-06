clear;
x1 = [-5:0.1:5];
g=exp(-2*pi^2*(0.5)^2.*(x1-0).^2);
figure,plot(x1,g),xlabel('ÆµÂÊ'),ylabel('·ù¶È');

w = [0:0.1:100];
for i = 1:length(w)
  	x(i) = -(log10( w(i)/0.8 ))^2;
   x(i) = x(i)/(2*(log10(0.55))^2);
end
lgt = exp(x);  % Log Gabor Transfer function
figure,plot(w, lgt); % log Gabor transfer function viewed on linear frequency scale
figure,plot(log10(w), lgt);   % log Gabor transfer function viewed on logarithmic frequency scales
%lgtt = ifft(lgt);
%figure,plot(w,lgtt);
%even = lgtt.*cos(w);
%odd = lgtt.*sin(w);
%figure,plot(w,even,'r',w,odd,'b');
