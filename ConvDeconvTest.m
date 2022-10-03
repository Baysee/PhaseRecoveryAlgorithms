% % convDeconv Test\

lent=2^10;

x1_=zeros(1,2^4);
x1_(1)=1;
x1=repmat(x1_,[1,lent/numel(x1_)]);

x2=singleGauss(2^2,lent/2,1:lent,0);
x2=circshift(x2,lent/2)

x3=conv(x1,x2);
[x4,b]=deconv(x3,x2);

figure;
subplot(3,2,1)
plot(x1); 
subplot(3,2,2)
plot(x2)
subplot(3,2,3:4)
plot(x3)
subplot(3,2,5)
plot(x4)
subplot(3,2,6)
plot(b)