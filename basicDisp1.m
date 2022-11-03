lent=2^12;                      % Signal length
tWind=800e-9;                   % Time window span

lent=2^18;
tWind=100e-9;
t=linspace(0,tWind,lent)-tWind/2;
dt=t(2)-t(1);Fs=1/dt;f=linspace(-Fs/2,Fs/2,lent);df=(f(2)-f(1));
fG=f*10^-9;tps=t*10^12;%GHz
scale=1;

beta2Km=-21.68e-24;
beta3Km=0.12661e-36;
nKm=210;

y=superGauss(0,10e-12/(2*sqrt(2*log(2))),1,t,0);

yDisp=nifft(nfft(y).*exp(1j*beta2Km/2*nKm*(2*pi*f).^2).*exp(1j*beta3Km*nKm/6*(2*pi*f).^3));
figure;
plot(t,y)
yyaxis right
plot(t, abs(yDisp)) 

hold on
plot(t, real(yDisp))
