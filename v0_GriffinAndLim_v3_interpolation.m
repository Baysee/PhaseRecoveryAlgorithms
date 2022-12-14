% Spectrogram function
 addpath( 'C:\Users\Lord Photon\Documents\MATLAB\library_repo\library' )
addpath( '/Users/ben/Documents/MATLAB/library_repo' )
%% Time frequency vectors definition

lent=2^13;                      % Signal length
tWind=800e-9;                   % Time window span


t=linspace(0,tWind,lent);
dt=t(2)-t(1);Fs=1/dt;f=linspace(-Fs/2,Fs/2,lent);df=(f(2)-f(1));
fG=f*10^-9;tps=t*10^12;%GHz
scale=1;



%% stft parameters

% Adjust these parameters as needed
winLen=2^8;
winInc=1;%winLen/2^4;%winLen-1;%/(2^2);
interpAmount_t=1; % For now, make this a power of 2 (or 1)!!
interpAmount_f=1; % For now, make this a power of 2 (or 1)!!

% win=hann(winLen+2).^3;win=win(2:end-1)';%ones(1,winLen);
win=ones(1,winLen);


% No need to change the ones below
nIncs=lent/winInc; % By making winInc a power of 2, we can make sure to have an integer number of windows.
windowCenters=(1:nIncs)*winInc;
winInds=(1:winLen)-winLen/2;


%% SUT generation

% SUT=exp(1j*pi*sin(2*pi/(40*winLen*dt)*t));



% % % % % % Linearly chirped signal
% fmax=Fs;
% fini=0; ffin=fmax/20;
% c=(ffin-fini)/tWind;
% SUT=exp(1j*pi*sin(2*pi*(c/2*t.^2+fini*t))).*superGauss(0,tWind/3,10,t,tWind/2);
% SUT=sin(2*pi*(c/2*t.^2+fini*t))+1;%.*superGauss(0,tWind/3,10,t,tWind/2);

fmax=Fs/30;
SUTf=superGauss(0,fmax,10,f,0).*(exp(1j*(tWind/4/(fmax*2*pi))*(2*pi*f).^2/2));
% SUTf=superGauss(0,sutBW,10,f,0).*(exp(1j*(tWind/4/(sutBW*2*pi))*(2*pi*f).^2/2))+...
%     superGauss(0,sutBW,10,f,0).*(exp(-1j*(tWind/4/(sutBW*2*pi))*(2*pi*f).^2/2));
SUT=nifft(SUTf,Fs);
% SUT=superGauss(0,tWind/20,1,t,2*tWind/3)+superGauss(0,tWind/10,2,t,tWind/3);

% % % % % ASE noise
% fmax=Fs/20;
%
%  SUT=(randn(1,lent)+1j*rand(1,lent));
%
%  SUTf=nfft(SUT,dt);
%  fShape=superGauss(0,fmax,100,f,0);
%  SUTf=SUTf.*fShape;
%  tShape=superGauss(0,tWind/3,4,t,tWind/2);
% SUT=nifft(SUTf,Fs).*tShape;
% SUT=SUT-mean(SUT);


%% Spectrogram Algorithm

% Get spgm from windowIncrease Above
stft=get_stft_winInds(nIncs,winInds,windowCenters,lent,win,dt,winLen,SUT);
spgmRaw=abs(stft).^2;
fspgm_raw=((1:winLen)-winLen/2)/winLen*Fs;
tspgm_raw=linspace(t(1),t(end),numel(stft(1,:)));

figure;imagesc(spgmRaw);

% Setup interpolation
nIncsInterp=nIncs*interpAmount_t;
windowCentersInterp=(1:nIncsInterp)*winInc/interpAmount_t;

tspgm=linspace(t(1),t(end),numel(stft(1,:))*interpAmount_t);
fspgm=linspace(f(1),f(end),numel(fspgm_raw)*interpAmount_f);%fspgm_raw;

[tspgm_rawM,fspgm_rawM]=meshgrid(tspgm_raw,fspgm);
[tspgmM,fspgmM]=meshgrid(tspgm,fspgm);


% spgm=interp2(tspgm_rawM,fspgm_rawM,spgmRaw,tspgmM,fspgmM,'spline');
% spgm=griddata(tspgm_rawM,fspgm_rawM,spgmRaw,tspgmM,fspgmM,'natural');


% spgm1=griddedInterpolant(tspgm_rawM,fspgm_rawM,spgmRaw,tspgmM,fspgmM,'linear');
spgmInterpolant=griddedInterpolant({tspgm_raw,fspgm_raw},spgmRaw','nearest');

spgmInterp=spgmInterpolant({tspgm,fspgm})';
% spgm=griddata(tspgm_rawM,fspgm_rawM,ab,tspgmM,fspgmM);

spgm=spgmRaw;%imgaussfilt(spgmInterp,interpAmount_t/2);
figure;subplot(2,1,1)
imagesc(spgmRaw)
subplot(2,1,2)
imagesc(spgm)

%t(windowCenters);

% 
% %% Spectrogram Plot
% h0=figure;
% h0.Position=[-1402 147 669 830];
% xlimsZoom=t(round(end/2))+winLen*dt*[-0.5 0.5];
% FS=16;
% 
% subplot(2,2,1)
% plot(t,abs(SUT).^2);
% ylabel('Intensity');
% yyaxis right
% plot(t,angle(SUT));
% xlabel('Time (a.u.)'); ylabel('phase (rad)')
% set(gca,'FontSize',FS)
% 
% subplot(2,2,2)
% plot(t,abs(SUT).^2);
% ylabel('Intensity');
% yyaxis right
% plot(t,angle(SUT));
% xlabel('Time (a.u.)'); ylabel('phase (rad)')
% xlim(xlimsZoom);
% set(gca,'FontSize',FS)
% 
% subplot(2,2,3:4)
% imagesc(tspgm,fspgm,spgm);
% ylim([-fmax fmax])
% xlabel('Time'); ylabel('Frequency');
% set(gca,'FontSize',FS)




%% Iterative Griffin and Lim algorithm


% ispgm=get_istft_winInds(lent,winInds,windowCenters,analysisWin,Fs,nIncs,S);
%
% figure;plot(real(ispgm)); hold on; %plot(imag(ispgm)); %plot(abs(ispgm));
% plot(SUT)
% legend('Real ispgm','real SUT')
winLenInterp=numel(fspgm);
winInterp=interp1(linspace(0,1,winLen),win,linspace(0,1,winLenInterp));
winIndsInterp=(1:numel(fspgm))-round(numel(fspgm)/2);
overlapAmount=numel(winIndsInterp)/(windowCentersInterp(2)-windowCentersInterp(1)); % This is the "overlapamount" AFTER interpolation
analysisWin=winInterp/(overlapAmount); % Analysis window for the inverse spgm

S0=sqrt(spgm);%.*(-1*(stft<0));%.*exp(1j*rand(size(spgm))*2*pi); % Seed stft
S0Mag=abs(sqrt(spgm)); % S0 has the correct amplitude, but not the correct phase
xt=get_istft_winInds(lent,winIndsInterp,windowCentersInterp,analysisWin,Fs,nIncsInterp,S0);
xt0=xt; % This is the first initial guess

maxIteration=200;
i=1;

% Convergence criterion
di=zeros(1,maxIteration);
diC=di;
diR=di;

h1=figure;
h1.Position=[55 117 990 861];%[-1528 112 821 865];
xlimsZoom=t(round(lent/2))+winLen*dt*[-2 2];

while i<maxIteration+1
    
    Si=get_stft_winInds(nIncsInterp,winIndsInterp,windowCentersInterp,lent,winInterp,dt,winLenInterp,xt);         % Get spgm of present signal guess xt
    Sip1=S0Mag.*Si./abs(Si);%.*exp(1j*angle(Si));%.*Si./abs(Si);                           % Enforce magnitude along with calculated phase from Si
    Sip1(isnan(Sip1))=0;
    xt=get_istft_winInds(lent,winIndsInterp,windowCentersInterp,analysisWin,Fs,nIncsInterp,Sip1);
    % di(i)=sqrt(sum(sum(abs(abs(Si)-abs(Sip1)).^2))
    % di(i)=norm(abs(abs(Si)-abs(Sip1)).^2)/norm(abs(Si).^2)
%     di(i)=(norm(abs(Si)-abs(Sip1))).^2;%)/norm(abs(Si).^2)
    diC(i)=sqrt(    sum(sum( abs(  abs(S0) - abs(Si)  ).^2))  /    sum(sum(  abs(S0).^2 )) );%)/norm(abs(Si).^2)
%     diC2(i)=sqrt(    sum(sum( abs(  sqrt(abs(S0).^2) - sqrt(abs(Si).^2)  ).^2 ))  /    sum(sum(  abs(S0).^2 )) );%)/norm(abs(Si).^2)
%     diC(i)=sqrt(    norm( abs(  sqrt(abs(S0).^2) - sqrt(abs(Si).^2)  ).^2, 'fro')  /    norm(  abs(S0).^2 , 'fro') );%)/norm(abs(Si).^2)
    diR(i)=sqrt(sum(abs(SUT-xt).^2)/sum(abs(SUT).^2));%)/norm(abs(Si).^2)
    xc=max(abs(xcorr(SUT,xt,50)))/(sqrt(sum(abs(SUT).^2)*sum(abs(xt).^2)));
    diR(i)=(1-xc);%)/norm(abs(Si).^2)
    
    % plot(t,real(xt)); drawnow();
    i=i+1;
    
    if mod(i,4)==1
    updatePlot(h1,xlimsZoom,t,tspgm,fspgm,xt0,SUT,xt,spgm,Si,diC,diR)
    end
end


    updatePlot(h1,xlimsZoom,t,tspgm,fspgm,xt0,SUT,xt,spgm,Si,diC,diR)











function updatePlot(h1,xlimsZoom,t,tspgm,fspgm,xt0,SUT,xt,spgm,Si,diC,diR)
figure(h1);
FS=16;
subplot(3,2,1)
imagesc(tspgm,fspgm,spgm)
title('SUT spgm')
subplot(3,2,2)
imagesc(tspgm,fspgm,abs(Si).^2)
title('current Iteration spgm')
subplot(3,2,3)
plotIniFin(t,xt0,xt,SUT,FS)
subplot(3,2,4)
plotIniFin(t,xt0,xt,SUT,FS)
ylims=ylim(); xlim(xlimsZoom); ylim(ylims)

subplot(3,2,5:6)
yyaxis left
plot(20*log10(diC))
ylabel('Normalized inconsistency (dB)')
yyaxis right
% plot(20*log10(diC)); hold on
plot(10*log10(diR))
ylabel('10*log10(1-[corr. coeff])')
legend('Convergence Curve (GL)','spectral convergence');%,'reconstruction error')
xlabel('Iteration'); 
set(gca,'FontSize',FS)
drawnow();
%
end



