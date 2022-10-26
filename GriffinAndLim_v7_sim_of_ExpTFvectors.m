% Spectrogram function
 addpath( 'C:\Users\Lord Photon\Documents\MATLAB\library_repo\library' )
addpath( '/Users/ben/Documents/MATLAB/library_repo' )
%% Time frequency vectors definition

lent=2^13;                      % Signal length
tWind=4e-9;                   % Time window span


t=linspace(0,tWind,lent);
dt=t(2)-t(1);Fs=1/dt;f=linspace(-Fs/2,Fs/2,lent);df=(f(2)-f(1));
fG=f*10^-9;tps=t*10^12;%GHz
scale=1;

fMaxTLS=56.4e9;

%% SUT generation

fmax=31e9/2;%Fs/10;
% SUTf=superGauss(0,fmax,10,f,0).*(exp(1j*(tWind/4/(fmax*2*pi))*(2*pi*f).^2/2));
SUTf=superGauss(0,fmax,10,f,0).*(exp(1j*(240*22e-24)*(2*pi*f).^2/2));%+superGauss(0,fmax,10,f,0).*(exp(-1j*(120*22e-24/2)*(2*pi*f).^2/2));
% SUTf=superGauss(0,sutBW,10,f,0).*(exp(1j*(tWind/4/(sutBW*2*pi))*(2*pi*f).^2/2))+...
%     superGauss(0,sutBW,10,f,0).*(exp(-1j*(tWind/4/(sutBW*2*pi))*(2*pi*f).^2/2));
SUT=nifft(SUTf,Fs);

% SUT=ones(size(SUT));



%% window Setup

% Adjust these parameters as needed
winLen=2^8;
winLent=winLen*dt
winInc=winLen;%winLen-1;%/(2^2);
interpAmount_t=4; % For now, make this a power of 2 (or 1)!!
interpAmount_f=1; % For now, make this a power of 2 (or 1)!!


winSec=ones(1,winLen);
windowInds=(1:lent)-lent/2;
% win=superGauss(0,winLen/2,100,windowInds,0);
win=zeros(1,lent); win(round(lent/2)-winLen/2:round(lent/2)+winLen/2-1)=winSec;

win=circshift(win,lent/2);


% No need to change the ones below
nIncs=round(lent/winInc); % By making winInc a power of 2, we can make sure to have an integer number of windows.
windowCenters=(1:nIncs)*winInc;






%% Spectrogram Algorithm

fspgm_raw=f;
tspgm_raw=(1:nIncs)*winLent;


% Get spgm from windowIncrease Above
stftTLSRaw=get_TLSstft_fullSigLen(f,lent,winLen,winLent,fMaxTLS,dt,SUT);spgmTLSRaw=abs(stftTLSRaw).^2;
spgmTLSRaw=spgmTLSRaw-mean(spgmTLSRaw(1,:));spgmTLSRaw=abs(spgmTLSRaw/(max(max(spgmTLSRaw))));

[a,nZeros]=find(f>-fMaxTLS/2,1);
f_tls1=[f(1:nZeros-2),linspace(-fMaxTLS/2,fMaxTLS/2,winLen),-flip(f(1:nZeros-2))];
spgmRaw=interp2fun(tspgm_raw,f_tls1,[zeros(round((numel(f_tls1)-winLen)/2),nIncs);spgmTLSRaw;zeros(round((numel(f_tls1)-winLen)/2),nIncs)],tspgm_raw,f);
% stft(1:lent/2-winLen/2,:)=0;stft(lent/2+winLen/2:end,:)=0;
% stft=get_stft_fullSigLen(nIncs,windowCenters,lent,win,dt,SUT);



%%% PLaying around with COLA
% sutRecon=get_istft_fullSigLen(lent,windowCenters,win/(sum(win)/winInc),Fs,nIncs,stft);
% figure;
% for i=1:winInc
% plot(circshift(win,windowCenters(i)),'--')
% hold on
% end
% figure;imagesc(spgmRaw);
% figure;plot(sutRecon); hold on; plot(circshift(win,lent/2))
% Setup interpolation

%% Setup interpolation

nIncsInterp=nIncs*interpAmount_t;
windowCentersInterp=(1:nIncsInterp)*winInc/interpAmount_t;

tspgm=linspace(t(1),t(end),numel(spgmRaw(1,:))*interpAmount_t);
fspgm=linspace(f(1),f(end),numel(fspgm_raw)*interpAmount_f);%fspgm_raw;

%  spgm=interp2fun(tspgm_raw,fspgm,spgmRaw,tspgm,fspgm);
 spgm=imresize(spgmRaw,[lent,interpAmount_t*nIncs]);


winLenInterp=numel(fspgm)*interpAmount_f;
winInterp=interp1(linspace(0,1,lent),win,linspace(0,1,winLenInterp));
winIndsInterp=(1:numel(fspgm))-round(numel(fspgm)/2);
overlapAmount=interpAmount_t*(sum(win))/winInc;%numel(winIndsInterp)/(windowCentersInterp(2)-windowCentersInterp(1)); % This is the "overlapamount" AFTER interpolation
analysisWin=winInterp/overlapAmount; % Analysis window for the inverse spgm


%% Iterative Griffin and Lim algorithm

S0=sqrt(spgm);%.*exp(1j*rand(size(spgm))*2*pi);%.*(-1*(stft<0));%.*exp(1j*rand(size(spgm))*2*pi); % Seed stft

xt0=get_istft_fullSigLen(lent,windowCentersInterp,analysisWin,Fs,nIncsInterp,S0);
% xt0=get_istft_fullSigLen(lent,windowCentersInterp,analysisWin,Fs,nIncsInterp,stft);

maxIteration=50;

% Convergence criterion
di=zeros(1,maxIteration);
diC=di;
diR=di;
[xt,Si]=phaseRecovLoop(nIncsInterp,windowCentersInterp,lent,winInterp,winLen,t,dt,xt0,tspgm,fspgm,spgm,SUT,analysisWin,Fs,maxIteration);
%





%% Phase analysis 

sec1t=find(abs(xt)>max(abs(xt))*0.6,1):find(abs(xt)>max(abs(xt))*0.6,1,'last');
% xt=filtSG_tf(xt1,t,f,30e9/df,2,1);
% sec1t=554:1391;
% sec1f=391:465;

t1=t(sec1t);
xt1=xt(sec1t);
fitObj2=fit(t1',unwrap(angle(xt1))','poly2');
figure;plot(t1,unwrap(angle(xt1)));hold on; plot(t1,fitObj2(t1))
beta2_t=1/(2*fitObj2.p1)
beta2_t/22e-24

f1=linspace(-Fs/2,Fs/2,numel(t1));
xf1=(nfft(xt1));
sec1f=find(abs(xf1)>max(abs(xf1))*0.6,1):find(abs(xf1)>max(abs(xf1))*0.6,1,'last');
fitObj1=fit(2*pi*f1(sec1f)',unwrap(angle(xf1(sec1f)))','poly2');

figure;plot(2*pi*f1(sec1f),unwrap(angle(xf1(sec1f)))); hold on; plot(2*pi*f1(sec1f),fitObj1(2*pi*f1(sec1f)));
beta2=fitObj1.p1*2
beta2/22e-24
