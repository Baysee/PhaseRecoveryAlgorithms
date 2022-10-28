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



%% SUT generation

fmax=180e9/2;%Fs/10;
% SUTf=superGauss(0,fmax,10,f,0).*(exp(1j*(tWind/4/(fmax*2*pi))*(2*pi*f).^2/2));
SUTf=superGauss(0,fmax,10,f,0).*(exp(1j*(240*22e-24/2)*(2*pi*f).^2/2));%+superGauss(0,fmax,10,f,0).*(exp(-1j*(120*22e-24/2)*(2*pi*f).^2/2));
% SUTf=superGauss(0,sutBW,10,f,0).*(exp(1j*(tWind/4/(sutBW*2*pi))*(2*pi*f).^2/2))+...
%     superGauss(0,sutBW,10,f,0).*(exp(-1j*(tWind/4/(sutBW*2*pi))*(2*pi*f).^2/2));
SUT=nifft(SUTf,Fs);

% SUT=ones(size(SUT));



%% window Setup


% Adjust these parameters as needed
fMaxTLS=200e9;
winLen=2^7;
winLent=winLen*dt
winInc=winLen;%winLen-1;%/(2^2);
interpAmount_t=1; % For now, make this a power of 2 (or 1)!!
interpAmount_f=1; % For now, make this a power of 2 (or 1)!!

winSec=ones(1,winLen);
windowInds=(1:lent)-lent/2;
% win=superGauss(0,winLen/2,100,windowInds,0);
win=zeros(1,lent); win(round(lent/2)-winLen/2:round(lent/2)+winLen/2-1)=winSec;
win=circshift(win,lent/2);

% No need to change the ones below
nIncs=round(lent/winInc); % By making winInc a power of 2, we can make sure to have an integer number of windows.
windowCenters=(1:nIncs)*winInc;


CL=2*pi*fMaxTLS/winLent;
phi2=1/CL;

phit=repmat( ((1:winLen)/winLen*winLent-winLent/2).^2*CL/2, [1,nIncs]) ;
phiw=phi2/2*(2*pi*f).^2;




%% Spectrogram Algorithm

% Get spgm from windowIncrease Above
[stft]=stft_TLS(SUT,winLen,phit,phiw);
% [xt]=istft_TLS(stft,winLen,phit,phiw);
% [stft2]=stft_TLS(xt,winLen,phit,phiw);

% stft=get_stft_fullSigLen(nIncs,windowCenters,lent,win,dt,SUT);
spgmRaw=abs(stft).^2;
% fspgm_raw=f;
% tspgm_raw=linspace(t(1),t(end),numel(stft(1,:)));


tspgm=linspace(t(1),t(end),numel(stft(1,:)));
% fspgm=linspace(f(1),f(end),numel(stft(:,1)));%fspgm_raw;
fspgm=linspace(0,fMaxTLS,numel(stft(:,1)))-fMaxTLS/2;%fspgm_raw;

 spgm=spgmRaw;


%% Iterative Griffin and Lim algorithm

S0=sqrt(spgm).*exp(1j*rand(size(spgm))*2*pi);%.*(-1*(stft<0));%.*exp(1j*rand(size(spgm))*2*pi); % Seed stft

xt0=istft_TLS(S0,winLen,phit,phiw);;%get_istft_fullSigLen(lent,windowCenters,analysisWin,Fs,nIncs,S0);
xt_sut=istft_TLS(stft,winLen,phit,phiw);;%get_istft_fullSigLen(lent,windowCenters,analysisWin,Fs,nIncs,S0);
% xt0=get_istft_fullSigLen(lent,windowCentersInterp,analysisWin,Fs,nIncsInterp,stft);

maxIteration=800;

% Convergence criterion
di=zeros(1,maxIteration);
diC=di;
diR=di;
% [xt,Si]=phaseRecovLoop(nIncs,windowCenters,lent,win,winLen,t,dt,xt0,tspgm,fspgm,spgm,SUT,analysisWin,Fs,maxIteration);
[xt,Si]=phaseRecovLoopTLS(nIncs,windowCenters,lent,winLen,t,f,dt,xt0,tspgm,fspgm,spgm,SUT,Fs,maxIteration,phit,phiw);
%


