% Spectrogram function
 addpath( 'C:\Users\Lord Photon\Documents\MATLAB\library_repo\library' )
%% Time frequency vectors definition

lent=2^12;                      % Signal length
tWind=800e-9;                   % Time window span


t=linspace(0,tWind,lent);
dt=t(2)-t(1);Fs=1/dt;f=linspace(-Fs/2,Fs/2,lent);df=(f(2)-f(1));
fG=f*10^-9;tps=t*10^12;%GHz
scale=1;



%% spgm parameters

% Adjust these parameters as needed
winLen=2^7
winInc=winLen/2^5;%winLen-1;%/(2^2);

% win=hann(winLen+2).';win=win(2:end-1); %%% Window shape
win=ones(1,winLen);


% No need to change the ones below
nIncs=lent/winInc; % By making winInc a power of 2, we can make sure to have an integer number of windows.
windowCenters=(1:nIncs)*winInc;
winInds=(1:winLen)-winLen/2;


%% SUT generation

% SUT=exp(1j*pi*sin(2*pi/(40*winLen*dt)*t));



% % % % % % Linearly chirped signal
fmax=Fs;
fini=0; ffin=fmax/10;
c=(ffin-fini)/tWind;
% SUT=exp(1j*pi*sin(2*pi*(c/2*t.^2+fini*t))).*superGauss(0,tWind/3,10,t,tWind/2);
% SUT=sin(2*pi*(c/2*t.^2+fini*t))+1;%.*superGauss(0,tWind/3,10,t,tWind/2);

sutBW=Fs/3;
SUTf=superGauss(0,sutBW,10,f,0).*(exp(1j*(tWind/4/(sutBW*2*pi))*(2*pi*f).^2/2));
SUT=nifft(SUTf,Fs)
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


spgm=get_spgm(nIncs,winInds,windowCenters,lent,win,dt,winLen,SUT);


stft=abs(spgm).^2;
tstft=t(windowCenters);
fstft=((1:winLen)-winLen/2)/winLen*Fs;


%% Spectrogram Plot
h0=figure;
h0.Position=[-1402 147 669 830];
xlimsZoom=t(round(end/2))+winLen*dt*[-0.5 0.5];
FS=16;

subplot(2,2,1)
plot(t,abs(SUT).^2);
ylabel('Intensity');
yyaxis right
plot(t,angle(SUT));
xlabel('Time (a.u.)'); ylabel('phase (rad)')
set(gca,'FontSize',FS)

subplot(2,2,2)
plot(t,abs(SUT).^2);
ylabel('Intensity');
yyaxis right
plot(t,angle(SUT));
xlabel('Time (a.u.)'); ylabel('phase (rad)')
xlim(xlimsZoom);
set(gca,'FontSize',FS)

subplot(2,2,3:4)
imagesc(tstft,fstft,stft);
ylim([-fmax fmax]/2)
xlabel('Time'); ylabel('Frequency');
set(gca,'FontSize',FS)




%% Iterative Griffin and Lim algorithm


% istft=get_ispgm(lent,winInds,windowCenters,analysisWin,Fs,nIncs,S);
%
% figure;plot(real(istft)); hold on; %plot(imag(istft)); %plot(abs(istft));
% plot(SUT)
% legend('Real istft','real SUT')


overlapAmount=numel(winInds)/(windowCenters(2)-windowCenters(1));
analysisWin=win/(overlapAmount); % Analysis window for the inverse STFT
fullWin=zeros(size(SUT)); fullWin(1:winLen)=1; fullWin=fullWin';

S0=sqrt(stft);%.*(-1*(spgm<0));%.*exp(1j*rand(size(stft))*2*pi); % Seed spgm
S0Mag=abs(sqrt(stft)); % S0 has the correct amplitude, but not the correct phase
xt=get_ispgm(lent,winInds,windowCenters,analysisWin,Fs,nIncs,S0);
xt0=xt; % This is the first initial guess

maxIteration=1000;
i=1;

% Convergence criterion
di=zeros(1,maxIteration);
diC=di;
diR=di;

h1=figure;
h1.Position=[55 117 990 861];%[-1528 112 821 865];
xlimsZoom=t(round(end/2))+winLen*dt*[-2 2];

while i<maxIteration+1
    
    Si=get_spgm(nIncs,winInds,windowCenters,lent,win,dt,winLen,xt);         % Get STFT of present signal guess xt
    Sip1=S0Mag.*sign(Si);%Si./abs(Si);%.*exp(1j*angle(Si));%.*Si./abs(Si);                           % Enforce magnitude along with calculated phase from Si
    Sip1(isnan(Sip1))=0;
    xt=get_ispgm(lent,winInds,windowCenters,analysisWin,Fs,nIncs,Sip1);
    % di(i)=sqrt(sum(sum(abs(abs(Si)-abs(Sip1)).^2))
    % di(i)=norm(abs(abs(Si)-abs(Sip1)).^2)/norm(abs(Si).^2)
%     di(i)=(norm(abs(Si)-abs(Sip1))).^2;%)/norm(abs(Si).^2)
    diC(i)=sqrt(    sum(sum( abs(  sqrt(abs(S0).^2) - sqrt(abs(Si).^2)  ).^2 ))  /    sum(sum(  abs(S0).^2 )) );%)/norm(abs(Si).^2)
%     diC(i)=sqrt(    norm( abs(  sqrt(abs(S0).^2) - sqrt(abs(Si).^2)  ).^2, 'fro')  /    norm(  abs(S0).^2 , 'fro') );%)/norm(abs(Si).^2)
    diR(i)=sqrt(sum(abs(SUT-xt).^2)/sum(abs(SUT).^2));%)/norm(abs(Si).^2)
    % plot(t,real(xt)); drawnow();
    i=i+1;
    updatePlot(h1,xlimsZoom,t,tstft,fstft,xt0,SUT,xt,stft,Si,diC,diR)
end



function updatePlot(h1,xlimsZoom,t,tstft,fstft,xt0,SUT,xt,stft,Si,diC,diR)
figure(h1);
FS=16;
subplot(3,2,1)
imagesc(tstft,fstft,stft)
title('SUT STFT')
subplot(3,2,2)
imagesc(tstft,fstft,abs(Si).^2)
title('current Iteration STFT')
subplot(3,2,3)
plotIniFin(t,xt0,xt,SUT,FS)
subplot(3,2,4)
plotIniFin(t,xt0,xt,SUT,FS)
ylims=ylim(); xlim(xlimsZoom); ylim(ylims)

subplot(3,2,5:6)
yyaxis left
plot(20*log10(diC))
yyaxis right
% plot(20*log10(diC)); hold on
plot(20*log10(diR))
legend('Convergence Curve (GL)','spectral convergence');%,'reconstruction error')
xlabel('Iteration'); ylabel('Normalized inconsistency (dB)')
set(gca,'FontSize',FS)
drawnow();
%
end
%
% 
% overlapAmount=numel(winInds)/(windowCenters(2)-windowCenters(1));
% checkCola_2D=zeros(winLen,overlapAmount);
% for i=1:overlapAmount
% %     checkCola_2D(:,i)=circshift(win,(i-1)*winInc);
% checkCola_2D(:,i)=analysisWin.*circshift(win,(i-1)*winInc);
% end
% figure;plot(sum(checkCola_2D.'))



function istft=get_ispgm(lent,winInds,windowCenters,analysisWin,Fs,nIncs,S)
% Reconstruct temporal waveform from stft by the overlap and add technique

istft=zeros(1,lent); % Initialize variable

for i=1:nIncs
    sutInds=mod(winInds+windowCenters(i),lent)+1;
    istft(sutInds)=istft(sutInds)+analysisWin.*nifft(S(:,i),Fs).';
end
% istft=istft;
end



function spgm=get_spgm(nIncs,winInds,windowCenters,lent,win,dt,winLen,SUT)

spgm=zeros(winLen,nIncs);

for i=1:nIncs
    sutInds=mod(winInds+windowCenters(i),lent)+1;
    spgm(:,i)=nfft(SUT(sutInds).*win,dt);
end

end



