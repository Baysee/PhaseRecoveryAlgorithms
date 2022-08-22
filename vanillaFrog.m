% % % Vanilla Frog Algorithm


% Spectrogram function
addpath( '/Users/ben/Documents/MATLAB/library' )
%% Time frequency vectors definition

lent=2^16;                      % Signal length
tWind=10e-9;                   % Time window span

% t=linspace(0,tWind,lent);
t=linspace(-tWind/2,tWind/2,lent);
dt=t(2)-t(1);Fs=1/dt;f=linspace(-Fs/2,Fs/2,lent);df=(f(2)-f(1));
fGHz=f*10^-9;tps=t*10^12; tns=t*1e9;
fG=1e-9; tn=1e9;
scale=1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% temporal phase signal definiton %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% These settings are used for both the T-TAI and the time-lens spectrogram
trTry=0.95e-9;%9.0761e-9;        % Time-lens apperture
eta=20;                % number of analysis points
fmax=Fs/2;         % full bandwidth that can be processed with the spectrogram

ntr=round(trTry/dt);    % tr is adjusted to be an integer amount of points
tr=ntr*dt;
nSamps=ceil(lent/ntr);
tSingleSamp=dt*(1:ntr);
lensInds=0:ntr:lent;




%%%%%%%%%%%%%%%%%%
%% Generate SUT %%
%%%%%%%%%%%%%%%%%%
SUT=ones(1,lent);



% Linearly chirped signal
fini=-fmax/2;0; ffin=fmax/2;
c=(ffin-fini)/tWind;
% SUT=sin(2*pi*(c/2*t.^2+fini*t));

sutWind=tWind/16; sutCenter=(t(end)+t(1))/2;
SUT=ones(1,lent).*exp(1j*2*pi*(c/2*t.^2+fini*t));
SUT=superGauss(0,sutWind,10,t,sutCenter);

% Help setup frequency range by carving slightly bigger than sut, then
% multiply by gaussain frequency chirp
% superGauss(0,sutWind*1.5,30,t,t(end)/2).*

% instantaneous frequency parameters
fEdge=fmax/4
fMid=-fmax/6
centerMid=sutCenter+tWind/50;
instF_sigma=sutWind/1.8;
SUTinstFIni=superGauss(0,instF_sigma,1,t,centerMid);
riseInd=find(SUT>0.5*max(SUT),1); fallInd=find(SUT>0.5*max(SUT),1,'last'); 
% rangeY=max(SUTinstFIni)-min([SUTinstFIni(riseInd), SUTinstFIni(fallInd)]);
% SUTinstF=SUTinstFIni*(abs(fEdge-fMid))/rangeY+min(fEdge,fMid);
extraBitFactor=20;
SUTregion=riseInd-round(sutWind/dt/extraBitFactor)+1:-1+fallInd+round(sutWind/dt/extraBitFactor);
extraBitFactorPlus=4;US=20;
SUTregP_US=riseInd-round(sutWind/dt/extraBitFactorPlus)+1:US:-1+fallInd+round(sutWind/dt/extraBitFactorPlus);

SUTinstFIni([1:SUTregion(1)-1,SUTregion(end)+1:end])=0;
SUTinstFIni(SUTregion)=SUTinstFIni(SUTregion)-min(SUTinstFIni(SUTregion)); SUTinstFIni=SUTinstFIni/max(SUTinstFIni);
SUTinstF=SUTinstFIni*(abs(fEdge-fMid))+min(fEdge,fMid);
SUTphase=cumtrapz(t,SUTinstF);

SUT=SUT.*exp(1j*2*pi*SUTphase);



SUT_f=nfft(SUT,dt);


figure;
subplot(3,1,1)
plot(t,abs(SUT).^2); hold on ; plot(t,SUTinstFIni); 
yyaxis right
plot(t,unwrap(angle(SUT)))
xlabel('Time')
legend('sut power','sutinstFini','anle SUT')
subplot(3,1,2)
plot(SUTinstF); 
yyaxis right; plot(SUTphase)
legend('SUTinstF','SUTphase')
xlabel('Time')

subplot(3,1,3)
plot(f,abs(nfft(SUT)))
xlabel('Freq')
xlim([-fmax/2 fmax/2])


%%%%%%%%%%%%%%%%
%% PG - FROG %%
%%%%%%%%%%%%%%%%
nDelays=1000;
tauInc_i=round(tWind/nDelays/dt);
tau=tauInc_i*dt;
tauDelays_i=circshift([1:nDelays]*tauInc_i,round(nDelays/2));
tauDelays=([1:nDelays]*tau-tWind/2);

shiftedSuts_i=[1:lent]'+tauDelays_i;
shiftedSuts_i=mod(shiftedSuts_i,lent+1); shiftedSuts_i(shiftedSuts_i==0)=1;
Psut=abs(SUT).^2;
shiftedSuts=Psut(shiftedSuts_i);
Esig=SUT'.*shiftedSuts;

EpgFrog=fftshift(fft(ifftshift(Esig)));
IpgFrog=abs(EpgFrog).^2;



%%% Reshape output signal to get 2D representation
spgm=EpgFrog;
f_spec=f;
t_spec=tauDelays;



%% Plotting section

ylims=([-0.1 1.1]);
ylimsNegPos=([-1.1 1.1]);
numtrWindLims=20;
tlims='auto';%(t(end)/2+[-numtrWindLims*tr numtrWindLims*tr])*1e9;
% tlims=[t(1) t(end)];

tps=1e12; tn=1e9; tus=1e6;
lbps='Time (ps)'; lbns='Time (ns)'; lbus='Time (us)';
fG=1e-9; lbG='Frequency (GHz)';
figure('Renderer', 'painters', 'Position', [50 50 1200 700])


subplot(4,4,1:4)
plot(t*tn,abs(SUT).^2/(max(abs(SUT).^2)))
yyaxis right
plot(t*tn,angle(SUT));%/(max(angle(SUT))))

xlim(tlims);%ylim(ylims);
xlabel(lbns)
legend('sampled SUT','Dispersed wvf','Time Lens')
ylabel('Phase (rad)')



subplot(4,4,5:8)

plot(fGHz,abs(SUT_f).^2/(max(abs(SUT).^2)))
yyaxis right
plot(fGHz,angle(SUT_f));%/(max(angle(SUT))))
% 
% plot(t*tn,abs(sampSUTdisp).^2/(max(abs(sampSUTdisp).^2)))
% yyaxis right
% plot(t*tn,sampSig)

ylabel('Intensity')
% yyaxis right
% plot(t*tns,sampSig)
xlim(tlims);%ylim(ylimsNegPos);
xlabel(lbns)
legend('abs','angle','spectrogram')
ylabel('Phase (rad)')

subplot(4,4,9:16)
% hcb=colorbar();
% xlim(tlims)
imagesc(t_spec*tn,f_spec*fG,abs(spgm).^2)
ylabel('Power (dB)')
xlabel('Time (ns)')



