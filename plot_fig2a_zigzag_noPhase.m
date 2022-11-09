% analyzeMenlo_simple

addpath('/Volumes/Extreme SSD/MATLAB_DataAnalysis/analysisDocs')
addpath( '/Users/ben/Documents/MATLAB/library' )
addpath( '/Users/ben/Documents/MATLAB/library/export_fig' )
addpath('/Users/ben/Documents/MATLAB/timeFrequencyAnalysis/experimentalData')
generalFolder_scope='spectrogram_20211109/';
apexFolder='/ApexControl/';

%% Bandwidth demonstration



%% ZIG ZAG
fnFolder='20211119/Menlo/200ps/';
datDir=[generalFolder_scope fnFolder];
% Relevant dates with menlo data: 14, 16,17, 18, 19

% % Run 1
% basefn_scope='_menlo56GHzBW_240&600km_ZigZag_27km200ps10h44m.mat';
% basefn_OSA='menlo56GHzBW_240&600km_ZigZag_2pm_';
% basefn_OSAtime='10h49m.mat';
% intOfInterest=3.6e4:4.4e4;%2e4:4e4;
% % Run 1
% basefn_scope='_menlo56GHzBW_240&600km_ZigZag_27km200ps10h43m.mat';
% basefn_OSA='menlo56GHzBW_240&600km_ZigZag_2pm_';
% basefn_OSAtime='10h49m.mat';
% intOfInterest=3.6e4:4.4e4;%2e4:4e4;
% % 


% Run 1
basefn_scope='_menlo40GHzBW_240&600km_ZigZag_27km200ps09h47m.mat';
basefn_OSA='menlo40GHzBW_240&600km_ZigZag_2_';
basefn_OSAtime='10h01m.mat';
intOfInterest=6e3:9e3;%2e4:4e4;



sigType_scope={'spectrogram', 'pmoff','SUT','SUT_TLS','SUT_trigger','TLS','trigger'};
sigType_OSA={'SUT','pmoff','pmon'};



% load([generalFolder_scope, fnFolder, sigType{1}, basefn_scope])



tl=200e-12;
fmax=56.4; %%%% The fact 26 was found by doing a fit to the SUT and adjusting for the right frequency axis.


[xsIniRaw,ysIniRaw]=loadCellData([generalFolder_scope, fnFolder, sigType_scope{1}, basefn_scope]);
[xTLSRaw,yTLSRaw]=loadCellData([generalFolder_scope, fnFolder, sigType_scope{6}, basefn_scope]);


% load([datDir fnSpec]);xsIni=x1-1010e-9;ysIni=y1;dts=mean(diff(xsIni));
% load([datDir fnSUToptical]);xSUTIni=x1;ySUTIni=y1;dtSUT=mean(diff(xSUTIni));
% load([datDir fnSUT]);xSUTIni=x2-1010e-9;ySUTIni=y2;dtSUT=mean(diff(xSUTIni));
% load([datDir fnTLS]); xTLS=x3-1010e-9; yTLS=y3;

noiseRegion=1:2e3;
% meanNoiseBackground=mean(ysIniRaw((noiseRegion))); stdNoiseBackground=std(ysIniRaw((noiseRegion)));

xsIni=xsIniRaw(intOfInterest); 

meanNoiseBackground=min(ysIniRaw(intOfInterest));%0.011;mean(ysIniRaw((noiseRegion)));
stdNoiseBackground=std(ysIniRaw((noiseRegion)));
ysIni=(ysIniRaw(intOfInterest)-meanNoiseBackground);
normFac=mean(findpeaks(ysIni,'MinPeakDistance',250));
ysIni=ysIni/normFac;%max(ysIni);%/stdNoiseBackground;%-min(ysIni(intOfInterest));

ysIni1=ysIni;


figure;plot((xsIniRaw(intOfInterest)-mean(xsIniRaw(intOfInterest)))*1e9,ysIniRaw(intOfInterest)*1e3)
xlabel('Time (ns)'); ylabel('Photovoltage (mV)')
% ysIni=real(filtSG(ysIni1',2e3,1,1));

%% Process Spectrogram Data
% Interpolate data and adjust length to
nInterp=1000;%1e1;
lent=numel(xsIni)*nInterp;
xInterp=linspace(xsIni(1), xsIni(end),lent);
yInterp=interp1(xsIni,ysIni,xInterp,'spline');
dtInterp=xInterp(2)-xInterp(1);
ntl=round(tl/dtInterp)
nLenses=floor(lent/ntl)
ys=yInterp(1:nLenses*ntl); xs=xInterp(1:nLenses*ntl)-mean(xInterp);
ys=ys-min(ys);
spgm=reshape(ys,ntl,nLenses);

% Center overal spectrogram
halfShift=1;thresh=10;
[spgm,ys]=centerSpectrogramF(spgm,ntl,nLenses,ys,halfShift,thresh);
spgm=flipud(spgm) ;
% ys=circshift(ys,-244);
% spgm=reshape(ys,ntl,nLenses);
% 
mf=sum(spgm,2);




mfNorm=mf-min(mf); 
mfNorm=mfNorm/max(mfNorm);
fSpecGHz=linspace(-fmax/2,fmax/2,numel(mf));
tSpecns=((1:nLenses)-nLenses/2)*tl*1e9;
figure;
imagesc(spgm)
% flipSpgm=1;
% if flipSpgm
%     spgm=flipud(spgm);
% end
% % 





% OSA 

[lambda,power]=loadCellData([generalFolder_scope, apexFolder fnFolder,basefn_OSA, sigType_OSA{1}, basefn_OSAtime]);
c=299792458;
freqOSA_G=c./lambda;
power_f=fliplr(power);
[~,maxFreqInd]=max(power_f);
powerf_lin=10.^(power_f/10)/max(10.^(power_f/10));
filtN=5e2
powerf_linNorm=real(filtSG(powerf_lin,filtN,1,1));powerf_linNorm=powerf_linNorm/max(powerf_linNorm)
centerFreq=freqOSA_G*powerf_linNorm'/(sum(powerf_linNorm));
freqCent=freqOSA_G-centerFreq;
percent=0.1;
rise=find(powerf_linNorm>percent,1); fall=find(powerf_linNorm>percent,1,'last');

risemfNorm=find(mfNorm>percent,1); fallmfNorm=find(mfNorm>percent,1,'last');

% reductionAmount=3;
% lenOSAfin=floor(numel(powerf_linNorm)/reductionAmount);
% powerf_linNorm2=zeros(lenOSAfin,1);
% for i=1:lenOSAfin
%     powerf_linNorm2(i)=sum(powerf_linNorm((1:reductionAmount)+(i-1)*reductionAmount));
% end


figure;
plot(freqCent,powerf_lin)
hold on
plot(freqCent,powerf_linNorm)
plot(freqCent(rise:fall),  powerf_linNorm(rise:fall),'.')
yyaxis right
plot(fSpecGHz,mfNorm)
hold on
plot(fSpecGHz(risemfNorm:fallmfNorm),mfNorm(risemfNorm:fallmfNorm),'.')




figure;
% subplot(2,1,1)
% imagesc(tSpecns,fSpecGHz,spgm/max(max(spgm)))
imagesc(tSpecns,fSpecGHz,10*log(spgm/max(max(spgm))))

map=jet(500);
colormap(map)
caxis([-30 0])

colorbar()

figure
% subplot(2,1,2)
plot(freqCent,powerf_linNorm)
hold on
plot(fSpecGHz,mfNorm)
xlimFrac=1.5;
xlim([xlimFrac*fSpecGHz(1) xlimFrac*fSpecGHz(end)])
legend(['OSA trace, 90%FW: ' num2str(abs(freqCent(rise))+abs(freqCent(fall))) ],...
['spgm Marginal trace, 90%FW: ' num2str(abs(fSpecGHz(risemfNorm))+abs(fSpecGHz(fallmfNorm)))])
