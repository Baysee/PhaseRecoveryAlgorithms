% analyzeMenlo_simple

% addpath('/Volumes/Extreme SSD/MATLAB_DataAnalysis/analysisDocs')
addpath( '/Users/ben/Documents/MATLAB/library' )

generalFolder='C:\Users\Lord Photon\Documents\MATLAB\time-frequency analysis\Spectrogram_20211108\experimentalData\spectrogram_20211109\';
% generalFolder='/Users/ben/Documents/MATLAB/timeFrequencyAnalysis/experimentalData/spectrogram_20211109/';
apexFolder='ApexControl/20211118/Menlo/62p5ps/';

%% Bandwidth demonstration -- 62.5 ps time lens

fnFolderOSO='20211118_OSO/62p5ps/';
fnFolderOSA=[];
% Relevant dates with menlo data: 14, 16,17, 18, 19

% % % Run 1
% basefn_scope1='408GHz_';
% basefn_scope2='ns_4nsSpan.csv';
% bf1=18; bf2=38; inc=4; 
% nFiles=(bf2-bf1)/inc;
% basefn_OSA='menlo408GHzBW_pm120kmCross_';
% basefn_OSAtime='20h37m.mat';
% alignShift=-78;
% % 
% % % Run 2
basefn_scope1='408GHz_pm120kmCross_';
basefn_scope2='ns_2p2nsSpan.csv';
bf1=16; bf2=28; inc=2; 
nFiles=(bf2-bf1)/inc+1;
basefn_OSA='menlo408GHzBW_pm120kmCross_';
basefn_OSAtime='20h37m.mat';
alignShift=140+15*799;

sigType_scope={'spectrogram', 'pmoff','SUT','SUT_TLS','SUT_trigger','TLS','trigger'};
sigType_OSA={'SUT','pmoff','pmon'};

osabasefn=[basefn_OSA sigType_OSA{1} basefn_OSAtime];

% load([generalFolder_scope, fnFolder, sigType{1}, basefn_scope])



tl=62.5e-12;
fmax=448; %%%% The fact 26 was found by doing a fit to the SUT and adjusting for the right frequency axis.

figure;
% [xsIniRaw,ysIniRaw]=loadCellData([generalFolder_scope, fnFolder, sigType_scope{1}, basefn_scope]);
xsIniRaw=[]; ysIniRaw=[];
for iFile=1: nFiles
[dat]=csvread([generalFolder, fnFolderOSO, basefn_scope1 num2str(bf1+(iFile-1)*inc) basefn_scope2],2,0);
dat(:,1)=dat(:,1);%+2.7*iFile;
plot(dat(:,1),dat(:,2))
hold on
if iFile~=1
[firstNonRepeatInd]=find(dat(:,1)>xsIniRaw(end),1);
dat(1:firstNonRepeatInd-1,:)=[];
end

xsIniRaw=[xsIniRaw;dat(:,1)];ysIniRaw=[ysIniRaw;dat(:,2)];

end
% [xTLSRaw,yTLSRaw]=loadCellData([generalFolder_scope, fnFolder, sigType_scope{6}, basefn_scope]);
% 
xsIniRaw=xsIniRaw*1e-12
% load([datDir fnSpec]);xsIni=x1-1010e-9;ysIni=y1;dts=mean(diff(xsIni));
% load([datDir fnSUToptical]);xSUTIni=x1;ySUTIni=y1;dtSUT=mean(diff(xSUTIni));
% load([datDir fnSUT]);xSUTIni=x2-1010e-9;ySUTIni=y2;dtSUT=mean(diff(xSUTIni));
% load([datDir fnTLS]); xTLS=x3-1010e-9; yTLS=y3;
% intOfInterest=1092:6285;
intOfInterest=1:numel(ysIniRaw);%2e4:4e4;
noiseRegion=1:2e3;
% meanNoiseBackground=mean(ysIniRaw((noiseRegion))); stdNoiseBackground=std(ysIniRaw((noiseRegion)));
meanNoiseBackground=0;mean(ysIniRaw((noiseRegion))); stdNoiseBackground=std(ysIniRaw((noiseRegion)));
xsIni=xsIniRaw(intOfInterest); ysIni=(ysIniRaw(intOfInterest)-meanNoiseBackground)/stdNoiseBackground;%-min(ysIni(intOfInterest));


ysIni1=ysIni;
dtIni=mean(diff(xsIni(1:1000))); FsIni=1/dtIni;
nptFilt=2e3
filtBW=6e3;%2*2e3*FsIni/numel(ysIni1)
ysIni=real(filtSG(ysIni1',filtBW,1,1));




%% Process Spectrogram Data
% Interpolate data and adjust length to
nInterp=42;%1e1;
lent=numel(xsIni)*nInterp;
xInterp=linspace(xsIni(1), xsIni(end),lent);
yInterp=interp1(xsIni,ysIni,xInterp,'spline');
dtInterp=xInterp(2)-xInterp(1);
ntl=round(tl/dtInterp)
nLenses=floor(lent/ntl)
ys=yInterp(1:nLenses*ntl); xs=xInterp(1:nLenses*ntl)-mean(xInterp);
spgm=reshape(ys,ntl,nLenses);


%% Enforce 16 GHz resolution

spgm_df=fmax/ntl*1e9;
npts_dfDesired=2*round(1/tl/spgm_df);
measuredIRF=superGauss(0,6.7e-12/2,1,xs,0);
measuredIRF=circshift(measuredIRF,round(numel(measuredIRF)/2));
[a,b]=deconv(measuredIRF,ys);


% [pks,locs,w,p]=findpeaks(ys,'MinPeakProminence',0.5);
% [pks,locs,w,p]=findpeaks(ys);%,'MinPeakProminence',0.4);
% 
% locsKeep=(locs+[1:npts_dfDesired]')-round(npts_dfDesired/2);
% locsKeep=sort(reshape(locsKeep,[numel(locsKeep),1]));
% removeThesePoints=true(size(ys)); removeThesePoints(locsKeep)=0;
% ysPeaksOnly=ys;
% ysPeaksOnly(removeThesePoints)=0;
% 
% figure;plot(ys); hold on; plot(locs,pks); plot(ysPeaksOnly)
% yyaxis right; plot(locs,p)
% ys=ysPeaksOnly;
% ys=abs(filtSG(ys,filtBW*nInterp/4,1,1));

% Center overal spectrogram
halfShift=0;thresh=10;
[spgm,ys]=centerSpectrogramF(spgm,ntl,nLenses,ys,halfShift,thresh);
% 
ys=circshift(ys,alignShift);
spgm=reshape(ys,ntl,nLenses);


mf=sum(spgm,2);
mt=sum(spgm,1);


mtNorm=mt-min(mt); mtNorm=mtNorm/max(mtNorm);
mfNorm=mf-min(mf); mfNorm=mfNorm/max(mfNorm);
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

%% Load OSA data
filtBW=500e2;
[freqCent,powerf_linNorm,lambda,power]=getOSAdat([generalFolder, apexFolder, fnFolderOSA, osabasefn],filtBW);
% powerf_linNorm2=smooth(powerf_linNorm,10);
% figure;plot(powerf_linNorm); hold on; plot(powerf_linNorm2)
freqCent=freqCent*1e9;
percent=0.1;
[rise,fall]=getRiseFall(powerf_linNorm,percent);

% 
% [lambda,power]=loadCellData([generalFolder_scope, apexFolder fnFolder,basefn_OSA, sigType_OSA{1}, basefn_OSAtime]);
% c=299792458;
% freqOSA_G=c./lambda;
% power_f=fliplr(power);
% [~,maxFreqInd]=max(power_f);
% powerf_lin=10.^(power_f/10)/max(10.^(power_f/10));
% powerf_linNorm=real(filtSG(powerf_lin,50e2,1,1));
% centerFreq=freqOSA_G*powerf_linNorm'/(sum(powerf_linNorm));
% freqCent=freqOSA_G-centerFreq;
% rise=find(powerf_linNorm>0.1,1); fall=find(powerf_linNorm>0.1,1,'last');

risemfNorm=find(mfNorm>percent,1); fallmfNorm=find(mfNorm>percent,1,'last');
risemtNorm=find(mtNorm>percent,1); fallmtNorm=find(mtNorm>percent,1,'last');

figure;
% subplot(2,1,1)
% imagesc(tSpecns,fSpecGHz,spgm/max(max(spgm)))
imagesc(tSpecns,fSpecGHz,10*log(spgm/max(max(spgm))))
% colormap('cmrmap')
caxis([-12 0])
colorbar()
xlabel('Time (ns)'); ylabel('Frequency (GHz)')

map=jet(500);
colormap(map)

figure
% subplot(2,1,2)
plot(freqCent*1e-9,powerf_linNorm)
hold on
plot(fSpecGHz,mfNorm)
xlimFrac=1.5;
xlim([xlimFrac*fSpecGHz(1) xlimFrac*fSpecGHz(end)])

xlabel('Frequency (GHz)')
legend(['OSA trace, 90%FW: ' num2str(1e-9*(abs(freqCent(rise))+abs(freqCent(fall)))) ],...
['spgm Marginal trace, 90%FW: ' num2str(abs(fSpecGHz(risemfNorm))+abs(fSpecGHz(fallmfNorm)))])
