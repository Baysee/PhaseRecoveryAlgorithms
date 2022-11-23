% manualInconsistencyCalc

figDir='/Users/ben/Desktop/PhD/manuscripts/Journals/TLS/phaseRecovery_TLS/figs/phaseRecovFigs/'

fnIn='oso_spgmIn.fig';
fnOut='oso_spgmOut.fig';


sIn=open([figDir fnIn]);
sIn=getimage(sIn);

sOut=open([figDir fnOut]);
sOut=getimage(sOut);

sIn=sIn/sum(sum(sIn));
sOut=sOut/sum(sum(sOut));


sqrt(sum(sum(abs(sqrt(sIn)-sqrt(sOut)).^2))/sum(sum(abs((sOut)))))