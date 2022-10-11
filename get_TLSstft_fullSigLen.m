function stft=get_TLSstft_fullSigLen(f,lent,winL,TL,fmax,dt,SUT)


CL=2*pi*fmax/TL;
nLenses=lent/winL;
phaseFun=repmat(exp(1j*((1:winL)/winL*TL-TL/2).^2/2*CL),[1,nLenses]);
dispFun=exp(1j*(2*pi*f).^2/(2*CL));

stft1D=nifft(dispFun.*nfft(phaseFun.*SUT,dt),1/dt);
stft=reshape(stft1D,[winL,nLenses]);
% stft=zeros(lent,nIncs);
% 
% for i=1:nIncs
% %     sutInds=mod(winInds+windowCenters(i),lent)+1;
%     stft(:,i)=nfft(SUT.*circshift(win,windowCenters(i)),dt);
% end

end