function [xt]=istft_TLS(stft,ntl,phit,phiw)

xt_stft=reshape(stft,[1 numel(stft)]);
xw_stft=nfft(xt_stft,1); 
xw_phit=xw_stft.*exp(-1j*phiw);
xt_phi1=nifft(xw_phit,1);
xt=xt_phi1.*exp(-1j*phit);
% 
% xt_phit=xt.*exp(1j*phit);
% xw_phit=nfft(xt_phit,1);
% xw_stft=xw_phit.*exp(1j*phiw);
% xt_stft=nifft(xw_stft,1);
% 
% stft=reshape(xt_stft,[ntl numel(xt)/ntl]);
