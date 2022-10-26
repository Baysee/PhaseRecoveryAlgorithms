function [stft]=stft_TLS(xt,ntl,phit,phiw)

xt_phit=xt.*exp(1j*phit);
xw_phit=nfft(xt_phit,1);
xw_stft=xw_phit.*exp(1j*phiw);
xt_stft=nifft(xw_stft,1);

stft=reshape(xt_stft,[ntl numel(xt)/ntl]);

