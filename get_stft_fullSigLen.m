function stft=get_stft_fullSigLen(nIncs,windowCenters,lent,win,dt,SUT)


stft=zeros(lent,nIncs);

for i=1:nIncs
%     sutInds=mod(winInds+windowCenters(i),lent)+1;
    stft(:,i)=nfft(SUT.*circshift(win,windowCenters(i)),dt);
end

end