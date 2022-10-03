function ispgm=get_istft_winInds(lent,winInds,windowCenters,analysisWin,Fs,nIncs,S)
% Reconstruct temporal waveform from spgm by the overlap and add technique

ispgm=zeros(1,lent); % Initialize variable

for i=1:nIncs
    sutInds=mod(winInds+windowCenters(i),lent)+1;
    ispgm(sutInds)=ispgm(sutInds)+analysisWin.*nifft(S(:,i),Fs).';
end
% ispgm=ispgm;
end

