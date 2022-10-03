function stft=get_stft_winInds(nIncs,winInds,windowCenters,lent,win,dt,winLen,SUT)
% By taking only the indeces of the widow width, the spectrum appears very
% rough (low frequency resolution). Althoguh this seems to be the most
% common approach online, it doesn't seem to work very well for me here. 
stft=zeros(winLen,nIncs);

for i=1:nIncs

    sutInds=mod(winInds+windowCenters(i),lent)+1;
    stft(:,i)=nfft(SUT(sutInds).*win,dt);
end

end



