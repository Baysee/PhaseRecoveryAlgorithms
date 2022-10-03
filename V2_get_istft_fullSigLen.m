
function ispgm=V2_get_istft_fullSigLen(lent,windowCenters,analysisWin,Fs,nIncs,S)
% Reconstruct temporal waveform from spgm by the overlap and add technique

ispgm=zeros(1,lent); % Initialize variable

% figure;
for i=1:nIncs
    
    ift=nifft(S(:,i),Fs);
%     sutInds=mod(winInds+windowCenters(i),lent)+1;
    ispgm=ispgm+circshift(analysisWin,windowCenters(i)).*ift.';%,;
    
%     
% %     
%     yyaxis left; hold off
% plot(real(ispgm)); hold off
%     yyaxis right; 
%     plot(real(ift)); hold on; 
%     plot(circshift(analysisWin,windowCenters(i))/max(analysisWin)*max(abs(real(ift))))
%     hold off
%     title(num2str(i))
%     drawnow()
end

% ispgm=circshift(ispgm,lent/2);
% ispgm=ispgm;
end
