datDir='C:\Users\Lord Photon\Documents\MATLAB\time-frequency analysis\outputFigs\'
datMatName='dataMatrix21h58m.mat';

load([datDir datMatName])
% % define loop parameters
clips_l=[20,40,60];
nptPerWin_l=[16,32,64];
bwIRF_l=[30:5:70]*1e9;
interpAmounts=[1,2,4,8,16];
% (icl, inptWin, ibwIRF, iInter);
b21=-12758; b22=5082;

dim1=1;dim2=2;
datToPlot1=squeeze(alldiRs(dim1,dim2,:,:));plotTitle1='Frog error';
datToPlot2=squeeze(100*abs(allMeanb2s1(dim1,dim2,:,:)-b21)/abs(b21));plotTitle2='b2 error';
datToPlot3=squeeze((allStdb2s1(dim1,dim2,:,:)));plotTitle3='std(b2)';
datToPlot4=squeeze(100*abs(allMeanb2s2(dim1,dim2,:,:)-b22)/abs(b22));plotTitle2='b2 error';
datToPlot5=squeeze((allStdb2s2(dim1,dim2,:,:)));plotTitle3='std(b2)';

figure;
subplot(3,2,1:2);
imagesc(interpAmounts,bwIRF_l,datToPlot1);colorbar(); title(plotTitle1)
subplot(3,2,3);
mesh(interpAmounts,bwIRF_l,datToPlot2);colorbar(); title(plotTitle2)
subplot(3,2,5);
mesh(interpAmounts,bwIRF_l,datToPlot3);colorbar(); title(plotTitle3)

subplot(3,2,4);
mesh(interpAmounts,bwIRF_l,datToPlot4);colorbar(); title(plotTitle2)
subplot(3,2,6);
mesh(interpAmounts,bwIRF_l,datToPlot5);colorbar(); title(plotTitle3)