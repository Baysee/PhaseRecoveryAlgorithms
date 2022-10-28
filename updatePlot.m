function updatePlot(h1,xlimsZoom,t,tspgm,fspgm,xt0,SUT,xt,spgm,SiIntensity,diC,diCLabel,diR,diRLabel)
figure(h1);
FS=12;
subplot(4,4,[1,2,5,6])
imagesc(tspgm,fspgm,spgm)
title('SUT spgm')
subplot(4,4,[3,4,7,8])
imagesc(tspgm,fspgm,SiIntensity)
title('current Iteration spgm')
subplot(4,4,[9,10])
plotIniFin(t,xt0,xt,SUT,FS)
subplot(4,4,[11,12])
plotIniFin(t,xt0,xt,SUT,FS)
ylims=ylim(); xlim(xlimsZoom); ylim(ylims)

subplot(4,4,[13:16])
yyaxis left
plot(diC)
ylabel(diCLabel)
yyaxis right
% plot(20*log10(diC)); hold on
plot(diR)
ylabel(diRLabel)
% legend('Convergence Curve (GL)','spectral convergence');%,'reconstruction error')
xlabel('Iteration'); 
set(gca,'FontSize',FS)
drawnow();
%
end












