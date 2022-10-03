function updatePlot(h1,xlimsZoom,t,tspgm,fspgm,xt0,SUT,xt,spgm,SiIntensity,diC,diCLabel,diR,diRLabel)
figure(h1);
FS=16;
subplot(3,2,1)
imagesc(tspgm,fspgm,spgm)
title('SUT spgm')
subplot(3,2,2)
imagesc(tspgm,fspgm,SiIntensity)
title('current Iteration spgm')
subplot(3,2,3)
plotIniFin(t,xt0,xt,SUT,FS)
subplot(3,2,4)
plotIniFin(t,xt0,xt,SUT,FS)
ylims=ylim(); xlim(xlimsZoom); ylim(ylims)

subplot(3,2,5:6)
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












