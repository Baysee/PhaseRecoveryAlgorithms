
function plotFin(t,xt,SUT,FS)

plot(t,real(xt));  hold on; 
plot(t,real(SUT),':')
legend('Converged solution','SUT')
ylabel('real part magnitude (a.u.)');
xlabel('Time (a.u.)')
set(gca,'FontSize',FS)

end
