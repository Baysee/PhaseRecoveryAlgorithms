

function plotIni(t,xt0,SUT,FS)
plot(t,real(xt0)); hold on; 
plot(t,real(SUT),':')
legend('initial guess','SUT')
ylabel('real part magnitude (a.u.)');
xlabel('Time (a.u.)')
set(gca,'FontSize',FS)
end
