

function plotIniFin(t,xt0,xt,SUT,FS)
yyaxis left;
hold off;
plot(t,real(xt0),'--'); hold on; 
plot(t,real(xt),'b');
% plot(t,real(SUT),':')
hold off
yyaxis right;
plot(t,unwrap(angle(xt)))
legend('initial guess','Last iteration','SUT')
ylabel('real part magnitude (a.u.)');
xlabel('Time (a.u.)')
set(gca,'FontSize',FS)
end
