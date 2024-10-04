clear all

dd = '../../dat/';
load([dd 'ml_fig6.mat'])
fs = 20;

figure(1); clf; hold on
imagesc(ps,T0s,devbar)
xlim([min(ps) max(ps)])
ylim([min(T0s) max(T0s)])
xlabel('Edge probability, p','fontsize',fs)
ylabel('Driving period, T_0 (s)','fontsize',fs)
title('Deviation score','fontsize',fs)
colorbar
set(gca,'ydir','normal','layer','top','fontsize',fs)
box on

dnorm = (devbar-min(devbar(:)))/(max(devbar(:))-min(devbar(:)));

figure(2); clf
bar3(T0s,dnorm)
xlabel('connectivity, $p$','interpreter','latex','fontsize',fs)
ylabel('period, $T$ (s)','interpreter','latex','fontsize',fs)
zlabel('norm.\ deviation score','interpreter','latex','fontsize',fs)
set(gca,'fontsize',fs,'xticklabels',ps)
box on

fd = '../../fig/';
print(gcf,'-depsc',[fd 'ml_fig6.eps'])
