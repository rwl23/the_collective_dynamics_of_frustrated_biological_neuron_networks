clear all

% figure settings
figure(1); clf
lw = 2;
ms = 6;
fs = 15; fs2 = 20;

% A: bifurcation diagram
% vdot = b + v^2 until v = vmax, then v = vreset

% bifurcations
b_sn = [0 0];
vreset_sn = [0 1];
b_snic = [0 0];
vreset_snic = [-1 0];
b_sho = linspace(-1,0,100);
vreset_sho = sqrt(-b_sho);
b_snho = 0;
vreset_snho = 0;

% plot: bifurcation diagram
subplot(2,2,1)
h = plot(b_sn,vreset_sn,'k-',b_snic,vreset_snic,'k-',...
    b_sho,vreset_sho,'k-',b_snho,vreset_snho,'ko',...
    'linewidth',lw);
set(h(4),'markersize',ms,'markerfacecolor','k')
xlim([-1 1])
ylim([-1 1])
set(gca,'fontsize',fs,'ytick',-1:1)
xlabel('$b$','interpreter','latex','fontsize',fs2)
ylabel('$v_{\rm reset}$','interpreter','latex','fontsize',fs2)

% B: coupled pair
load('../../dat/fig4b.mat')
lw = 1.5;
fs = 15; fs2 = 20;

subplot(4,2,2)
h = plot(t1/T,V1,'-','linewidth',lw);
set(h(1),'color','b')
set(h(2),'color','r','linestyle','--')
xlim([0 8])
ylim([vreset-.1 vmax+.1])
set(gca,'fontsize',fs,'xticklabel',[])
ylabel('$v_1$, $v_2$','interpreter','latex','fontsize',fs)
title('Weak coupling, $g < 1/\tau$','interpreter','latex','fontsize',fs)

subplot(4,2,4)
h = plot(t2/T,V2,'-','linewidth',lw);
set(h(1),'color','b')
set(h(2),'color','r','linestyle','--')
xlim([0 8])
ylim([vreset-.1 vmax+.1])
set(gca,'fontsize',fs)
xlabel('Time, $t/\tau$','interpreter','latex','fontsize',fs)
ylabel('$v_1$, $v_2$','interpreter','latex','fontsize',fs)
title('Strong coupling, $g > 1/\tau$','interpreter','latex','fontsize',fs)

fd = '../../fig/';
print(gcf,'-depsc',[fd 'fig4ab.eps'])

figure(2); clf

% C: Uncoupled driven pair
load('../../dat/fig4c.mat')
lw = 2;
fs = 15; fs2 = 20;

T0_ = linspace(0,20,100);
d_ = .25./T0_.^.5;

subplot(2,2,3)
% correct deviation score: should have divided by 2
h = plot(T0s,mean(d,2)/2,'k-',T0_,d_,'b--');
set(h,'linewidth',lw)
ylim([.05 .15])
set(gca,'fontsize',fs,'ytick',.1:.1:.3,'ytick',.05:.05:.15)
xlabel('Driving period, $T$','interpreter','latex','fontsize',fs)
ylabel('Deviation score','interpreter','latex','fontsize',fs)
legend(h(2),{'$\sim1/\sqrt{T}$'},'location','ne',...
    'interpreter','latex','fontsize',fs)

% D: Coupled driven pair
load('../../dat/fig4d.mat')
lw = 2;
fs = 15; fs2 = 20;

subplot(2,2,4)
% correct deviation score: should have divided by 2
plot(gs*T0,mean(d,2)/2,'k-','linewidth',lw)
ylim([.05 .12])
set(gca,'fontsize',fs)
xlabel('Coupling strength, $g$ ($\times T$)',...
    'interpreter','latex','fontsize',fs)
ylabel('Deviation score','interpreter','latex','fontsize',fs)

fd = '../../fig/';
print(gcf,'-depsc',[fd 'fig4cd.eps'])