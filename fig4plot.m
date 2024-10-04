clear all
dd = '../../dat/';
fd = '../../fig/';
blueT = [37 90 146]/256;
bluep = [128 128 247]/256;
redp = [239 135 132]/256;
brownq = [223 199 159]/256;
blueq = [192 228 245]/256;

% A: coupled pair
load([dd 'ml_fig4a.mat'])
figure(1); clf
lw = 2; lw2 = .5;
fs = 15; fs2 = 18;

Toff = .8*T % offset time (post-equilibration)

[ig,j] = min(abs(ta-Toff));
ta_ = ta(j:end)-ta(j);
VAa_ = VAa(j:end);
VBa_ = VBa(j:end);

subplot(3,2,2)
plot(ta_/T0,VAa_,'b-',ta_/T0,VBa_,'r--','linewidth',lw)
xlim([0 max(ta_)/T0])
ylim([min(VAa_)-5 max(VAa_)+5])
set(gca,'fontsize',fs2,'xtick',0:4,'xticklabel',[])
ylabel('$V_1$, $V_2$','interpreter','latex','fontsize',fs2)
title('Weak coupling, $\gamma < 1/T$','interpreter','latex','fontsize',fs2)

[ig,j] = min(abs(tb-Toff));
tb_ = tb(j:end)-tb(j);
VAb_ = VAb(j:end);
VBb_ = VBb(j:end);

subplot(3,2,4)
plot(tb_/T0,VAb_,'b-',tb_/T0,VBb_,'r--','linewidth',lw)
xlim([0 max(tb_)/T0])
ylim([min(VAb_)-5 max(VAb_)+5])
set(gca,'fontsize',fs2,'xtick',0:4)
xlabel('Time, $t/T$','interpreter','latex','fontsize',fs2)
ylabel('$V_1$, $V_2$','interpreter','latex','fontsize',fs2)
title('Strong coupling, $\gamma > 1/T$','interpreter','latex','fontsize',fs2)

print(gcf,'-depsc',[fd 'ml_fig4a.eps'])

% B: 18x18 square grid
load([dd 'ml_fig4b.mat'])
lw = 1.5; lw2 = .5;
fs = 15; fs2 = 18;

dT = mean(squeeze(dev(:,1,:)),2);
dTerr = std(squeeze(dev(:,1,:)),[],2);

figure(2); clf
subplot(2,2,1); hold on
bar(T0s,dT,'linewidth',lw2,'facecolor',blueT)
er = errorbar(T0s,dT,dTerr,'k','linewidth',lw2);
er.LineStyle = 'none';
ylim([0 .25])
xlim([10 230])
xlabel('period, $T$ (s)','interpreter','latex','fontsize',fs)
ylabel('deviation score','interpreter','latex','fontsize',fs)
set(gca,'fontsize',fs,'xtick',T0s)
box on

print(gcf,'-depsc',[fd 'ml_fig4b1.eps'])

ps = [0 .75];
d40 = mean(squeeze(dev(1,[1 4],:)),2);
d40err = std(squeeze(dev(1,[1 4],:)),[],2);

figure(3); clf
subplot(2,3,2); hold on
bar(ps(1),d40(1),.6,'linewidth',lw2,'facecolor',bluep)
bar(ps(2),d40(2),.6,'linewidth',lw2,'facecolor',redp)
er = errorbar(ps,d40,d40err,'k','linewidth',lw2);
er.LineStyle = 'none';
xlim([-.5 1.25])
ylim([0 .25])
xlabel('$\qquad\qquad\quad$ connectivity, $p$','interpreter','latex','fontsize',fs)
ylabel('deviation score','interpreter','latex','fontsize',fs)
title('$T = 40$ s','interpreter','latex','fontsize',fs)
set(gca,'xtick',ps,'fontsize',fs)
box on

d200 = mean(squeeze(dev(5,[1 4],:)),2);
d200err = std(squeeze(dev(5,[1 4],:)),[],2);

subplot(2,3,3); hold on
bar(ps(1),d200(1),.6,'linewidth',lw2,'facecolor',bluep)
bar(ps(2),d200(2),.6,'linewidth',lw2,'facecolor',redp)
er = errorbar(ps,d200,d200err,'k','linewidth',lw2);
er.LineStyle = 'none';
xlim([-.5 1.25])
ylim([0 .25])
%xlabel('connectivity, $p$','interpreter','latex','fontsize',fs)
%ylabel('deviation score','interpreter','latex','fontsize',fs)
title('$T = 200$ s','interpreter','latex','fontsize',fs)
set(gca,'xtick',ps,'fontsize',fs,'yticklabels',[])
box on

print(gcf,'-depsc',[fd 'ml_fig4b2.eps'])

% C: 18x18 triangular lattice (parallelogram)
load([dd 'ml_fig4c.mat'])
lw = 1.5; lw2 = .5;
fs = 15; fs2 = 18;

qs = [0 .5];
d = mean(squeeze(dev(4,[1 3],:)),2);
derr = std(squeeze(dev(4,[1 3],:)),[],2);

figure(4); clf
subplot(2,3,2); hold on
bar(qs(1),d(1),.4,'linewidth',lw2,'facecolor',brownq)
bar(qs(2),d(2),.4,'linewidth',lw2,'facecolor',blueq)
er = errorbar(qs,d,derr,'k','linewidth',lw2);
er.LineStyle = 'none';
xlim([-.3 .8])
ylim([0 .25])
xlabel('disruption, $q$','interpreter','latex','fontsize',fs)
ylabel('deviation score','interpreter','latex','fontsize',fs)
title('$T = 200$ s','interpreter','latex','fontsize',fs)
set(gca,'xtick',qs,'fontsize',fs)
box on

print(gcf,'-depsc',[fd 'ml_fig4c.eps'])
