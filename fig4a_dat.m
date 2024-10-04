clear all
fs = 20;
lw = 3;

CM = 20; % uF/cm^2
gCa = 4; % mS/cm^2
gK = 8; % mS/cm^2
gL = 2; % mS/cm^2
VCa = 120; % mV
VK = -80; % mV
VL = -60; % mV
V1 = -1.2; % mV
V2 = 18; % mV
V3 = 12; % mV
V4 = 17.4; % mV

% Choose phi near SNHO point
phi = 0.184; % 1/s

% Choose I in excitable regime
I = 39.5; % uA/cm^2

% Driving parameters
Vp = 5; % mV
T0 = 80; % s
Z = 20;
T = Z*T0;

% Initial conditions
V0A = -20; % mV
N0A = 0;
V0B = -10; % mV
N0B = 0;
x0 = [V0A N0A V0B N0B];

% Case a: weak coupling
gammaa = .003;
ta = []; xa = []; x0z = x0;
for z = 1:Z
    [tz,xz] = ode45(@(t,x) ...
        ml_pair_rhs(t,x,CM,gCa,gK,gL,VCa,VK,VL,V1,V2,V3,V4,...
        phi,I,gammaa),[0 T0],x0z);
    ta = [ta; tz + (z-1)*T0];
    xa = [xa; xz];
    x0z = xz(end,:) + [Vp 0 Vp 0];
end
VAa = xa(:,1);
NAa = xa(:,2);
VBa = xa(:,3);
NBa = xa(:,4);

figure(1); clf
plot(ta,VAa,'b-',ta,VBa,'r--',(1:Z)*T0,max([VAa;VBa]),'kv',...
    'linewidth',lw)
xlim([.7*T T])
%ylim([-40 30])
xlabel('Time, t (s)')
ylabel('Voltage, V (mV)')
title(['Weak coupling: \gammaT_0 = ' num2str(gammaa*T0)])
set(gca,'fontsize',fs)

% Case b: strong coupling
gammab = .05;
tb = []; xb = []; x0z = x0;
for z = 1:Z
    [tz,xz] = ode45(@(t,x) ...
        ml_pair_rhs(t,x,CM,gCa,gK,gL,VCa,VK,VL,V1,V2,V3,V4,...
        phi,I,gammab),[0 T0],x0z);
    tb = [tb; tz + (z-1)*T0];
    xb = [xb; xz];
    x0z = xz(end,:) + [Vp 0 Vp 0];
end
VAb = xb(:,1);
NAb = xb(:,2);
VBb = xb(:,3);
NBb = xb(:,4);

figure(2); clf
plot(tb,VAb,'b-',tb,VBb,'r--',(1:Z)*T0,max([VAb;VBb]),'kv',...
    'linewidth',lw)
xlim([.7*T T])
%ylim([-40 30])
xlabel('Time, t (s)')
ylabel('Voltage, V (mV)')
title(['Strong coupling: \gammaT_0 = ' num2str(gammab*T0)])
set(gca,'fontsize',fs)

save('../../dat/ml_fig4a.mat')
