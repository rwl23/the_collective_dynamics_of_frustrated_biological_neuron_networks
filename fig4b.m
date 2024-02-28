function qif_pair
% v1dot = b + v1^2 + g*(v2-v1) until v1 = 1, then v1 = 0
% v2dot = b + v2^2 + g*(v1-v2) until v2 = 1, then v2 = 0

T = 20;
opt2 = optimoptions('fsolve','display','none');
b = fsolve(@(b)period_fun(b,T),1/T^2,opt2);
g1 = 1/3/T;
g2 = 3/T;
vreset = 0;
vmax = 1;
V0_1 = [0; .02];
V0_2 = [0; .1];

Tt = T*10;
opt = odeset('events',@(t,V)max_fun(t,V,vmax));

% low g
t1 = []; V1 = []; V0i = V0_1; tmax = 0; te = 0;
while tmax < Tt
    [ti,Vi,te,Ve,j] = ...
        ode45(@(t,V)qif_fun(t,V,b,g1),...
        [te(end) Tt],V0i,opt);
    if ~isempty(j)
        t1 = [t1; ti];
        V1 = [V1; Vi];
        V0i = Vi(end,:);
        V0i(j) = vreset;
        tmax = t1(end);
    else
        tmax = Tt;
    end
end

% high g
t2 = []; V2 = []; V0i = V0_2; tmax = 0; te = 0;
while tmax < Tt
    [ti,Vi,te,Ve,j] = ...
        ode45(@(t,V)qif_fun(t,V,b,g2),...
        [te(end) Tt],V0i,opt);
    if ~isempty(j)
        t2 = [t2; ti];
        V2 = [V2; Vi];
        V0i = Vi(end,:);
        V0i(j) = vreset;
        tmax = t2(end);
    else
        tmax = Tt;
    end
end

figure(2); clf
lw = 1; fs = 15;

subplot(2,2,1:2)
h = plot(t1/T,V1,'-','linewidth',lw);
xlim([0 8])
ylim([vreset-.1 vmax+.1])
ylabel('v1, v2')
title(['gT = ' num2str(g1*T)])
set(gca,'fontsize',fs)

subplot(2,2,3:4)
h = plot(t2/T,V2,'-','linewidth',lw);
xlim([0 8])
ylim([vreset-.1 vmax+.1])
xlabel('t/T')
ylabel('v1, v2')
title(['gT = ' num2str(g2*T)])
set(gca,'fontsize',fs)

save('../../dat/fig4b.mat')

function F = period_fun(b,T)
F = T - atan(1./sqrt(b))./sqrt(b);

function dVdt = qif_fun(t,V,b,g)
v1 = V(1);
v2 = V(2);
dv1dt = b + v1.^2 + g*(v2-v1);
dv2dt = b + v2.^2 + g*(v1-v2);
dVdt = [dv1dt; dv2dt];

function [value,isterminal,direction] = max_fun(t,V,vmax)
value = V - vmax;
isterminal = [1; 1];
direction = [1; 1];

