% Triangular lattice (6 neighbors)
clear all
fs = 16;
lw = 3;

% Number of cells per dimension
K = 18; % rows
M = 18; % columns

Y = 4; % number of looped values
R = 10; % number of random samplings (disrupted nodes, initial conditions)

% Base parameters
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

% Choose I barely in excitable regime
I = 39.5; % uA/cm^2

% Parameter vector and tensor of cell matrix's parameters
theta0 = [CM gCa gK gL VCa VK VL V1 V2 V3 V4 phi I];
for i = 1:length(theta0)
    theta(:,:,i) = theta0(i)*ones(K,M);
end

% Choose gamma in the weak coupling regime (<~ 1/T0 ~ .02)
gamma_ = .002;

% Base pulse strength (smallest that produces excitation, ~5 mV)
Vp0 = 5; % mV
H = 1; % max percent change for heterogeneity

% Driving parameters
Z = 20; % total number of periods
Zr = 10; % number periods for relaxation time
T0s = linspace(40,200,Y); % s

% Fraction of disrupted nodes
qs = linspace(0,.75,Y);

% Loop over T0
for i = 1:length(T0s)
    T0 = T0s(i)
    T = Z*T0;
    trelax = Zr*T0;

    % Loop over q
    for j = 1:length(qs)
        q = qs(j)

        % Loop over random samplings (nodes, ICs, Vp values)
        for k = 1:R
            
            % Heterogeneous pulse strengths
            Vp = Vp0*(1+H*(2*rand(K,M)-1));

            % Random node disruption (with probability q)
            Dr = ones(K,M-1); % rightward edges (bidirectional)
            Dd = ones(K-1,M); % downward edges (bidirectional)
            Ddr = ones(K-1,M-1); % down-and-rightward edges (bidirec.)
            Q = rand(K,M) < q; % disrupted nodes
            Dr(Q(:,1:end-1)) = 0;
            Dr(Q(:,2:end)) = 0;
            Dd(Q(1:end-1,:)) = 0;
            Dd(Q(2:end,:)) = 0;
            Ddr(Q(1:end-1,1:end-1)) = 0;
            Ddr(Q(2:end,2:end)) = 0;
            
            % Random initial conditions V ~ -20 mV (N = 0)
            V0mat = -20*ones(K,M) + 5*rand(K,M);
            N0mat = zeros(K,M);
            x0 = [reshape(V0mat,1,K*M) reshape(N0mat,1,K*M)];

            % Dynamics
            t = []; x = []; x0z = x0;
            for z = 1:Z
                [tz,xz] = ode45(@(t,x) ...
                    ml_triang_rhs(t,x,theta,gamma_,Dr,Dd,Ddr,K,M),...
                    [0 T0],x0z);
                t = [t; tz + (z-1)*T0];
                x = [x; xz];
                x0z = xz(end,:) + [Vp(:)' zeros(1,K*M)];
            end
            NT = length(t);
            Vvec = x(:,1:K*M);
            Vmat = reshape(Vvec,NT,K,M);
            Nvec = x(:,K*M+1:end);
            Nmat = reshape(Nvec,NT,K,M);

            dev(i,j,k) = dev_grid(Vmat,t,trelax,K,M);
        end
    end
end
devbar = mean(dev,3);

figure(1); clf; hold on
imagesc(qs,T0s,devbar)
xlim([min(qs) max(qs)])
ylim([min(T0s) max(T0s)])
xlabel('Disruption fraction, q','fontsize',fs)
ylabel('Driving period, T_0 (s)','fontsize',fs)
title('Deviation score','fontsize',fs)
colorbar
set(gca,'ydir','normal','layer','top','fontsize',fs)
box on

dd = '../../dat/';
save([dd 'ml_fig4c.mat'])