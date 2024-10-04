function devbar = dev_grid(Vmat,t,trelax,K,M)

% Normalize V
for i = 1:K
    for j = 1:M
        Vmin = min(Vmat(:,i,j));
        Vmax = max(Vmat(:,i,j));
        U(:,i,j) = (Vmat(:,i,j)-Vmin)/(Vmax-Vmin);
    end
end

% Get population average and find deviation post-relaxation
Upop = mean(mean(U,3),2);
for i = 1:K
    for j = 1:M
        dU(:,i,j) = abs(U(:,i,j)-Upop);
    end
end
tr = t(t > trelax);
Tr = tr(end)-tr(1);
dUr = dU(t > trelax,:,:);

% Deviation score: average deviation over time
dt = diff(tr);
for i = 1:K
    for j = 1:M
        dUbar = (dUr(2:end,i,j)+dUr(1:end-1,i,j))/2;
        dev(i,j) = dUbar'*dt/Tr;
    end
end
devbar = mean(dev(:));
