function  [x_new, y_new] = timestepping_FEuler(x_old,  y_old, dt, veloc_fldStruct, kappa_fldStruct)
%% Forward Euler Timestepping

% Specifications:

% Input:
% x_old, y_old: initial positions of all particles (Row vector)
% dt: increment in time
% nts: number of timesteps

% Output:
% [x_new, y_new]: Positions of all particles

%% computation
nparticles = length(x_old);
% assert(all(size(x_old) == [1, nparticles]))

u = veloc_fldStruct.u;
kappa = kappa_fldStruct.kappa;
Divkappa = kappa_fldStruct.Divkappa;

% Part 1: Deterministic drift
% Advective drift size([2, part]); 2 stands for x and y direction
advective_drift = u(x_old, y_old);
% assert(all(size(advective_drift) == [2, nparticles]))

% Diffusive drift due to Anti-Ito
Divkappa_drift = Divkappa(x_old, y_old);
% assert(all(size(Divkappa_drift) == [2, nparticles]))

deterministic_drift = (advective_drift + Divkappa_drift)*dt;

% Part 2: Stochastic drift
% stochastic_drift = sqrt(dt)*sqrt(2*kappa(x_old, y_old)).*randn([2, nparticles]);
% Reminder/Be careful:
% Dimension of sqrt2kappa(x_old, y_old) = row vector of length nparticles
% and randn([2, nparticles])=matrix of size(2, nparticles) are NOT the same.

stochastic_drift = zeros(2, nparticles);
for part = 1:nparticles
    kappa_tensor = kappa(x_old(part), y_old(part));
    K_sigma = 2*dt*kappa_tensor;
    
    stochastic_drift(:,part) = mvnrnd([0,0],K_sigma,1)';
end

% Total_drift
drift = deterministic_drift + stochastic_drift;

% Time stepping: advect and diffuse the particle
x_new = x_old + drift(1,:);
y_new = y_old + drift(2,:);


end
