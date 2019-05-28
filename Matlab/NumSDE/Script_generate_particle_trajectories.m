% Generate Simulated Trajectories
addpath(genpath('../../../Bayesian_SDE_Inference/Two_Dimensional/Scripts'))
addpath(genpath('../../../Bayesian_SDE_Inference/Two_Dimensional/StructConversion'))
%% Estimation of a variable diffusivity from a Brownian trajectory

% RULE OF THUMB:
%  (t, dim, *sample*)
% 1) Increase in time = Increase in row
% 2) If stationary, number of row <=> along x,y,z

% 3) Gradient of one scalar is ALWAYS row vectors at first. If it involves mU1tiple sampes then treated in case by case basis.

%clear all
close all
%addpath('../StructConversion/')

rng(0)

%% Ideal testcase
T = 12*365.25*24*3600;
%T = 120*24*3600;
Nts = T/(24*3600);  % 12 years, daily record
Nts_sub = (24*3600)/(2*3600);  % sub timestep = 2 hours
% disp('One year data!')
% T = 365*24*3600;
% Nts = 12;  % 10 years, 6 record in each of the 12 months 
L = 3840 * 1000;  % in metre

PART_ranges_min = [0, 0];
PART_ranges_max = [L, L];
kappa_max = 10000;  % m^2/s
U_max = 0.5; % m/s

% Rescale to [0, 1] domain
sc = 1/L;
t_0 = 0; t_1 = T*sc;
kappa_scale = kappa_max*sc;
PART_ranges_min = PART_ranges_min*sc;
PART_ranges_max = PART_ranges_max*sc;

% Initial allocation of particles
npart_pc_DIR = [1,1];
ncell_DIR = [25, 25];
%ncell_DIR = ncell_DIR*4;

t_range = [t_0, t_1];

%% Choose configuration
% kappa_Profile = 'TG_iso';
kappa_Profile = 'sinusoidal';
%kappa_Profile = 'const';
veloc_Profile = 'TTG';
RunProfile = ['K_', kappa_Profile];
RunProfile = [veloc_Profile, '_', kappa_Profile];

% Reminder:
% kappa: diffusivity for FK
% sigma: diffusion coeff in SDE

param.Pe = 20;
if (contains(veloc_Profile,'10') == 1)
  param.Pe = 20/10;  
end
if (contains(veloc_Profile,'100') == 1)
  param.Pe = 20/100;  
end
param.kappa_scale = kappa_scale;
param.U_scale = param.Pe*kappa_scale/(L*sc);
param.GradU_scale = param.Pe*kappa_scale/(L*sc);

kappa_fldStruct = Prescribed_DiffusivityTensorField(kappa_Profile, param);
veloc_fldStruct = Prescribed_VelocityField(veloc_Profile, param);

%% Visualisatoin of velocity field and diffusivty field
fld_vis = 1;
if (fld_vis == 1)
    % Visualise velocity field
    x_c = linspace(PART_ranges_min(1), PART_ranges_max(1), 128);
    y_c = linspace(PART_ranges_min(2), PART_ranges_max(2), 128);
    [X, Y] = meshgrid(x_c, y_c);

    psi = zeros(size(X));
    U = zeros(size(X));
    V = zeros(size(X));
    K11 = zeros(size(X));
    K22 = zeros(size(X));
    K12 = zeros(size(X));
    for pt = 1:length(X(:))
        u_c = veloc_fldStruct.u(X(pt), Y(pt));
        U(pt) = u_c(1);
        V(pt) = u_c(2);
        psi(pt) = veloc_fldStruct.psi(X(pt), Y(pt));
    
        K_c = kappa_fldStruct.kappa(X(pt), Y(pt));
        K11(pt) = K_c(1,1);
        K22(pt) = K_c(2,2);
        K12(pt) = K_c(1,2);
    end
      
    figure(1)
    subplot(2,3,1)
    contourf(X, Y, U)
    caxis([min(U(:)), max(U(:))]);
    colorbar;
    subplot(2,3,2)
    contourf(X, Y, V)
    caxis([min(V(:)), max(V(:))]);
    colorbar;
    subplot(2,3,3)
    contourf(X, Y, psi)
    caxis([min(psi(:)), max(psi(:))]);
    colorbar;
    
    subplot(2,3,4)
    contourf(X, Y, K11)
    subplot(2,3,5)
    contourf(X, Y, K22)
    subplot(2,3,6)
    contourf(X, Y, K12)
end

stop

%% Simulations
disp('Numerical Scheme: Forward Euler');
[x, y, ts_list] = generate_particle_trajectories(t_range, Nts, Nts_sub, ...
    npart_pc_DIR, ncell_DIR, PART_ranges_min, PART_ranges_max, ...
    veloc_fldStruct, kappa_fldStruct);

%% Preview trajectories
fld_vis = 0;
if (fld_vis == 1)
    figure(2)
    subplot(2,2,1)
    scatter(x(1, :),y(1, :))
    
    subplot(2,2,2)
     plot(x(:,1:160:end),y(:,1:160:end))
%    plot(x(:,1:4:end),y(:,1:4:end))
    
    subplot(2,2,3)
    diffx = x(end,:)-x(1,:);
    diffy = y(end,:)-y(1,:);
    histogram(diffx)
    hold on
    histogram(diffy)
    
    subplot(2,2,4)
    h2 = histogram2(x(end, :), y(end, :), 10, 'DisplayStyle','tile','ShowEmptyBins','on');
    h2.DisplayStyle = 'tile';
    view(2)
    
    K = cov(diffx, diffy)/(2*t_1);
    disp(K)
    disp(kappa_fldStruct.kappa(0, 0))
end


%% Output
disp('Saving the trajectories!')
filepath = ['/data/tying/NumSDE/', RunProfile, '/PART_TRAJ/'];
[status, msg, msgID] = mkdir(filepath);

Npart_file = prod(npart_pc_DIR .* ncell_DIR);
TempIntv_days_file = T/(24*3600);  % Hardcoded
SamplingInterval_days_file = 1;  % Hardcoded

traj_filename = ['traj_Npart', num2str(Npart_file), ...
    '_TempIntv', num2str(TempIntv_days_file), ...
    '_h',num2str(SamplingInterval_days_file), ...
    '_nReal1', '.mat'];

save([filepath, traj_filename], 'x', 'y', 'ts_list', 'kappa_fldStruct', ...
                                'veloc_fldStruct', 'param', 't_1','-v7.3')


%% Downsample output
run('Script_LoadTrajData');
disp('Finished Loading Raw Trajectory Data.')

filepath = ['../../trajdat/',RunProfile, '/'];
[status, msg, msgID] = mkdir(filepath);

for h_days = [5, 15, 30, 60, 90]
    h = h_days* (24*3600) * sc;
    
    dt_ind = h/(ts_list(2)- ts_list(1));
    
    xS = x(1:round(dt_ind):end, :);
    yS = y(1:round(dt_ind):end, :);
    ts_listS = ts_list(1:round(dt_ind):end, :);
    
    filenameS = [filepath, 'h', num2str(h_days),'d'];
    write_traj(xS, yS, ts_listS, t_1, filenameS);
end


