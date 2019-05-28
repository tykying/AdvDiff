clear all
unittest = 0;

if unittest == 1
    K_fn = @(x) ones(size(x))*0.05;
    dK_fn = @(x) zeros(size(x));
    q0_fn = @(x) exp(-0.5*((x-0.25)/0.05).^2);
    
    n = 16;
    dx = 1/n;
    nstep = 128*n;
    
    q = FVD_solver(n, nstep, K_fn, q0_fn);
    
    plot(q);
    
    q = linspace(-1, 10, 12)+0.5;
    X = 0.32:0.02:0.64';
    assert(all(abs(10*X - intpl_likelihood(X, q)) < 1e-13));
    
    loglik = eval_loglikelihood(X, K_fn, q0_fn);
end

%% Inference
% Truth
Kscale = 0.005;
K0 = 2;
K_fn = @(x) (sin(2*pi*x) + K0)*Kscale;
dK_fn = @(x) 2*pi*cos(2*pi*x)*Kscale;
% K_fn = @(x) exp(-0.5*((x-0.5)/0.2).^2) *Kscale;
% dK_fn = @(x) ((x-0.5)/0.2^2).*exp(-0.5*((x-0.5)/0.2).^2) *Kscale;

q0_fnL = @(x) 1/3.*(x <= 1/3);
q0_fnC = @(x) 1/3.*(x <= 2/3).*(x > 1/3);
q0_fnR = @(x) 1/3.*(x > 2/3);

q0_fnL = @(x) 2.*(x <= 1/2);
q0_fnC = @(x) 2.*(x <= 3/4).*(x > 1/4);  % Not in use
q0_fnR = @(x) 2.*(x > 1/2);
%q0_fn = @(x) (x < 0.55).*(x > 0.45);

% Simulate particles
nstep = 2^12;
npart = 100000;
X0L = init_particle(q0_fnL, npart);
XL = sim_particle(nstep, K_fn, dK_fn, X0L);

X0C = init_particle(q0_fnC, npart);
XC = sim_particle(nstep, K_fn, dK_fn, X0C);

X0R = init_particle(q0_fnR, npart);
XR = sim_particle(nstep, K_fn, dK_fn, X0R);


% Part 1: Infer only K0
X = XL;
q0_fn = q0_fnL;

K_arr = 1:0.25:3;
loglik_arr = zeros(size(K_arr));
for K_ind = 1:length(K_arr(:))
    K = K_arr(K_ind);
    K_param = @(x) (sin(2*pi*x) + K)*Kscale;
    
    loglik_arr(K_ind) = eval_loglikelihood(X, K_param, q0_fn);
end

% Locate MAP
MAP_ind = find(loglik_arr == max(loglik_arr(:)));
K0_MAP = K_arr(MAP_ind);
K_MAP = @(x) (sin(2*pi*x) + K0_MAP)*Kscale;
[loglik_MAP, q_MAP, xmid] = eval_loglikelihood(X, K_MAP, q0_fn);

%% Result visualisation
figure
subplot(2,2,1)
plot(K_arr, loglik_arr)
title("Likelihood vs K0")
subplot(2,2,2)
plot(xmid, q_MAP)
title("q_{MAP}")
subplot(2,2,3)
histogram(X, 40)
title("Particle positions at time T")
subplot(2,2,4)
plot(xmid, K_MAP(xmid), '--')
hold on
plot(xmid, K_fn(xmid))
title("Inferred and Exact K")

%% Sampling
nDOF = 16;
Data = [XL, XC, XR];
theta_Init = 2*Kscale*ones([1, nDOF]);
theta_SSD = theta_Init*0.05;
Nsteps = 50000;

% Load MCMC
addpath('/home/s1046972/opt/Bayesian_SDE_Inference/Two_Dimensional/MCMC');
compute_logLik = @(Data, theta_Canon) compute_logLik_q0(Data, theta_Canon, q0_fnL, q0_fnC, q0_fnR);

[theta_store, logPost_ts, acceptance_ratio] = mh_Gibbs(Data, compute_logLik, @compute_logPrior, theta_Init, theta_SSD, Nsteps);

%%
K_mean = mean(squeeze(theta_store(1, :, Nsteps*0.8:end)), 2);

xDOF = linspace(0, 1, nDOF+1);
xDOF = 0.5*(xDOF(1:end-1) + xDOF(2:end));

figure(2)
plot(xDOF, K_mean, "--")
hold on
plot(xDOF, K_fn(xDOF))


%%
function logLik = compute_logLik_q0(Data, theta_Canon, q0_fnL, q0_fnC, q0_fnR)
    XL = Data(:, 1);
    XC = Data(:, 2);
    XR = Data(:, 3);
    K_dof = theta_Canon;
    K_param = @(x) K_disc(x, K_dof);
    
    logLik_L = eval_loglikelihood(XL, K_param, q0_fnL);
%     logLik_C = eval_loglikelihood(XC, K_param, q0_fnC);
    logLik_R = eval_loglikelihood(XR, K_param, q0_fnR);
    
%     logLik = logLik_L + logLik_C + logLik_R;
    logLik = logLik_L + logLik_R;
end

function logPrior = compute_logPrior(theta_Canon)
    K_dof = theta_Canon;
    
    logPrior = 0;
    if any(K_dof < 0)
        logPrior = -inf;
    end  
end

function K = K_disc(x, K_dof)
gl = 1;
K_dof_gl = [K_dof(end), K_dof, K_dof(1)];

K = zeros(size(x));
xl = x.*length(K_dof);

for i = 1:length(xl)
    xi = xl(i);
    if mod(xl(i), 1) == 0
        K(i) = 0.5*(K_dof_gl(xi+gl) + K_dof_gl(xi+1+gl));
    else
        K(i) = K_dof(ceil(xi));
    end
end

end

%% Evaluate Likelihood
function [loglik, q, xmid] = eval_loglikelihood(X, K_fn, q0_fn)
n = 64;
nstep = 64*n;

[q, xmid] = FVD_solver(n, nstep, K_fn, q0_fn);
loglik = sum(log(intpl_likelihood(X, q)));

end

function X_intpl = intpl_likelihood(X, q)
X_intpl = zeros(size(X));
gl = 1;
n = length(q)-2*gl;  % Two ghost cells

for i = 1:length(X(:))
    % Linear interpolation
    ind = floor(X(i)*n);  % Relative to edges
    alpha = X(i)*n - ind;
    if alpha > 0.5
        beta = alpha - 0.5;  % Shift: Relative to mid-points
        q_eval = (1-beta)*q(gl + ind+1) + beta*q(gl + ind+2);
    else
        beta = 0.5+alpha;
        q_eval = (1-beta)*q(gl + ind) + beta*q(gl + ind+1);
    end
    
    X_intpl(i) = q_eval;
end
end

%% Particles solvers
function X0 = init_particle(q0_fn, npart)
X0 = [];
% Rejection sampling
M = 10;

while (size(X0, 1) < npart)
    Y = rand;
    alpha = q0_fn(Y)/M;
    
    if rand < alpha
        X0 = [X0; Y];
    end
end
end

function X = sim_particle(nstep, K_fn, dK_fn, X0)
imposeBC = @(X) X - floor(X);   % Periodic on [0, 1]

T = 1;
dt = T/nstep;

X = X0;
for k = 1:nstep
    dX = dt* dK_fn(X) + sqrt(2*K_fn(X)*dt) .* randn(size(X));
    X = X + dX;
end

X = imposeBC(X);
end

%% FVD solver
function [q, xmid] = FVD_solver(n, nstep, K_fn, q0_fn)
q_BC = @(q) [q(end-1), q(2:end-1), q(2)];
dq_BC = @(dq) [dq(end), dq, dq(1)];

L = 1;
T = 1;
dx = L/n;
dt = T/nstep;

xedge = linspace(0, L, n+1);

% % Append ghost cells
xmid = 0.5*(xedge(1:end-1) + xedge(2:end));
xmid = [2*xmid(1)-xmid(2), xmid, 2*xmid(end)-xmid(end-1)];

Kedge = K_fn(xedge);
q = q0_fn(xmid);
q = q/(sum(q)*dx);

q = q_BC(q);

for k = 1:nstep
    dq = dt/(dx*dx) * (Kedge(1:end-1).*q(1:end-2) ...
        - (Kedge(2:end)+Kedge(1:end-1)) .* q(2:end-1) ...
        + Kedge(2:end) .* q(3:end));
    dq = dq_BC(dq);
    
    q = q + dq;
end

end

