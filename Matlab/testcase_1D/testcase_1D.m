% clear all
% nIND = 2;
% nDOF = 16;
% PriorID = 18;
% DS_rate = 4*1024;
% % DS_rate = 1;

%% Inference
% Truth
Kscale = 0.005;
K0 = 2;

% Sinusodial
K_fn = @(x) (sin(2*pi*x) + K0)*Kscale;
dK_fn = @(x) 2*pi*cos(2*pi*x)*Kscale;
d2K_fn = @(x) -2*pi*2*pi*sin(2*pi*x)*Kscale;

% Mixed Gaussian
K0 = 1;
dnormpdf = @(x,mu,sigma) -(x-mu)/(sqrt(2*pi)*sigma^3) .* exp(-0.5*((x-mu)/sigma).^2);
K_profile = @(x) (normpdf(x,1/3,0.1) + normpdf(x,2/3,0.1))/2;
dK_profile = @(x) (dnormpdf(x,1/3,0.1) + dnormpdf(x,2/3,0.1))/2;
K_fn = @(x)  (K0 + (K_profile(x)-K_profile(0))/(K_profile(0.5)-K_profile(0)))*Kscale;
dK_fn = @(x)  dK_profile(x)/(K_profile(0.5)-K_profile(0))*Kscale;

% Initialise indicator functions
nIND_d = 128;  %_d: data
q0_fn_List = cell(nIND_d, 1);
for IND = 1:nIND_d
    q0_fn_List{IND} = @(x) nIND_d.*(x <= IND/nIND_d).*(x > (IND-1)/nIND_d);
end

filename = sprintf('Gauss_NPart24_2T.mat', nIND_d);
if exist(filename, 'file') == 2
    disp("Loading simulated particles")
    load(filename)
else
    disp("Simulating particles")
    % Simulate particles
    nstep = 2^12;
    npart = 2^24/nIND_d;
    
    X2_List = cell(nIND_d, 1);
    X1_List = cell(nIND_d, 1);
    X0_List = cell(nIND_d, 1);
    
    parfor IND = 1:nIND_d
        rng(127*IND)
        X0 = init_particle(q0_fn_List{IND}, npart)
        X1 = sim_particle(nstep, K_fn, dK_fn, X0);
        X2 = sim_particle(nstep, K_fn, dK_fn, X1);
        
        X2_List{IND} = X2;
        X1_List{IND} = X1;
        X0_List{IND} = X0;
        disp(['Finisihed IND=', num2str(IND)]);
    end
    
    save(filename, 'X2_List','X0_List','X1_List', 'K_fn','dK_fn','K0', 'Kscale')
end

%% Re-grouping indicator functions and data
%nIND = 16;
%nIND = 16;
X_List = X1_List;

% AD-HOC: amendment for T=2
% X_List = X2_List;
% END AD-HOC: amendment

assert(mod(length(X_List), nIND)==0)
Group = length(X_List)/nIND;

q0_fnSet = cell(nIND, 1);
DataX = cell(nIND, 1);
Data0 = cell(nIND, 1);
for IND = 1:nIND
    q0_fnSet{IND} = @(x) nIND.*(x <= IND/nIND).*(x > (IND-1)/nIND);
    
    Data0{IND} = X0_List{(IND-1)*Group + 1};
    DataX{IND} = X_List{(IND-1)*Group + 1};
    for group_id = 2:Group
        DataX{IND} = cat(1, DataX{IND}, X_List{(IND-1)*Group + group_id});
        Data0{IND} = cat(1, Data0{IND}, X0_List{(IND-1)*Group + group_id});
    end
end

for DS_rate = [4096, 1]
    for PriorID = [18, 0, 17]
    %% Downsample data
    DataDS = cell(size(DataX));
    
    for IND = 1:nIND
        X_full = DataX{IND};
        DataDS{IND} = X_full(1:DS_rate:end);
        
        X0_full = Data0{IND};
        Data0{IND} = X0_full(1:DS_rate:end);
    end
    
    %% Recast data
    Data_QW = cell(size(DataDS));
    Data_ind0 = cell(size(DataDS));
    Data_gamma0 = cell(size(DataDS));
    
    for IND = 1:nIND
        X = DataDS{IND};
        
        nx = 80;
        % Set up ghost layer
        gl = 1;
        QW = zeros(1, nx+ 2*gl);  % q weight (including ghost layer)
        ind0 = zeros(size(X));  % Where (1-beta) = gamma is
        gamma0 = zeros(size(X));
        
        for i = 1:length(X(:))
            % Linear interpolation
            ind = floor(X(i)*nx);  % Relative to edges
            alpha = X(i)*nx - ind;
            if alpha > 0.5
                beta = alpha - 0.5;  % Shift: Relative to mid-points
                QW(gl+ind+1) = QW(gl+ind+1) + (1-beta);
                QW(gl+ind+2) = QW(gl+ind+2) + beta;
                
                ind0(i) = gl+ind+1;
                gamma0(i) = (1-beta);
            else
                beta = 0.5+alpha;
                QW(gl+ind) = QW(gl+ind) + (1-beta);
                QW(gl+ind+1) = QW(gl+ind+1) + beta;
                
                ind0(i) = gl+ind;
                gamma0(i) = (1-beta);
            end
        end
        
        Data_QW{IND} = QW;  % Only for Q is linear
        Data_ind0{IND} = ind0;
        Data_gamma0{IND} = gamma0;
    end
    
    Data_Bin = {Data_ind0, Data_gamma0};
    
    %% Unittest for interpolation operator
    % IND = 1
    % q = 1:80;
    % q = [q(end), q, q(1)];
    % tic
    % for k = 1:1
    % % intpl1 = sum(log(intpl_likelihood(Data{IND}, q)));
    % % intpl2 = sum(Data_QW{IND}.*log(q))  % if not for log(q) then perfect
    % intpl3 = sum(log(q(Data_ind0{IND})'.*Data_gamma0{IND} + q(Data_ind0{IND}+1)'.*(1-Data_gamma0{IND})));
    % end
    % q = toc;
    % q
    % stop
    
    %% Sampling
    %nDOF = 32;
    xDOF = linspace(0, 1, nDOF+1);
    xDOF = 0.5*(xDOF(1:end-1) + xDOF(2:end));
    
    theta_Init = Kscale*ones([1, nDOF]);
    theta_SSD = 0.20*Kscale*ones([1, nDOF]);
    
    Nsteps_pv = 1000;
    Nsteps = Nsteps_pv*nDOF;
    
    % %dx = 1/nDOF;
    % LapK_exact = d2K_fn(xDOF(2:end-1));
    % LapK = (theta_Init(3:end)+theta_Init(1:end-2)-2*theta_Init(2:end-1))*(nDOF^2);
    
    % Load MCMC
    addpath('/home/s1046972/opt/Bayesian_SDE_Inference/Two_Dimensional/MCMC');
    
    
    % Data = DataX;
    % compute_logLik = @(Data, theta_Canon) compute_logLik_q0(Data, theta_Canon, q0_fnSet);
    
    % Data = Data_QW;
    % compute_logLik = @(Data, theta_Canon) compute_logLik_qW(Data, theta_Canon, q0_fnSet);
    
    Data = Data_Bin;
    compute_logLik = @(Data, theta_Canon) compute_logLik_Bin(Data, theta_Canon, q0_fnSet);
    
    compute_logPrior = @(theta_Canon) compute_logPrior_ID(theta_Canon, PriorID);
    
    % Sampling process
    for k = 1:5
        filepath = "./output/";
        filename = sprintf("MCMC_Gauss_DOF%i_IND%i_dx80_DS%i_Prior%i_%i.mat", nDOF, nIND, DS_rate, PriorID, k);
        filename = strcat(filepath, filename)
        
        rng(1039*k)
        
        if exist(filename, 'file') == 2
            load(filename)
            
            theta_Init = theta_store(1, :, end);
            theta_SSD(acceptance_ratio < 0.15) = 0.5*theta_SSD(acceptance_ratio < 0.15);
            theta_SSD(acceptance_ratio > 0.35) = 1.5*theta_SSD(acceptance_ratio > 0.35);
        else
            [theta_store, logPost_ts, acceptance_ratio] = mh_Gibbs(Data, compute_logLik, compute_logPrior, theta_Init, theta_SSD, Nsteps);
            save(filename, 'theta_store', 'logPost_ts', 'acceptance_ratio', 'theta_SSD')
            
            theta_Init = theta_store(1, :, end);
            theta_SSD(acceptance_ratio < 0.15) = 0.5*theta_SSD(acceptance_ratio < 0.15);
            theta_SSD(acceptance_ratio > 0.35) = 1.5*theta_SSD(acceptance_ratio > 0.35);
        end
        
        if exist('theta_store_ALL') == 0
            theta_store_ALL = theta_store;
            logPost_ts_ALL = logPost_ts;
        else
            theta_store_ALL = cat(3, theta_store_ALL, theta_store);
            logPost_ts_ALL = cat(3, logPost_ts_ALL, logPost_ts);
        end
    end
    end
end

%% Visualisation
addpath('/home/s1046972/opt/Bayesian_SDE_Inference/Two_Dimensional/SampleAnalysis');
addpath('/home/s1046972/opt/Bayesian_SDE_Inference/Two_Dimensional/data_visualisation_nonkernel/toolbox/plottools');
vis = 0;

if vis == 1
    % Part 1: Infer only K0
    X = Data{1};
    q0_fn = q0_fnSet{1};
    
    K_arr = 0:0.1:3;
    loglik_arr = zeros(size(K_arr));
    for K_ind = 1:length(K_arr(:))
        K = K_arr(K_ind);
        K_param = @(x) (K + (K_profile(x)-K_profile(0))/(K_profile(0.5)-K_profile(0)))*Kscale;
        
        loglik_arr(K_ind) = eval_loglikelihood(X, K_param, q0_fn);
    end
    
    % Locate MAP
    MAP_ind = find(loglik_arr == max(loglik_arr(:)));
    K0_MAP = K_arr(MAP_ind);
    K0_MAPfn = @(x) (K0_MAP + (K_profile(x)-K_profile(0))/(K_profile(0.5)-K_profile(0)))*Kscale;
    [loglik_MAP, q_MAP, xmid] = eval_loglikelihood(X, K0_MAPfn, q0_fn);
    
    
    figure
    subplot(2,2,1)
    plot(K_arr, loglik_arr)
    title("Likelihood vs K0")
    subplot(2,2,2)
    plot(xmid, q_MAP)
    title("q_{MAP}")
    subplot(2,2,3)
    NBin = 64;
    Edges = linspace(0, 1, NBin+1);
    hdata = histogram(X, Edges, 'Normalization', 'pdf');
    q_empirical = hdata.Values;
    title("Particle positions at time T")
    subplot(2,2,4)
    plot(xmid, K0_MAPfn(xmid), '--')
    hold on
    plot(xmid, K_fn(xmid))
    title("Inferred and Exact K")
    
    % Part 2: Infer only K0
    
    %assert(Nsteps == size(theta_store,3)-1);
    assert(nDOF == size(theta_store,2));
    
    theta_store = theta_store_ALL;
    logPost_ts = logPost_ts_ALL;
    
    % Histogram of selected cells
    % figure(10)
    % for i = 3:6
    %     subplot(1,4,i-2)
    %     histogram(theta_store(1, i, Nsteps*0.5:end))
    %     hold on
    % end
    
    % Analyse the data
    prctile_p = [0:1:100];
    theta_Survey = survey_theta_store_kernel(theta_store, prctile_p);
    
    pdf_y = squeeze(theta_Survey.theta_pdf_x)';
    pdf = squeeze(theta_Survey.theta_pdf)';
    
    K_mean = theta_Survey.theta_Moments(1,:,1);
    K_MAP = theta_Survey.theta_MAP;
    
    figure(11)
    theta_sample = squeeze(theta_store(1, 1:3, :))';
    plot(theta_sample(1:5:end, :))
    xlabel('Sample')
    legend('Left', 'Of interest', 'Right')
    hold on
    
    %
    % xDOF = linspace(0, 1, nDOF+1);
    % xDOF = 0.5*(xDOF(1:end-1) + xDOF(2:end));
    
    xvis = linspace(0, 1, 1000);
    K_exact = K_fn(xDOF);
    K_MAP_f = @(x) K_disc(x, K_MAP);
    K_mean_f = @(x) K_disc(x, K_mean);
    K_exact_f = @(x) K_disc(x, K_exact);
    
    %
    figure(20)
    ax = gca
    hold(ax, 'on')
    plot_pdf_heatmap(ax, xDOF, pdf_y, pdf)
    h_mean = plot(xvis, K_mean_f(xvis), "--")
    hold on
    h_MAP = plot(xvis, K_MAP_f(xvis), ".")
    h_exact = plot(xvis, K_fn(xvis))
    h_exact_pw = plot(xvis, K_exact_f(xvis), ":");
    
    
    xlabel('x')
    ylabel('\kappa(x)')
    legend([h_mean, h_MAP, h_exact, h_exact_pw], {'Mean', 'MAP', 'Exact', 'Exact PW'})
    
    %
    FigObj = figure('Units', 'inches', 'Position', [5 5 7.25 7.25]);
    ax_cur = gca
    hold(ax_cur, 'on')
    plot_pdf_heatmap(ax_cur, xDOF, pdf_y, pdf)
    h_MAP = plot(ax_cur, xvis, K_MAP_f(xvis), ".")
    h_exact = plot(ax_cur, xvis, K_fn(xvis))
    
    tilte_string = sprintf("bins = %i; DOF = %i", nIND, nDOF )
    title(ax_cur, tilte_string)
    xlabel(ax_cur, '$x$', 'Interpreter', 'latex')
    ylabel(ax_cur, '$\kappa(x)$','Interpreter', 'latex')
    legend([h_MAP, h_exact], {'MAP', 'Exact'})
    
    set(ax_cur.XLabel, 'FontSize', 24, 'Interpreter', 'latex')
    set(ax_cur.YLabel, 'FontSize', 24, 'Interpreter', 'latex')
    set(ax_cur.Title, 'FontSize', 24, 'Interpreter', 'latex')
    set(ax_cur.XAxis, 'FontSize', 18)
    set(ax_cur.YAxis, 'FontSize', 18)
    set(ax_cur,'TickLabelInterpreter', 'latex');
    set(ax_cur.Legend, 'FontSize', 18, 'Interpreter', 'latex')
    
    % Overfitting
    n = NBin;
    q0_fn = q0_fnSet{1};
    q_particles = [q_empirical(end) q_empirical q_empirical(1)];
    
    nstep = 32*n;
    [q_MAP, xmid] = FVD_solver(n, nstep, K_MAP_f, q0_fn);
    
    nstep = 32*n;
    [q_exact, xmid] = FVD_solver(n, nstep, K_fn, q0_fn);
    
    
    FigObj = figure('Units', 'inches', 'Position', [5 5 7.25 7.25*1.25]);
    
    ax_cur = subplot(2,1,1)
    plot(ax_cur, xmid, q_MAP, '--', xmid, q_particles, '.', xmid, q_exact)
    title('q/ particle distribution')
    legend('MAP', 'Particle', 'Exact')
    
    xlabel(ax_cur, '$x$', 'Interpreter', 'latex')
    ylabel(ax_cur, '$q$', 'Interpreter', 'latex')
    set(ax_cur.XLabel, 'FontSize', 24, 'Interpreter', 'latex')
    set(ax_cur.YLabel, 'FontSize', 24, 'Interpreter', 'latex')
    set(ax_cur.Title, 'FontSize', 24, 'Interpreter', 'latex')
    set(ax_cur.XAxis, 'FontSize', 18)
    set(ax_cur.YAxis, 'FontSize', 18)
    set(ax_cur,'TickLabelInterpreter', 'latex');
    set(ax_cur.Legend, 'FontSize', 18, 'Interpreter', 'latex')
    
    ax_cur = subplot(2,1,2)
    plot(ax_cur, xmid, q_MAP-q_exact,'--', xmid, q_particles-q_exact, '-')
    title('$q - q_{exact}$')
    legend('MAP', 'Particle')
    
    xlabel(ax_cur, '$x$', 'Interpreter', 'latex')
    ylabel(ax_cur, '$q - q_{exact}$', 'Interpreter', 'latex')
    set(ax_cur.XLabel, 'FontSize', 24, 'Interpreter', 'latex')
    set(ax_cur.YLabel, 'FontSize', 24, 'Interpreter', 'latex')
    set(ax_cur.Title, 'FontSize', 24, 'Interpreter', 'latex')
    set(ax_cur.XAxis, 'FontSize', 18)
    set(ax_cur.YAxis, 'FontSize', 18)
    set(ax_cur,'TickLabelInterpreter', 'latex');
    set(ax_cur.Legend, 'FontSize', 18, 'Interpreter', 'latex')
    corrcoef(q_MAP-q_exact, q_particles-q_exact)
end

%% Functions
function logLik = compute_logLik_q0(Data, theta_Canon, q0_fnSet)
NDataSet = length(Data);
logLik_Set = zeros(NDataSet, 1);

K_dof = theta_Canon;
K_param = @(x) K_disc(x, K_dof);

parfor DataID = 1:NDataSet
    X = Data{DataID};
    q0_fn = q0_fnSet{DataID};
    logLik_Set(DataID) = eval_loglikelihood(X, K_param, q0_fn);
end

logLik = sum(logLik_Set);
end

function logLik = compute_logLik_qW(Data_QW, theta_Canon, q0_fnSet)
NDataSet = length(Data_QW);
logLik_Set = zeros(NDataSet, 1);

K_dof = theta_Canon;
K_param = @(x) K_disc(x, K_dof);

for DataID = 1:NDataSet
    qW = Data_QW{DataID};
    q0_fn = q0_fnSet{DataID};
    logLik_Set(DataID) = eval_loglikelihood_QW(qW, K_param, q0_fn);
end

logLik = sum(logLik_Set);
end

function logLik = compute_logLik_Bin(Data_Bin, theta_Canon, q0_fnSet)
Data_ind0 = Data_Bin{1};
Data_gamma0 = Data_Bin{2};

NDataSet = length(Data_ind0);
logLik_Set = zeros(NDataSet, 1);

K_dof = theta_Canon;
K_param = @(x) K_disc(x, K_dof);

parfor DataID = 1:NDataSet
    ind0 = Data_ind0{DataID};
    gamma0 = Data_gamma0{DataID};
    
    q0_fn = q0_fnSet{DataID};
    logLik_Set(DataID) = eval_loglikelihood_Bin(ind0, gamma0, K_param, q0_fn);
end

logLik = sum(logLik_Set);
end

function logPrior = compute_logPrior_ID(theta_Canon, PriorID)
K_dof = [theta_Canon, theta_Canon(1)];  % Periodic conditions
nDOF = length(theta_Canon);

logPrior = 0;

if PriorID == 0
    if any(K_dof > 1)
        logPrior = -inf;
    end
elseif PriorID == 17
    dKdx_scale = 0.06*4;
    dx = 1/nDOF;
    GradK = diff(K_dof)/dx;
    logPrior = -sum(GradK.^2)/(2*(dKdx_scale^2));
elseif PriorID == 18
    d2Kdx2_scale = 1*4;
    dx = 1/nDOF;
    LapK = (K_dof(3:end)+K_dof(1:end-2)-2*K_dof(2:end-1))/(dx^2);
    logPrior = -sum(LapK.^2)/(2*(d2Kdx2_scale^2));
end

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

%% FVD solver
function [q, xmid] = FVD_solver(n, nstep, K_fn, q0_fn)
q_BC = @(q) [q(end-1), q(2:end-1), q(2)];
dq_BC = @(dq) [dq(end), dq, dq(1)];

L = 1;
T = 1;
dx = L/n;
dt = T/nstep;

% AD-HOC: amendment for T=2
% T = 2;
% nstep = 2*nstep;
% END AD-HOC: amendment

xedge = linspace(0, L, n+1);

% % Append ghost cells
xmid = 0.5*(xedge(1:end-1) + xedge(2:end));
xmid = [2*xmid(1)-xmid(2), xmid, 2*xmid(end)-xmid(end-1)];

Kedge = K_fn(xedge);
q = q0_fn(xmid);
q = q_BC(q);

% Normalise
q = q/(sum(q(2:end-1))*dx);


% q = q0;

% Forward Euler
for k = 1:nstep
    dq = dt/(dx*dx) * (Kedge(1:end-1).*q(1:end-2) ...
        - (Kedge(2:end)+Kedge(1:end-1)) .* q(2:end-1) ...
        + Kedge(2:end) .* q(3:end));
    dq = dq_BC(dq);
    
    q = q + dq;
end

% Output: included q

%     % Heun's method
% for k = 1:nstep
%     % Evaluate correction q
%     dq = dt/(dx*dx) * (Kedge(1:end-1).*q(1:end-2) ...
%         - (Kedge(2:end)+Kedge(1:end-1)) .* q(2:end-1) ...
%         + Kedge(2:end) .* q(3:end));
%     dq = dq_BC(dq);
%
%     qc = q + dq;
%
%     % Evaulate correctoin dq
%     dqc = dt/(dx*dx) * (Kedge(1:end-1).*qc(1:end-2) ...
%         - (Kedge(2:end)+Kedge(1:end-1)) .* qc(2:end-1) ...
%         + Kedge(2:end) .* qc(3:end));
%     dqc = dq_BC(dqc);
%
%     q = q + 0.5*(dq+dqc);
% end

end



%% Evaluate Likelihood
function [loglik, q, xmid] = eval_loglikelihood_Bin(ind0, gamma0, K_fn, q0_fn)
nx = 80;
% n = 256;
nstep = 32*nx;

[q, xmid] = FVD_solver(nx, nstep, K_fn, q0_fn);
loglik = sum(log(q(ind0)'.*gamma0 + q(ind0+1)'.*(1-gamma0)));

end

function [loglik, q, xmid] = eval_loglikelihood_QW(QW, K_fn, q0_fn)
gl = 1;
n = length(QW)-2*gl;
% n = 256;
nstep = 32*n;

[q, xmid] = FVD_solver(n, nstep, K_fn, q0_fn);
loglik = sum(log(q).*QW);

end

% Old, inefficient method
function [loglik, q, xmid] = eval_loglikelihood(X, K_fn, q0_fn)
nx = 80;
% n = 256;
nstep = 32*nx;

[q, xmid] = FVD_solver(nx, nstep, K_fn, q0_fn);
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
    
    X = imposeBC(X);
end

end

