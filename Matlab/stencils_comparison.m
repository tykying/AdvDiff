Fig_outputfolder = '/home/s1046972/opt/AdvDiff_FVD/Matlab/FVD/';


stencil_N = [-1/2, 2, 1/2; 2, -8, 2; 1/2, 2, -1/2];
stencil_MPFA = [-1/4, 3/2, 3/4; 3/2, -7, 3/2; 3/4, 3/2, -1/4];

t = 1;
t = 20;

%a = 2; c = 2; b = 4; N=128; k_step = 50; dt = 0.1; %% Observe the instability

param_set = 4;

if param_set == 1
    a = 2; c = 1.5; b = 4; dt = 0.1;  %% Observe the instability in Naive
elseif param_set == 2
    a = 2; c = 1.5; b = 4; dt = 0.075;  %% No instability
elseif param_set == 3
    a = 2; c = 1.5; b = 4; dt = 0.125;  %% Observe the instability in both
elseif param_set == 4
    a = 3; c = -1.5; b = 4; dt = 0.05;  %% Observe the instability in both
end

k_step = t/dt;

K = [a, c; c, b];
eig(K)
assert(all(eig(K)>0))

stencil_N = [-c/2, b, c/2; a, -2*(a+b), a; c/2, b, -c/2];
stencil_MPFA = [c^2/(4*a)+c^2/(4*b)-(1/2)*c, -c^2/(2*b)+(2*a*b-c^2)/(2*a), (1/2)*c+c^2/(4*a)+c^2/(4*b); ...
               -c^2/(2*a)+(2*a*b-c^2)/(2*b), -(2*a*b-c^2)/b-(2*a*b-c^2)/a, -c^2/(2*a)+(2*a*b-c^2)/(2*b); ...
               (1/2)*c+c^2/(4*a)+c^2/(4*b), -c^2/(2*b)+(2*a*b-c^2)/(2*a), c^2/(4*a)+c^2/(4*b)-(1/2)*c];


N = 128;
q = zeros(N, N);
Ncell=0;

IC = 1
if IC == 1
    % Uniform across grid
    q(N/2-Ncell:N/2+1+Ncell, N/2-Ncell:N/2+1+Ncell) = 1/(2*(Ncell+1))^2;
elseif IC == 2
    % Impose initial Gaussian profile
    x = linspace(1,N,N);
    [X, Y] = meshgrid(x,x);
    Sigma = [16, 0; 0, 16];
    
    for i = 1:(size(q,1))
        for j = 1:(size(q,2))
            d = [X(i,j)- (N+1)/2; Y(i,j)- (N+1)/2];
            q(i, j) = 1/sqrt(det(2*pi*Sigma))*exp(-1/2*d'*inv(Sigma)*d);
        end
    end
end

figure(100)
contourf(q)
title('q0')
pbaspect([1,1,1])


q0 = q;

dq_step_MPFA = dq_stencil(q0, stencil_MPFA);
dq_step_N = dq_stencil(q0, stencil_N);

figure(1)
subplot(1,2,1)
contourf(dq_step_N)
title('dq: Naive')
pbaspect([1,1,1])
subplot(1,2,2)
contourf(dq_step_MPFA)
title('dq: MPFA')
pbaspect([1,1,1])

figure(2)
contourf(dq_step_N-dq_step_MPFA)
title('Naive-MPFA')
pbaspect([1,1,1])


% % Forward Euler-method
% q_N = q0;
% q_MPFA = q0;
% for k = 1:k_step
%     q_MPFA = q_MPFA + dt*dq_stencil(q_MPFA, stencil_MPFA);
%     q_N = q_N + dt*dq_stencil(q_N, stencil_N);
% end
% 
% figure(3)
% subplot(1,2,1)
% contourf(q_N)
% title(['q Euler: Naive; dt = ', num2str(dt)])
% pbaspect([1,1,1])
% subplot(1,2,2)
% contourf(q_MPFA)
% title(['q Euler: MPFA; dt = ', num2str(dt)])
% pbaspect([1,1,1])
% 
% figure(4)
% contourf(q_N-q_MPFA)
% title('q Euler: Naive-MPFA')
% pbaspect([1,1,1])


% Heun's method
q_N = q0;
q_MPFA = q0;
for k = 1:k_step
    q_MPFA_h = q_MPFA + dt*dq_stencil(q_MPFA, stencil_MPFA);
    q_MPFA = q_MPFA + 0.5*dt*(dq_stencil(q_MPFA, stencil_MPFA)+dq_stencil(q_MPFA_h, stencil_MPFA));
    
    q_N_h = q_N + dt*dq_stencil(q_N, stencil_N);
    q_N = q_N + 0.5*dt*(dq_stencil(q_N, stencil_N)+dq_stencil(q_N_h, stencil_N));
end


%%
figure(13)
subplot(1,2,1)
contourf(q_N)
title({['q Heun: Naive; dt = ', num2str(dt)]; ['q int = ', num2str(sum(q_N(:)))]})
pbaspect([1,1,1])
subplot(1,2,2)
contourf(q_MPFA)
title({['q Heun: MPFA; dt = ', num2str(dt)]; ['q int = ', num2str(sum(q_MPFA(:)))]})
pbaspect([1,1,1])

figure(14)
contourf(q_N-q_MPFA)
title('q Heun: Naive-MPFA')
pbaspect([1,1,1])


%% Probability
MU = 0.5*[N+1, N+1];
K0 = [2.0000    1.5000;
    1.5000    4.0000];
Sigma_t = 2*K0*t;

rng default  % For reproducibility
R = mvnrnd(MU,Sigma_t,15^2);

FigObj = figure('Name', 'Fields_TrajVis', 'NumberTitle','off');
set(FigObj, 'Units', 'inches');
set(FigObj, 'Position', [5 5 7.25 7.25/2])

subplot(1,2,1)
scatter(R(:,1)*0+MU(1),R(:,2)*0+MU(2),'+')
xlim([1, N]);
ylim([1, N]);
pbaspect([1,1,1])
title({'SDE: 225 particles'; 't = 0'})


subplot(1,2,2)
contourf(q0)
pbaspect([1,1,1])
xlim([1, N]);
ylim([1, N]);
title({'PDE'; 't = 0'})

Fig_outputpath = [ Fig_outputfolder, 'SDE_FP_duality_t0', '.eps'];
%print(FigObj, '-depsc2', Fig_outputpath)


FigObj = figure('Name', 'Fields_TrajVis', 'NumberTitle','off');
set(FigObj, 'Units', 'inches');
set(FigObj, 'Position', [5 5 7.25 7.25/2])

subplot(1,2,1)
scatter(R(:,1),R(:,2),'+')
xlim([1, N]);
ylim([1, N]);
pbaspect([1,1,1])
title('t = 20')


subplot(1,2,2)
contourf(q_MPFA)
pbaspect([1,1,1])
xlim([1, N]);
ylim([1, N]);
title('t = 20')


Fig_outputpath = [ Fig_outputfolder, 'SDE_FP_duality_tf', '.eps'];
%print(FigObj, '-depsc2', Fig_outputpath)

%%
FigObj = figure('Name', 'Fields_TrajVis', 'NumberTitle','off');
set(FigObj, 'Units', 'inches');
set(FigObj, 'Position', [5 5 7.25 7.25/2])

% subplot(2,2,1)
% scatter(R(:,1)*0+MU(1),R(:,2)*0+MU(2),'+')
% xlim([1, N]);
% ylim([1, N]);
% pbaspect([1,1,1])
% title({'SDE: 225 particles'; 't = 0'})
% hold on
% scatter(R(161,1)*0+MU(1),R(161,2)*0+MU(2),'r+')
% 
% subplot(2,2,2)
% contourf(q0)
% pbaspect([1,1,1])
% xlim([1, N]);
% ylim([1, N]);
% title({'PDE'; 't = 0'})

subplot(1,2,1)
scatter(R(:,1),R(:,2),'+')
xlim([1, N]);
ylim([1, N]);
pbaspect([1,1,1])
title('t = 20')
hold on
scatter(R(161,1),R(161,2),'r+')


subplot(1,2,2)
contourf(q_MPFA)
pbaspect([1,1,1])
xlim([1, N]);
ylim([1, N]);
title('t = 20')
hold on
scatter(R(161,1),R(161,2),'r+')

Fig_outputpath = [ Fig_outputfolder, 'SDE_FP_duality_HL_tf', '.eps'];
Fig_outputpath = [ Fig_outputfolder, 'SDE_FP_duality_HL_tf2', '.eps'];
print(FigObj, '-depsc2', Fig_outputpath)


%%
% Applying stencil
function dq_out = dq_stencil(q, stencil_in)
q_copy = q;
dq_out = zeros(size(q));

for i = 2:(size(q,1)-1)
    for j = 2:(size(q,2)-1)
        dq_out(i,j) = sum(sum(stencil_in.*q_copy(i-1:i+1, j+1:-1:j-1)));
    end
end
end
