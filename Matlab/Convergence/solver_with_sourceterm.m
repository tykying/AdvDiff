clear all

m = [];
error = [];

% Adaptive Markov Chain Monte Carlo: Theory and Method Atchad
% An Adaptive Metropolis Algorit Haario
%% dt = C dx
% Lax-Wendoff (temporal_error_LaxWendoff_dx.txt)
nts(           1 ) =        2048 ; error(           1 ) =    2.3354105808755108E-002
nts(           2 ) =        4096 ; error(           2 ) =    6.0544307904864892E-003
nts(           3 ) =        8192 ; error(           3 ) =    1.5237828301798892E-003
nts(           4 ) =       16384 ; error(           4 ) =    3.7994672647586025E-004
nts(           5 ) =       32768 ; error(           5 ) =    9.4248868737989135E-005
nts(           6 ) =       65536 ; error(           6 ) =    2.3323910947910241E-005

 m(1) =    16; nts(1) =  2048; error(1) = 1.0150621962693467E-002
 m(2) =    32; nts(2) =  4096; error(2) = 2.6503883346978827E-003
 m(3) =    64; nts(3) =  8192; error(3) = 6.7604298641262982E-004
 m(4) =   128; nts(4) = 16384; error(4) = 1.7340954961623441E-004
 m(5) =   256; nts(5) = 32768; error(5) = 4.6382300954223657E-005
 m(6) =   512; nts(6) = 65536; error(6) = 1.3936601350631727E-005
 m(7) =  1024; nts(7) =131072; error(7) = 5.1131865209441011E-006


nts_LW = nts;
error_LW = error;


% % MC (temporal_error_MC_dx.txt)
nts(           1 ) =        2048 ; error(           1 ) =    2.2934563130978993E-002
nts(           2 ) =        4096 ; error(           2 ) =    4.5426445050427361E-003
nts(           3 ) =        8192 ; error(           3 ) =    6.8441461963404280E-004
nts(           4 ) =       16384 ; error(           4 ) =    1.6300199072385039E-004
nts(           5 ) =       32768 ; error(           5 ) =    2.7267157288033321E-005
nts(           6 ) =       65536 ; error(           6 ) =    7.8651262620956072E-006

 m(1) =    16; nts(1) =  2048; error(1) = 1.4321477890812081E-002
 m(2) =    32; nts(2) =  4096; error(2) = 4.3610788135775517E-003
 m(3) =    64; nts(3) =  8192; error(3) = 5.4447675033657370E-004
 m(4) =   128; nts(4) = 16384; error(4) = 1.6211866985435639E-004
 m(5) =   256; nts(5) = 32768; error(5) = 1.8206682967721743E-005
 m(6) =   512; nts(6) = 65536; error(6) = 4.8927217999051925E-006
 m(7) =  1024; nts(7) =131072; error(7) = 1.8353616411431197E-006


nts_MC = nts;
error_MC = error;


% % Minmod (temporal_error_MM_dx.txt)
nts(           1 ) =        2048 ; error(           1 ) =    4.6924047429488878E-002
nts(           2 ) =        4096 ; error(           2 ) =    1.4603795905765285E-002
nts(           3 ) =        8192 ; error(           3 ) =    4.0386627222084703E-003
nts(           4 ) =       16384 ; error(           4 ) =    1.0290227962327962E-003
nts(           5 ) =       32768 ; error(           5 ) =    2.5777592878034329E-004
nts(           6 ) =       65536 ; error(           6 ) =    6.4318176701857685E-005

nts_MM = nts;
error_MM = error;

% % Sweby (temporal_error_MM_dx.txt)
 nts(           1 ) =        2048 ; error(           1 ) =    4.5905653529079972E-002
 nts(           2 ) =        4096 ; error(           2 ) =    1.2328063395295605E-002
 nts(           3 ) =        8192 ; error(           3 ) =    3.6329316605049327E-003
 nts(           4 ) =       16384 ; error(           4 ) =    1.0136263252839908E-003
 nts(           5 ) =       32768 ; error(           5 ) =    2.6390448045991737E-004
 nts(           6 ) =       65536 ; error(           6 ) =    6.7355050364254862E-005

nts_SB = nts;
error_SB = error;

% % LW under 2D advection + Strang splitting
 nts(           1 ) =        2048 ; error(           1 ) =    2.3379569003564547E-002
 nts(           2 ) =        4096 ; error(           2 ) =    6.0694405574502090E-003
 nts(           3 ) =        8192 ; error(           3 ) =    1.5314861725730785E-003
 nts(           4 ) =       16384 ; error(           4 ) =    3.8372944499193054E-004
 nts(           5 ) =       32768 ; error(           5 ) =    9.6035577174224414E-005

nts_LW2 = nts;
error_LW2 = error;



%%
% addpath('/home/s1046972/opt/Bayesian_SDE_Inference/Two_Dimensional/data_visualisation_nonkernel/scripts')

nts_List = {nts_LW, nts_MC, nts_MM, nts_SB};
error_List = {error_LW, error_MC, error_MM, error_SB};
label_List = {'Lax-Wendoff', 'MC', 'Minmod', 'Sweby'};

% Verify same dt
nts = nts_List{1};
for k = 1:length(nts_List)
    assert(all(nts == nts_List{k}));
end
m = nts/128;

FigObj = figure('Units', 'inches', 'Position', [5 5 7.25 7.25*1.25]);

for k = 1:length(error_List)
    h = 1./nts;
    error = error_List{k};
    
    e_vis = error(2:end)./error(1:end-1);
    h_vis = h(2:end)./h(1:end-1);
    
    ax = subplot(2,1,1)
    plot(ax, m(2:end), log(e_vis)./log(h_vis))
    hold(ax, 'on')
    
    ax2 = subplot(2,1,2)
    loglog(ax2, h,error)
    hold(ax2, 'on')
end

legend(ax2, label_List, 'Location','southeast');

% Polish the figure
ax_cur = subplot(2,1,1)
ylim(ax_cur, [1.5, 3])
xlabel(ax_cur, '$N=1/h$', 'Interpreter', 'latex')
ylabel(ax_cur, '$\log(\Delta e / \Delta h)$', 'Interpreter', 'latex')
set(ax_cur.XLabel, 'FontSize', 24, 'Interpreter', 'latex')
set(ax_cur.YLabel, 'FontSize', 24, 'Interpreter', 'latex')
set(ax_cur.XAxis, 'FontSize', 18)
set(ax_cur.YAxis, 'FontSize', 18)
set(ax_cur,'TickLabelInterpreter', 'latex');

ax_cur = subplot(2,1,2)
ylim(ax_cur, [1e-6, 1e-1])

xlabel(ax_cur, '$h$', 'Interpreter', 'latex')
ylabel(ax_cur, '$e_h$', 'Interpreter', 'latex')
set(ax_cur.XLabel, 'FontSize', 24, 'Interpreter', 'latex')
set(ax_cur.YLabel, 'FontSize', 24, 'Interpreter', 'latex')
set(ax_cur.XAxis, 'FontSize', 18)
set(ax_cur.YAxis, 'FontSize', 18)
set(ax_cur,'TickLabelInterpreter', 'latex');
set(ax_cur.Legend, 'FontSize', 18, 'Interpreter', 'latex')


%N = kg * m^2/2