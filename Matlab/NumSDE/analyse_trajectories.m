addpath('/home/s1046972/opt/Bayesian_SDE_Inference/Two_Dimensional/')
addpath('/home/s1046972/opt/Bayesian_SDE_Inference/Two_Dimensional/Scripts/')
addpath('/home/s1046972/opt/Bayesian_SDE_Inference/Two_Dimensional/Utilities/')
addpath('/home/s1046972/opt/Bayesian_SDE_Inference/Two_Dimensional/StructConversion/')
addpath('/home/s1046972/opt/Bayesian_SDE_Inference/Two_Dimensional/PlotUtilities//')

clear all
close all
%%
save_data = 0

kappa_Profile = 'TG_iso';
kappa_Profile = 'sinusoidal';
%kappa_Profile = 'const';
% veloc_Profile = 'zero';
veloc_Profile = 'TTG';
veloc_Profile = 'TTG100';
RunProfile = ['K_', kappa_Profile];
RunProfile = [veloc_Profile, '_', kappa_Profile];

Nx_bin_ARG = 16;
SamplingInterval_days = 30;
DS_rate = 1;
InfScheme = 'test';
Nsteps_pcv = 0;

run('Diffusivity_inference_from_MCMC_PreProcess');


%% Visualise RemainRatio at different SIv
L_km = 1000;
grid_vis = setup_grid_vis(Mesh, RectMesh_Param, L_km);

SamplingInterval_days_List_vis = [SamplingInterval_days];
SamplingInterval_days_List = SamplingInterval_days_List_vis;


RatioRemainNeigh_List = RatioRemainNeigh;
RatioRemain_List = RatioRemain;

%%
% ------- Remain Ratio -------- %
[Nrow, Ncol] = size(SamplingInterval_days_List);
ax_cur_sp = cell(Nrow, Ncol);
Data_ij_vis_List = cell(Nrow, Ncol);

for k = 1:length(SamplingInterval_days_List(:))
    SamplingInterval_vis = SamplingInterval_days_List(k);
    
    % Reorganise the data
    RatioRemainNeigh_vis = convert_bink_to_binij(RatioRemain_List(:, :, k), RectMesh_Param);
    
    Data_ij_vis_List{k} = RatioRemainNeigh_vis;
end

% Log10 scale
row_CLim = {[0, 1], [0, 1]};

% To be appead later
YLabelString_List = {'Fraction'};
TitleString_List =  {''};

FigObj = figure('Name', 'Fields_RemainRatio', 'NumberTitle','off');
set(FigObj, 'Units', 'inches');
set(FigObj, 'Position', [5 5 7.25 7.25/2.5])


ax_sp = imagesc_RClabel(grid_vis, Data_ij_vis_List, row_CLim, YLabelString_List, TitleString_List);

for k = 1:length(SamplingInterval_days_List(:))
    ax_cur = ax_sp{k};
    SamplingInterval_vis = SamplingInterval_days_List(k);
    
    if SamplingInterval_vis > 1
        title_str = sprintf('$s=%d$ days', SamplingInterval_vis);
    else
        title_str = sprintf('$s=%d$ day', SamplingInterval_vis);
    end
    title(ax_cur, title_str);
end


% Remove Ticks
run('Script_RemoveTicks');

% Ad-hoc: systematically enlarge the plots
run('Script_Enlargesubaxis');