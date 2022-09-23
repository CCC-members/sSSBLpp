clear all;
close all;
tic
addpath('common_functions');
addpath('simulation_data');
addpath('concurrency_evaluation');
addpath(genpath('\mnt\Cloud\OneDrive\Tools\brainstorm3'));
addpath('D:\Tools\fieldtrip-master');
ft_defaults

%% Su EEG/ECoG preprocessed file
preproced_data_path = 'E:\OneDrive - Neuroinformatics Collaboratory\Papers\sSSBLpp\MonkeyData\cleaned\002';

%% Brainstorm protocol files for EEG/ECoG Lead Fields and EEG/ECoG channel info
protocol_path = 'E:\OneDrive - Neuroinformatics Collaboratory\Papers\sSSBLpp\BS_protocol_EEGvsECoG';

%% Parameters
method_label  = {'sSSBL' 'sSSBLparcelled' 'sSSBLlaplacian' 'sSSBL2Disotropic' 'sSSBL3Disotropic' 'sSSBLunwrapped' 'sSSBL++' 'eLORETA' 'LCMV'};
% sSSBL with priors of smoothness in 2D cartesian space of samples and band frequencies
% sSSBLparcelled with priors of smoothness in 3D cartesian space of samples, band frequencies and parcells
% sSSBL3Dfield with priors of rotational invariance of 3D source fields 
% sSSBLunwrapped with priors of compensation for surface curvature
% sSSBLplus with all prior information

bands         = [0.1 3; 5.1 7; 9.1 14; 16.1 31; 33.1 90];
band_label    = {'delta' 'theta' 'alpha' 'beta' 'gamma'};
% bands         = [0.1 4; 8 16];
% band_label    = {'low' 'high'};

%% Prepare data for ECG/EEG concurrency evaluation
sim_data      = struct;
[sim_data]    = eecg_sim_structural(sim_data,protocol_path);
[sim_data]    = eecg_sim_functional(sim_data,preproced_data_path);
[sim_data]    = eecg_sim_landmark(sim_data,bands);

%% Concurrency evaluation
disp('-->> performing concurrency evaluation');
% outputs: 
% J     - 1xlength(bands) cell array with MEG and EEG source activity vectors concatenated along the second dimension
% stat  - 1xlength(bands) cell array with MEG and EEG source statistic vectors concatenated along the second dimension   
% indms - 1xlength(bands) cell array with MEG and EEG active indexes in concatenated in 1x2 cell array 

%% sSSBL++
disp('-->> Spectral SSBL with all prior information');
ismethod  = 1;
isparcel  = 1; % 0 (no smoothness) 1 (parcel smoothness)
isneigh   = 1; % 0 (no neighbor structure) 1 (Laplacian neighbor structure)
isfield   = 2; % 1 (projected Lead Field) 2 (3D Lead Field with 2D isotropic rotational invariance) 3 (3D Lead Field with 3D isotropic rotational invariance)

% giri compensation
iscurv    = 1; % 0 (no compensation) 1 (giri curvature compensation) -1 (sulci curvature compensation)
[J_giri,stat_giri,indms_giri,T_giri,data_giri] = concurrency(sim_data,bands,ismethod,isparcel,isneigh,isfield,iscurv);

% sulci compensation
iscurv    = -1; % 0 (no compensation) 1 (giri curvature compensation) -1 (sulci curvature compensation)
[J_sulc,stat_sulc,indms_sulc,T_sulc,data_sulc] = concurrency(sim_data,bands,ismethod,isparcel,isneigh,isfield,iscurv);

J_sSSBLpp          = {J_giri J_sulc};
stat_sSSBLpp       = {stat_giri stat_sulc};
indms_sSSBLpp      = {indms_giri indms_sulc};
T_sSSBLpp          = {T_giri T_sulc};
data_sSSBLpp       = {data_giri data_sulc};

%% eLORETA
disp('-->> Spectral eLORETA');
ismethod  = 2;
isparcel  = 0; % 0 (no smoothness) 1 (parcel smoothness)
isneigh   = 0; % 0 (no neighbor structure) 1 (Laplacian neighbor structure)
isfield   = 1; % 1 (projected Lead Field) 2 (3D Lead Field with 2D isotropic rotational invariance) 3 (3D Lead Field with 3D isotropic rotational invariance)
iscurv    = 0; % 0 (no compensation) 1 (giri curvature compensation) -1 (sulci curvature compensation)

[J,stat,indms,T,data] = concurrency(sim_data,bands,ismethod,isparcel,isneigh,isfield,iscurv);
J_eLORETA          = J;
stat_eLORETA       = stat;
indms_eLORETA      = indms;
T_eLORETA          = T;
data_eLORETA       = data;

%% LCMV
disp('-->> Spectral LCMV');
ismethod  = 3;
isparcel  = 0; % 0 (no smoothness) 1 (parcel smoothness)
isneigh   = 0; % 0 (no neighbor structure) 1 (Laplacian neighbor structure)
isfield   = 1; % 1 (projected Lead Field) 2 (3D Lead Field with 2D isotropic rotational invariance) 3 (3D Lead Field with 3D isotropic rotational invariance)
iscurv    = 0; % 0 (no compensation) 1 (giri curvature compensation) -1 (sulci curvature compensation)
 
[J,stat,indms,T,data] = concurrency(sim_data,bands,ismethod,isparcel,isneigh,isfield,iscurv);
J_LCMV          = J; 
stat_LCMV       = stat;
indms_LCMV      = indms;
T_LCMV          = T;
data_LCMV       = data;
%% Band filtered time series
restoredefaultpath
addpath('common_functions');
addpath('simulation_data');
addpath('concurrency_evaluation');
[sim_data] = filt_tseries(sim_data,bands);
addpath(genpath('D:\Tools\brainstorm3'));
addpath('D:\Tools\fieldtrip-master');
ft_defaults

%% Performance
%% sSSBLpp
[space_freq_sSSBLpp,time_freq_sSSBLpp]                     = performance_concurrency(sim_data,J_sSSBLpp,T_sSSBLpp,bands,band_label,2,1,method_label{7});
figure_concurrency(sim_data,data_sSSBLpp,J_sSSBLpp,space_freq_sSSBLpp,time_freq_sSSBLpp,band_label,method_label{7},1);
figure_blur(sim_data,space_freq_sSSBLpp,band_label,method_label{7});

%% eLORETA
[space_freq_eLORETA,time_freq_eLORETA]                     = performance_concurrency(sim_data,J_eLORETA,T_eLORETA,bands,band_label,1,0,method_label{8});
figure_concurrency(sim_data,data_eLORETA,J_eLORETA,space_freq_eLORETA,time_freq_eLORETA,band_label,method_label{8},0);
figure_blur(sim_data,space_freq_eLORETA,band_label,method_label{8});

%% LCMV
[space_freq_LCMV,time_freq_LCMV]                           = performance_concurrency(sim_data,J_LCMV,T_LCMV,bands,band_label,1,0,method_label{9});
figure_concurrency(sim_data,data_LCMV,J_LCMV,space_freq_LCMV,time_freq_LCMV,band_label,method_label{9},0);
figure_blur(sim_data,space_freq_LCMV,band_label,method_label{9});
