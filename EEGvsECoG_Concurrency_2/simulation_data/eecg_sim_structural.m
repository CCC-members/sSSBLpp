function [sim_data] = eecg_sim_structural(sim_data,protocol_path)
%% Loading files
leadfield_data1_file = fullfile(protocol_path,'data','MonkeyECoG','Implantation','headmodel_surf_openmeeg_short.mat');
leadfield_data1 = load(leadfield_data1_file);

leadfield_data2_file = fullfile(protocol_path,'data','MonkeyEEG','@intra','headmodel_surf_openmeeg_short.mat');
leadfield_data2 = load(leadfield_data2_file);

data1_channel_file = fullfile(protocol_path,'data','MonkeyECoG','Implantation','channel.mat');
data1_channel = load(data1_channel_file);

data2_channel_file = fullfile(protocol_path,'data','MonkeyEEG','@intra','channel.mat');
data2_channel = load(data2_channel_file);

surface_file = fullfile(protocol_path,'anat','MonkeyEEG','tess_Cortex-mid_Su_gii_6000V.mat');
surface = load(surface_file);

surface_file_short = fullfile(protocol_path,'anat','MonkeyEEG','tess_Cortex-mid_Su_gii_6000V.mat');
surface_short = load(surface_file_short);

%% Split cortical surface into Left and Right Hemispheres
[indL,indR,VerticesL,VerticesR,FacesL,FacesR,Vertices] = split_hemispheres_monkey(surface);
ind3DL                = [3*indL (3*indL - 1) (3*indL - 2)];
ind3DL                = ind3DL';
ind3DL                = ind3DL(:);
ind3DR                = [3*indR (3*indR - 1) (3*indR - 2)];
ind3DR                = ind3DR';
ind3DR                = ind3DR(:);
surface.VerticesL     = VerticesL;
surface.FacesL        = FacesL;
surface.indL          = indL;
surface.ind3DL        = ind3DL;
surface.VerticesR     = VerticesR;
surface.FacesR        = FacesR;
surface.indR          = indR;
surface.ind3DR        = ind3DR;
surface.Vertices      = Vertices;

%% Surface curvature compensator
disp('-->> Creating curvature compensator');
aSulc                 = 5; % baseline of sulci curvature factor
aGiri                 = 5; % baseline of giri curvature factor
bSulc                 = 3; % scale of sulci curvature factor
bGiri                 = 3; % scale of giri curvature factor

Curv                  = surface.Curvature;
Sulc                  = surface.SulciMap;
Curv                  = abs(Curv);
CurvSulc              = zeros(length(Curv),1);
CurvGiri              = zeros(length(Curv),1);
CurvSulc(Sulc == 1)   = aSulc + bSulc.*Curv(Sulc == 1);
CurvSulc(Sulc == 0)   = 1;
CurvGiri(Sulc == 0)   = aGiri + bGiri.*Curv(Sulc == 0);
CurvGiri(Sulc == 1)   = 1;

Sulc3D                = zeros(1,3*length(Sulc));
CurvSulc3D            = zeros(1,3*length(Curv));
CurvGiri3D            = zeros(1,3*length(Curv));

node3 = 1;
for node = 1:length(Curv)
    CurvSulc3D([node3 node3+1 node3+2]) = repmat(CurvSulc(node),1,3);
    CurvGiri3D([node3 node3+1 node3+2]) = repmat(CurvGiri(node),1,3);
    Sulc3D([node3 node3+1 node3+2])     = repmat(Sulc(node),1,3);
    node3                               = node3 + 3;
end

% Compensated Lead Fields
disp('-->> Creating field information');
L_data13D             = leadfield_data1.Gain;
rej_chan              = [22 25 94 101 102 114];
L_data13D(rej_chan,:) = [];
L_data23D             = leadfield_data2.Gain;

GridOrient            = leadfield_data1.GridOrient;
GridAtlas             = leadfield_data1.GridAtlas;

L_data1               = bst_gain_orient(L_data13D, GridOrient,GridAtlas);
L_data2               = bst_gain_orient(L_data23D, GridOrient,GridAtlas);

clearvars leadfield_data1 leadfield_data2

L_data13Dsulc         = L_data13D.*repmat(CurvSulc3D,size(L_data13D,1),1);
L_data23Dsulc         = L_data23D.*repmat(CurvSulc3D,size(L_data23D,1),1);
L_data13Dgiri         = L_data13D.*repmat(CurvGiri3D,size(L_data13D,1),1);
L_data23Dgiri         = L_data23D.*repmat(CurvGiri3D,size(L_data23D,1),1);

L_data1sulc           = bst_gain_orient(L_data13Dsulc,GridOrient,GridAtlas);
L_data2sulc           = bst_gain_orient(L_data23Dsulc,GridOrient,GridAtlas);
L_data1giri           = bst_gain_orient(L_data13Dgiri,GridOrient,GridAtlas);
L_data2giri           = bst_gain_orient(L_data23Dgiri,GridOrient,GridAtlas);

%% Parcells
disp('-->> Creating parcel smoother');
parcellation_none     = cell(length(L_data1),1);
for area = 1:length(L_data1)
    parcellation_none{area}                 = area;
end

parcellation_none3D   = cell(length(L_data1),1);
for area = 1:length(L_data1)
    q0                                     = 3*(area-1);
    parcellation_none3D{area}              = [q0+1;q0+2;q0+3];
end

Atlas                 = surface.rois;
parcellation          = cell(2*length(Atlas),1);
for area = 1:length(Atlas)
    parcellation{area}                     = intersect(Atlas{area},indL,'stable');
    parcellation{area+length(Atlas)}       = intersect(Atlas{area},indR,'stable');
end

parcellation3D        = cell(2*length(Atlas),1);
for area = 1:length(Atlas)
    tmpL                                   = intersect(Atlas{area},indL,'stable');
    for node = 1:length(tmpL)
        q0                                 = 3*(tmpL(node)-1);
        tmp_parcellation3D                 = [q0+1;q0+2;q0+3];
        parcellation3D{area}               = cat(1,parcellation3D{area},tmp_parcellation3D );
    end
    tmpR                                   = intersect(Atlas{area},indR,'stable');
    for node = 1:length(tmpR)
        q0                                 = 3*(tmpR(node)-1);
        tmp_parcellation3D                 = [q0+1;q0+2;q0+3];
        parcellation3D{area+length(Atlas)} = cat(1,parcellation3D{area+length(Atlas)},tmp_parcellation3D );
    end
end

%% Laplacian and Normals
disp('-->> Creating Laplacian&Normals');
Faces                = surface.Faces;
[D,D3D]              = graph_laplacian(Faces);
Dinv                 = speye(length(D))/D;
Dinv                 = (Dinv + Dinv)/2;
D3Dinv               = speye(length(D3D))/D3D;
D3Dinv               = (D3Dinv + D3Dinv)/2;
Ninv                 = blk_diag(GridOrient',1);
N                    = Ninv';
DN                   = D*N;
DNinv                = Ninv*Dinv;
nEigs                = 100;
FEM                  = firstOrderFEM(VerticesL,FacesL);
[evecs,evals]        = eigs(FEM.laplacian,FEM.vtxInnerProds,nEigs,'sm');
evalsd               = diag(evals);
indemd               = find(abs(evalsd)<10^(-12)*abs(evalsd(end)));
evecs(:,indemd)      = [];
stremd               = precomputeEarthMoversADMM(VerticesL,FacesL,evecs);

%% Saving data
sim_data.structural.L_data13D              = L_data13D;
sim_data.structural.L_data13Dgiri          = L_data13Dgiri;
sim_data.structural.L_data13Dsulc          = L_data13Dsulc;
sim_data.structural.L_data23D              = L_data23D;
sim_data.structural.L_data23Dgiri          = L_data23Dgiri;
sim_data.structural.L_data23Dsulc          = L_data23Dsulc;
sim_data.structural.L_data1                = L_data1;
sim_data.structural.L_data1giri            = L_data1giri;
sim_data.structural.L_data1sulc            = L_data1sulc;
sim_data.structural.L_data2                = L_data2;
sim_data.structural.L_data2giri            = L_data2giri;
sim_data.structural.L_data2sulc            = L_data2sulc;
sim_data.structural.surface                = surface;
sim_data.structural.surface_short          = surface_short;
sim_data.structural.data1_channel          = data1_channel;
sim_data.structural.rej_chan               = rej_chan;
sim_data.structural.data2_channel          = data2_channel;
sim_data.structural.parcellation_none      = parcellation_none;
sim_data.structural.parcellation_none3D    = parcellation_none3D;
sim_data.structural.parcellation           = parcellation;
sim_data.structural.parcellation3D         = parcellation3D;
sim_data.structural.GridOrient             = GridOrient;
sim_data.structural.GridAtlas              = GridAtlas;
sim_data.structural.Sulc                   = Sulc;
sim_data.structural.CurvSulc               = CurvSulc;
sim_data.structural.CurvGiri               = CurvGiri;
sim_data.structural.DN                     = DN;
sim_data.structural.DNinv                  = DNinv;
sim_data.structural.D                      = D;
sim_data.structural.Dinv                   = Dinv;
sim_data.structural.D3D                    = D3D;
sim_data.structural.D3Dinv                 = D3Dinv;
sim_data.structural.N                      = N;
sim_data.structural.Ninv                   = Ninv;
sim_data.structural.stremd                 = stremd;
