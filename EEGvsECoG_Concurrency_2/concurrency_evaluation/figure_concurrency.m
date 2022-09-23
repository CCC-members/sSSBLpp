function [mapsfig1,mapsfig2,mapsfig3] = figure_concurrency(sim_data,data0,J0,space_freq,time_freq,band_label,method_label,iscurv)
%% Parameters
%% Colormap
load('mycolor.mat');
load('cmap_pearson');
load('cmap_pval');
surface              = sim_data.structural.surface;
Vertices             = surface.Vertices;
Ngen                 = length(Vertices);
Faces                = surface.Faces;
VertConn             = surface.VertConn;
smoothValue          = 0.66;
SurfSmoothIterations = 10;
Vertices             = tess_smooth(Vertices, smoothValue, SurfSmoothIterations, VertConn, 1);
Sulc                 = surface.SulciMap;
ind                  = surface.indL;
if iscurv == 0
    J       = J0;
else
    for band = 1:length(band_label)
        J{band} = (J0{1}{band} + J0{2}{band})/2;
    end
end
%%
mapsfig1             = figure('Color','w','Name',strcat(method_label,'  ECoG/EEG activation'));
for band = 1:length(band_label)
    %% Source 1
    J1                     = J{band}(:,1);
    map1                   = zeros(Ngen,1);
    map1(ind)              = J1(ind);
    map1                   = map1/max(map1);
    %% Source 2
    J2                     = J{band}(:,2);
    map2                   = zeros(Ngen,1);
    map2(ind)              = J2(ind);
    map2                   = map2/max(map2);
    %% Plot Source1/Source2 (8 subplots)
    %1-----------------------------------------------------------------------
    subplot(10,4,8*(band-1)+1)
    patch('Faces',Faces,'Vertices',Vertices,'FaceVertexCData',...
        log10(map1),'FaceColor','interp','EdgeColor','none','FaceAlpha',.99);
    colormap(gca,cmap);
    axis off
    set(gca,'xcolor','w','ycolor','w','zcolor','w');
    view(180,0);
    text(max(Vertices(:,1)),max(Vertices(:,2)),max(Vertices(:,3)),'L')
    %1-----------------------------------------------------------------------
    if band == 1
        title([method_label ' ' 'ECoG'])
    end
    
    %2-----------------------------------------------------------------------
    subplot(10,4,8*(band-1)+2)
    patch('Faces',Faces,'Vertices',Vertices,'FaceVertexCData',...
        log10(map1),'FaceColor','interp','EdgeColor','none','FaceAlpha',.99);
    colormap(gca,cmap);
    axis off
    set(gca,'xcolor','w','ycolor','w','zcolor','w');
    view(180,90);
    text(max(Vertices(:,1)),min(Vertices(:,2)),max(Vertices(:,3)),'D')
    %2-----------------------------------------------------------------------
    if band == 1
        title([method_label ' ' 'ECoG'])
    end
    
    %3-----------------------------------------------------------------------
    subplot(10,4,8*(band-1)+3)
    patch('Faces',Faces,'Vertices',Vertices,'FaceVertexCData',...
        log10(map2),'FaceColor','interp','EdgeColor','none','FaceAlpha',.99);
    colormap(gca,cmap);
    axis off
    set(gca,'xcolor','w','ycolor','w','zcolor','w');
    view(180,0);
    text(max(Vertices(:,1)),max(Vertices(:,2)),max(Vertices(:,3)),'L')
    %3-----------------------------------------------------------------------
    if band == 1
        title([method_label ' ' 'EEG'])
    end
    
    %4-----------------------------------------------------------------------
    subplot(10,4,8*(band-1)+4)
    patch('Faces',Faces,'Vertices',Vertices,'FaceVertexCData',...
        log10(map2),'FaceColor','interp','EdgeColor','none','FaceAlpha',.99);
    colormap(gca,cmap);
    axis off
    set(gca,'xcolor','w','ycolor','w','zcolor','w');
    view(180,90);
    text(max(Vertices(:,1)),min(Vertices(:,2)),max(Vertices(:,3)),'D')
    %4-----------------------------------------------------------------------
    if band == 1
        title([method_label ' ' 'EEG'])
    end
    
    %5-----------------------------------------------------------------------
    subplot(10,4,8*(band-1)+5)
    patch('Faces',Faces,'Vertices',Vertices,'FaceVertexCData',...
        log10(map1),'FaceColor','interp','EdgeColor','none','FaceAlpha',.99);
    colormap(gca,cmap);
    axis off
    set(gca,'xcolor','w','ycolor','w','zcolor','w');
    view(-90,0);
    text(min(Vertices(:,1)),max(Vertices(:,2)),max(Vertices(:,3)),'P')
    %5------------------------------------------------------------------------
    
    %6------------------------------------------------------------------------
    subplot(10,4,8*(band-1)+6)
    patch('Faces',Faces,'Vertices',Vertices,'FaceVertexCData',...
        log10(map1),'FaceColor','interp','EdgeColor','none','FaceAlpha',.99);
    colormap(gca,cmap);
    axis off
    set(gca,'xcolor','w','ycolor','w','zcolor','w');
    view(180,-90);
    text(max(Vertices(:,1)),max(Vertices(:,2)),min(Vertices(:,3)),'V')
    %6------------------------------------------------------------------------
    
    %7------------------------------------------------------------------------
    subplot(10,4,8*(band-1)+7)
    patch('Faces',Faces,'Vertices',Vertices,'FaceVertexCData',...
        log10(map2),'FaceColor','interp','EdgeColor','none','FaceAlpha',.99);
    colormap(gca,cmap);
    axis off
    set(gca,'xcolor','w','ycolor','w','zcolor','w');
    view(-90,0);
    text(min(Vertices(:,1)),max(Vertices(:,2)),max(Vertices(:,3)),'P')
    %7------------------------------------------------------------------------
    
    %8------------------------------------------------------------------------
    subplot(10,4,8*(band-1)+8)
    patch('Faces',Faces,'Vertices',Vertices,'FaceVertexCData',...
        log10(map2),'FaceColor','interp','EdgeColor','none','FaceAlpha',.99);
    colormap(gca,cmap);
    axis off
    set(gca,'xcolor','w','ycolor','w','zcolor','w');
    view(180,-90);
    text(max(Vertices(:,1)),max(Vertices(:,2)),min(Vertices(:,3)),'V')
    %8------------------------------------------------------------------------
end
file_name = strcat(method_label,'_EECG_activation','.fig');
savefig(mapsfig1,file_name);
%%
mapsfig2             = figure('Color','w','Name',strcat(method_label,'  Temporal P/Corr'));
for band = 1:length(band_label)
    %% Partial Correlations
    pcr                     = time_freq.pcr(:,band);
    map1(ind)               = pcr;
    map1(map1 <= 0)         = -0.05;
    map1(map1 >= 0.5)       = 0.5;
    %% P-value
    pvp                     = time_freq.pvp(:,band);
    map2(ind)               = pvp;
    map2(map2 > 0.05)       = 0.052;
    map2(map1 <= 0)         = 0.052;   
    map1(map2 > 0.05)       = -0.05;
    %% Plot PartialCorrelations/P-value
    %1-----------------------------------------------------------------------
    subplot(10,4,8*(band-1)+1)
    patch('Faces',Faces,'Vertices',Vertices,'FaceVertexCData',...
        map1,'FaceColor','interp','EdgeColor','none','FaceAlpha',.99);
    colormap(gca,cmap_pearson);
    caxis([0 0.5]);
    axis off
    set(gca,'xcolor','w','ycolor','w','zcolor','w');
    view(180,0);
    text(max(Vertices(:,1)),max(Vertices(:,2)),max(Vertices(:,3)),'L')
    %1-----------------------------------------------------------------------
    if band == 1
        title('PCorr(p-val<=0.05*)')
    end
    
    %2-----------------------------------------------------------------------
    subplot(10,4,8*(band-1)+2)
    patch('Faces',Faces,'Vertices',Vertices,'FaceVertexCData',...
        map1,'FaceColor','interp','EdgeColor','none','FaceAlpha',.99);
    colormap(gca,cmap_pearson);
    caxis([0 0.5]);
    axis off
    set(gca,'xcolor','w','ycolor','w','zcolor','w');
    view(180,90);
    text(max(Vertices(:,1)),min(Vertices(:,2)),max(Vertices(:,3)),'D')
    %2-----------------------------------------------------------------------
    if band == 1
        title('PCorr(p-val<=0.05*)')
    end
    
    %3-----------------------------------------------------------------------
    subplot(10,4,8*(band-1)+3)
    patch('Faces',Faces,'Vertices',Vertices,'FaceVertexCData',...
        map2,'FaceColor','interp','EdgeColor','none','FaceAlpha',.99);
    colormap(gca,cmap_pval);
    caxis([0 0.05]);
    axis off
    set(gca,'xcolor','w','ycolor','w','zcolor','w');
    view(180,0);
    text(max(Vertices(:,1)),max(Vertices(:,2)),max(Vertices(:,3)),'L')
    %3-----------------------------------------------------------------------
    if band == 1
        title('p-val<=0.05*')
    end
    
    %4-----------------------------------------------------------------------
    subplot(10,4,8*(band-1)+4)
    patch('Faces',Faces,'Vertices',Vertices,'FaceVertexCData',...
        map2,'FaceColor','interp','EdgeColor','none','FaceAlpha',.99);
    colormap(gca,cmap_pval);
    caxis([0 0.05]);
    axis off
    set(gca,'xcolor','w','ycolor','w','zcolor','w');
    view(180,90);
    text(max(Vertices(:,1)),min(Vertices(:,2)),max(Vertices(:,3)),'D')
    %4-----------------------------------------------------------------------
    if band == 1
        title('p-val<=0.05*')
    end
    
    %5-----------------------------------------------------------------------
    subplot(10,4,8*(band-1)+5)
    patch('Faces',Faces,'Vertices',Vertices,'FaceVertexCData',...
        map1,'FaceColor','interp','EdgeColor','none','FaceAlpha',.99);
    colormap(gca,cmap_pearson);
    caxis([0 0.5]);
    axis off
    set(gca,'xcolor','w','ycolor','w','zcolor','w');
    view(-90,0);
    text(min(Vertices(:,1)),max(Vertices(:,2)),max(Vertices(:,3)),'P')
    %5------------------------------------------------------------------------
    
    %6------------------------------------------------------------------------
    subplot(10,4,8*(band-1)+6)
    patch('Faces',Faces,'Vertices',Vertices,'FaceVertexCData',...
        map1,'FaceColor','interp','EdgeColor','none','FaceAlpha',.99);
    colormap(gca,cmap_pearson);
    caxis([0 0.5]);
    axis off
    set(gca,'xcolor','w','ycolor','w','zcolor','w');
    view(180,-90);
    text(max(Vertices(:,1)),max(Vertices(:,2)),min(Vertices(:,3)),'V')
    %6------------------------------------------------------------------------
    
    %7------------------------------------------------------------------------
    subplot(10,4,8*(band-1)+7)
    patch('Faces',Faces,'Vertices',Vertices,'FaceVertexCData',...
        map2,'FaceColor','interp','EdgeColor','none','FaceAlpha',.99);
    colormap(gca,cmap_pval);
    caxis([0 0.05]);
    axis off
    set(gca,'xcolor','w','ycolor','w','zcolor','w');
    view(-90,0);
    text(min(Vertices(:,1)),max(Vertices(:,2)),max(Vertices(:,3)),'P')
    %7------------------------------------------------------------------------
    
    %8------------------------------------------------------------------------
    subplot(10,4,8*(band-1)+8)
    patch('Faces',Faces,'Vertices',Vertices,'FaceVertexCData',...
        map2,'FaceColor','interp','EdgeColor','none','FaceAlpha',.99);
    colormap(gca,cmap_pval);
    caxis([0 0.05]);
    axis off
    set(gca,'xcolor','w','ycolor','w','zcolor','w');
    view(180,-90);
    text(max(Vertices(:,1)),max(Vertices(:,2)),min(Vertices(:,3)),'V')
    %8------------------------------------------------------------------------
end
file_name = strcat(method_label,'_Temporal_PCorr','.fig');
savefig(mapsfig2,file_name);
%%
mapsfig3             = figure('Color','w','Name',strcat(method_label,'  Spatial P/Corr'));
for band = 1:length(band_label)
    J1                   = J{band}(:,1);
    J2                   = J{band}(:,2);
    map1                 = J1(ind);
    map2                 = J2(ind);
    map12                = map1.*map2;
    map1(map12 < 1E-8)   = [];
    map2(map12 < 1E-8)   = [];
    map1                 = log(map1);
    map2                 = log(map2);

%     J1                   = space_freq.J1(:,band);
%     map1                 = J1(ind);
%     map1                 = log(map1 - min(map1) + 1);
%     J2                   = space_freq.J2(:,band);
%     map2                 = J2(ind);
%     map2                 = log(map2 - min(map2) + 1);
    mdl                  = LinearModel.fit(map1,map2,'VarNames',{'ecog','eeg'});
    subplot(5,1,band); plot(mdl);
    cor                  = num2str(round(space_freq.cor(band),4));
    pvc                  = space_freq.pvc(band);
    if pvc > 0.05
        p_val = ">0.05ns";
    elseif (pvc > 0.01) && (pvc <= 0.05)
        p_val = "<=0.05*";
    elseif (pvc > 0.001) && (pvc <= 0.01)
        p_val = "<=0.01**";
    elseif (pvc > 0.0001) && (pvc <= 0.001)
        p_val = "<=0.001***";
    elseif pvc <= 0.0001
        p_val = "<=0.0001****";
    end
    ttl = strcat('Corr:', cor, '  P-val:', p_val);
    title(ttl);
end
file_name = strcat(method_label,'_Spatial_PCorr','.fig');
savefig(mapsfig3,fullfile(file_name));

end
