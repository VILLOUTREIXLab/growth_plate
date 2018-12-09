function [] = draw_spatial_profiles(features,grid,opt)
%DRAW_SPATIAL_PROFILES Draws spatial profiles
%
names = fieldnames(features);

for f = 2:length(names),

    h = figure('units','normalized','outerposition',[0 0 0.5 1])
    valtemp = squeeze(features.(names{f}).vals(:,2,2));
    indtemp = find(valtemp);
    plot(features.(names{f}).vals(indtemp,2,2),grid.x_bins(indtemp),'LineWidth',5,'Color','k');
    xlim([min(features.(names{f}).vals(indtemp,2,2)) max(features.(names{f}).vals(indtemp,2,2))]);
    title(features.(names{f}).title,'Interpreter','None');
    if opt.save_figs,
        mkdir([opt.save_folder]);
        saveas(h,[opt.save_folder,names{f},'_spatial_profile']);
        saveas(h,[opt.save_folder,names{f},'_spatial_profile','.png']);
    end
    close all
    
end

% drawing log spatial profiles for volume and surface area or centroid 
% shift
if isfield(features,'centroid_shift')
    fi = 5;
else
    fi = (3:4); 
end

for f = fi
     h = figure('units','normalized','outerposition',[0 0 0.5 1])
    valtemp = squeeze(features.(names{f}).vals(:,2,2));
    indtemp = find(valtemp);
    semilogx(features.(names{f}).vals(indtemp,2,2),grid.x_bins(indtemp),'LineWidth',5,'Color','k');
%     xlim([min(log(features.(names{f}).vals(indtemp,2,2))) max(log(features.(names{f}).vals(indtemp,2,2)))]);
    title(['Log ',features.(names{f}).title],'Interpreter','None');
    if opt.save_figs,
        mkdir([opt.save_folder]);
        saveas(h,[opt.save_folder,names{f},'_spatial_profile_log']);
        saveas(h,[opt.save_folder,names{f},'_spatial_profile_log','.png']);
    end
    close all
end

% superimposing volume, surface area and sphericity
% non log and log
if ~isfield(features,'centroid_shift')
        h = figure('units','normalized','outerposition',[0 0 0.5 1])
    for f = 3:5
        valtemp = squeeze(features.(names{f}).vals(:,2,2));
    indtemp = find(valtemp);
    hold on,
    plot(features.(names{f}).vals(indtemp,2,2)/max(features.(names{f}).vals(indtemp,2,2)),grid.x_bins(indtemp),'LineWidth',5);
%     xlim([min(log(features.(names{f}).vals(indtemp,2,2))) max(log(features.(names{f}).vals(indtemp,2,2)))]);
    hold off,
    end
    legend({'Vol','Surf','Spher'});
    title('Volume, surface area, sphericity','Interpreter','None');
    if opt.save_figs,
        mkdir([opt.save_folder]);
        saveas(h,[opt.save_folder,'combined_vol_surf_spher_spatial_profile']);
        saveas(h,[opt.save_folder,'combined_vol_surf_spher_spatial_profile','.png']);
    end
    close all    
end

end

