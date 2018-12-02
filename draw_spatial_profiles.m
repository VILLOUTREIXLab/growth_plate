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

end

