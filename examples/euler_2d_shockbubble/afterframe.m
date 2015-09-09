colormap(jet)

axis([0 2.0 0 0.5])
axis equal; 
axis tight

fprintf('%10s : %12.4f\n','qmin',qmin);
fprintf('%10s : %12.4f\n','qmax',qmax);

if (PlotType == 3)
    caxis([0 200]);
else
    if (mq == 5)
        caxis([0 1]);
    else
        caxis([0.1 2.81]);
    end
end

showpatchborders
setpatchborderprops('linewidth',1);
hidegridlines(1:3);
delete(get(gca,'title'));
delete(get(gca,'xlabel'));
delete(get(gca,'ylabel'));
tstr = sprintf('Clawpack 5.0 : t = %8.4f\n',t);
title(tstr,'fontsize',16);

prt = false;
NoQuery = 0;
if (prt)
    axis([0 2 0 0.5]);
    axis off;
    delete(get(gca,'title'));
    figsize = [8,2];  % Should match size set in valout_tikz.f
    
    % Need to set the figure size (not paper size) with
    % 'export_fig'
    set(gca,'position',[0 0 1 1]);
    set(gcf,'units','inches');
    set(gcf,'position',[1 8 figsize]);
    
    hidepatchborders;
    hidegridlines;

    id = input('Input id to use : ');
    if (~isempty(id) | (id == 999))    
        fname_prefix = 'amr_sb';
        if (PlotType == 3)
            fname_prefix = [fname_prefix,'_schlrn'];
        end
        yn = 'y';
        fname_png = sprintf('results_%03d/%s_%04d_%02d.png',...
            id,fname_prefix,Frame,plot_level);
        if (exist(fname_png,'file'))
            str = sprintf('Overwrite file %s (y/[n]) ? ',fname_png);
            yn = input(str,'s');
            if (isempty(yn))
                yn = 'n';
            end
        end
                
        if (strcmp(yn,'y') == 1)
            % We have to use 'export_fig' here to get the 
            % transparency right.  If I just set the figure
            % color to 'none', space between the patches
            % shows through (why?)
            fprintf('Printing %s\n',fname_png);            
            export_fig('-dpng','-transparent','-r1024',...
                '-a1','-nocrop',fname_png);
            plotgrid = 1;
            create_tikz_plot(id,Frame,fname_prefix,1);
        end                
    end        
end

clear afterframe
clear mapc2m
