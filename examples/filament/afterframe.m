s = 0;
axis([-s 2+s -s 2+s])
daspect([1 1 1]);

if (PlotParallelPartitions == 0)
    yrbcolormap;
end
showpatchborders(1:9);
caxis([0 1]);

if (t > 0)
    hold on;
    N = 1e4;
    [xout,yout] = filament_soln(N,t);
    plot(xout,yout,'k','linewidth',2);
    hold off;
end

% colormap(white);
str = sprintf('AMRClaw : t = %6.2f',t);
title(str,'fontsize',14);

shg;

NoQuery = 0;
prt = false;
if (prt)
    MaxFrames = 31;
    axis([0 2 0 2]);
    axis off;
    delete(get(gca,'title'));
    figsize = [4,4];  % Should match size set in options

    % Use this with 'export_fig'
    set(gca,'position',[0 0 1 1]);
    set(gcf,'units','inches');
    set(gcf,'position',[1 7 figsize]);
    
    % Start printing
    id = input('Input id to use : ');
    if (~isempty(id) | id == 999)

        fname_soln_tikz = sprintf('filament_tikz_%04d.tex',Frame);
        fname_soln_tikz_dir = sprintf('results_%03d/filament_tikz_%04d.tex',id,Frame);
        create_filament_soln_tikz(fname_soln_tikz_dir,xout,yout,figsize,1024,1024);

        
        % No mesh
        hidegridlines;
        hidepatchborders;
        if (PlotType == 3)
            fname_prefix = sprintf('amr_adv_schlrn',Frame);            
        else
            fname_prefix = sprintf('amr_adv',Frame);
        end
        yn = 'y';
        fname_png = sprintf('results_%03d/%s_%04d_%02d.png',id,fname_prefix,Frame,plot_level);
        if (exist(fname_png))
            str = sprintf('Overwrite file %s (y/[n]) ? ',fname_png);
            yn = input(str,'s');
            if (isempty(yn))
                yn = 'n';
            end
        end
                
        if (strcmp(lower(yn),'y') == 1)
            fprintf('Printing %s\n',fname_png);
            % print('-dpng','-r256',fname_png);
            export_fig('-dpng','-transparent','-r256',...
               '-a1','-p0','-nocrop',fname_png);
            amrclaw = 1;
            create_tikz_plot(id,Frame,fname_prefix,amrclaw,fname_soln_tikz);
        end
                
    end
    
end

clear afterframe
clear mapc2m
clear mapc2m_squareddisk
clear mapc2m_pillowdisk

