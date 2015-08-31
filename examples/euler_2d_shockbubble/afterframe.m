if (PlotType ~= 3)
    colormap(jet)
end

axis([0 2.0 0 0.5])
axis equal; 
axis tight

fprintf('%10s : %12.4f\n','qmin',qmin);
fprintf('%10s : %12.4f\n','qmax',qmax);

if (PlotType == 3)
    caxis([0 200]);
else
    caxis([0.1 2.81]);
end

showpatchborders
setpatchborderprops('linewidth',1);
hidegridlines(1:3);

tstr = sprintf('Clawpack 5.0 : t = %8.4f\n',t);
title(tstr,'fontsize',16);
set(gca,'fontsize',16);


prt = true;
NoQuery = 0;
if (prt)
    MaxFrames = 31;
    axis off;

    id = input('Input id to use : ');
    if (~isempty(id))
        axis off;
        delete(get(gca,'title'));
        
        set(gcf,'papersize',[4,1]);
        set(gca,'position',[0 0 1 1]);
        set(gcf,'paperposition',[0 0 4 1]);
    
        % With mesh
        setpatchborderprops('linewidth',1);
        hidegridlines;
        if (PlotType == 3)
            fname = sprintf('results_%03d/amr_sb_schlrn_%02d.png',id,Frame);            
        else
            fname = sprintf('results_%03d/amr_sb_%02d.png',id,Frame);
        end
        yn = 'y';
        if (exist(fname))
            str = sprintf('Are you sure you want to overwrite %s (y/[n]) ? ',...
                fname);            
            yn = input(str,'s');               
            if isempty(yn)
                yn = 'n';
            end
        end      
        if (strcmp(lower(yn),'y') == 1)
            fprintf('Printing %s\n',fname);
            print('-r512','-dpng',fname);
        end
        
        % no mesh
        hidegridlines;
        hidepatchborders;
        if (PlotType == 3)
            fname = sprintf('results_%03d/amr_sb_schlrn_nomesh_%02d.png',id,Frame);            
        else
            fname = sprintf('results_%03d/amr_sb_nomesh_%02d.png',id,Frame);
        end
        if (strcmp(lower(yn),'y') == 1)
            fprintf('Printing %s\n',fname);
            print('-r512','-dpng',fname);
        end
    end
    
end

clear afterframe
clear mapc2m
