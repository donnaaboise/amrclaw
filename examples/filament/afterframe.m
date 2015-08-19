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
    filename = framename(Frame,'filament0000','png');    
    print('-dpng','-r800',filename);
end;

clear afterframe
clear mapc2m
clear mapc2m_squareddisk
clear mapc2m_pillowdisk

