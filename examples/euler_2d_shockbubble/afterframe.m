colormap(jet)

axis([0 2.0 0 0.5])
axis equal; 
axis tight

fprintf('%10s : %12.4f\n','qmin',qmin);
fprintf('%10s : %12.4f\n','qmax',qmax);

caxis([0.1 2.81]);

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
    MaxFrames = 11;
    fname = sprintf('sbamr_nomesh20_%2.2d.png',Frame);
    disp(fname);
    print('-dpng',fname);
end

clear afterframe
clear mapc2m
