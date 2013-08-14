axis([0 1 0 1]);
daspect([1 1 1]);

% yrbcolormap;

showpatchborders;

hidegridlines();
showpatchborders();
hidecontourlines;
setpatchborderprops(1:5,'linewidth',2);
showgridlines(5);

caxis([0 1]);

fprintf('%10s : %24.16f\n','qmin',qmin);
fprintf('%10s : %24.16f\n','qmax',qmax);

yrbcolormap;
% Add extra colors to highlight undershoots and overshoots.
colormap([[0 1 1]; colormap; [1 0 1]]);

colorbar;
o = findobj('Tag','Colorbar');

set(o,'ytick',[qmin qmax])
set(o,'yticklabel',[qmin qmax]);
set(o,'fontsize',16,'fontweight','bold')


clear afterframe;
