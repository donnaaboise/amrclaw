global square

setviews;

alpha = 0.4;
s = 1e-2;    
alim = [-1,1];
axis([alim alim]);
daspect([1 1 1]);
view(vtop)

showpatchborders(1:10);
setpatchborderprops('linewidth',1)
caxis([0,1])

fprintf('qmin = %24.16e\n',qmin);
fprintf('qmax = %24.16e\n',qmax);

if (mq == 3)
    % Plot the error
    ca = [-max([qmin,qmax]),max([qmin,qmax])];
    colormap(jet);
else    
    % Plot the solution
    yrbcolormap
    ca = [0, 1];
end

% colorbar;
caxis(ca)
hidegridlines
if square
    axis([0,1,0,1]);
end

showgridlines(2)

colormap(parula);
ca = [-max([qmin,qmax]),max([qmin,qmax])];
caxis(ca/1e5);
% caxis([-1,1]*1e-12);

%
NoQuery = 0;
prt = false;
if (prt)
  MaxFrames = 8;
  axis([0 1 0 1]);
  filename = sprintf('annulus_%04d.png',Frame)
  print('-dpng',filename);
end

shg

clear afterframe;
clear mapc2m;
