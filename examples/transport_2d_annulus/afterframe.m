setviews;

alpha = 0.4;
s = 1e-2;    
alim = [-1-alpha,1+alpha];
axis([alim alim]);
daspect([1 1 1]);
view(vtop)

% Use this for the caxis
ca = [-max(abs([qmin,qmax])),max(abs([qmin,qmax]))];

showpatchborders(1:10);
setpatchborderprops('linewidth',1)
caxis([0,1])

fprintf('qmin = %24.16e\n',qmin);
fprintf('qmax = %24.16e\n',qmax);

hidepatchborders(9)
showpatchborders;

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

caxis([-1,1]*1e-12);
% axis([-1,0,-1,0])
showgridlines

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
