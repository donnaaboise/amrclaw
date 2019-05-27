function [xp,yp,zp] = mapc2m(xc,yc)

global square

map = 'twisted_annulus';
square = false;

beta = 0.2;
twist = -0.02;
twist = load('twist.dat');

switch map
    case 'annulus'
       L = eye(2);
    case 'twisted_annulus'
       L = [1 twist; 0 1];
end

a = @(x,y) L(1,1)*x + L(1,2)*y;
b = @(x,y) L(2,1)*x + L(2,2)*y;


[xp,yp,zp] = mapc2m_annulus(a(xc,yc),b(xc,yc),beta);
% square = true;
% xp = xc;
% yp = yc;
% zp = 0*xc;

% xp = a(xc,yc);
% yp = b(xc,yc);
% zp = 0*xp;

end
