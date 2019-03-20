function [xp,yp,zp] = mapc2m(xc,yc)

map = 'twisted_annulus';

switch map
    case 'annulus'
       L = eye(2);
    case 'twisted_annulus'
       L = [1 -0.2; 0 1];
end

a = @(x,y) L(1,1)*x + L(1,2)*y;
b = @(x,y) L(2,1)*x + L(2,2)*y;


beta = 0.2;

[xp,yp,zp] = mapc2m_annulus(a(xc,yc),b(xc,yc),beta);
% xp = a(xc,yc);
% yp = b(xc,yc);
% zp = 0*xp;

end
