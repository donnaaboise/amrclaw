% Compute area per grid cell.

% Cartesian grid :

apgc_cart = 128^2/4;

apgc_pillow = 128^2/pi;

apgc_squareddisk = 5*128^2/pi;

fprintf('\n');
fprintf('Wall clock time \n');
fprintf('Pillow disk     %12.2f\n',212);
fprintf('Squared disk    %12.2f\n',18);
fprintf('Cart disk       %12.2f\n',203);

fprintf('\n');
fprintf('Grid cells per unit area\n');
fprintf('Pillow disk     %12.2f\n',apgc_pillow);
fprintf('Squared disk    %12.2f\n',apgc_squareddisk);
fprintf('Cart disk       %12.2f\n',apgc_cart);

fprintf('\n');
fprintf('Area per grid cell\n');
fprintf('Pillow disk     %12.4e\n',1/apgc_pillow);
fprintf('Squared disk    %12.4e\n',1/apgc_squareddisk);
fprintf('Cart disk       %12.4e\n',1/apgc_cart);
