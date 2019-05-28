xl = linspace(0, 1, 16+1);
yl = xl;

[X, Y] = meshgrid(xl, yl);

K = (sin(2*pi*X).*sin(2*pi*Y)).^2 + sin(pi*X).*cos(pi*Y);

%contourf(X, Y, K.^2);
imagesc(K.^2)