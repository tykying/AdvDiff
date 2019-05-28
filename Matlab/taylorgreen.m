% Grid space
m = 128

x = linspace(0, m, m+1)+1
y = linspace(0, m, m+1)+1

[X, Y] = meshgrid(x,y)

u = sin(X).*cos(Y);
v = -cos(X).*sin(Y);

% subplot(1,2,1)
% contourf(X, Y, u)
% 
% subplot(1,2,2)
% contourf(X, Y, v)


figure(2)
subplot(1,2,1)
psi = (X/m+2*Y/m).*sin((X-1)*pi/(m)).*sin((Y-1)*pi/(0.5*m));
contourf(X, Y, psi)

subplot(1,2,2)
K = sin((X-1)*pi/m).*sin((Y-1)*pi/m);
contourf(X, Y, K)

%% Logical space [0, 1]
m = 128

x = linspace(0, 1, m+1)
y = linspace(0, 1, m+1)

[X, Y] = meshgrid(x,y)

figure(3)
subplot(1,2,1)
psi = (X+2*Y).*sin(X*pi).*sin(Y*2*pi);
contourf(X, Y, psi)

subplot(1,2,2)
K = sin(X*pi).*sin(Y*pi);
contourf(X, Y, K)
