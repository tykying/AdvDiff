% Infer a global isotropic diffusivity using m grid point
% K_const30d: s =    3.000E+01 days ; dt =    1.200E+01 hours
% niter = 10

m = [8, 16, 20, 32];

kappa_exact = 0.00260416666667;
%_mn: without interpolation; %_m: with interpolation
kappa_mn = [0.005390367249492142, 0.003413746843747649, 0.003073960340380079, 0.00278136345991069];
kappa_m = [0.0036811954758914985, 0.0028584867806513573, 0.002732610583200991, 0.0026547685699128097];
error_mn = abs(kappa_mn-kappa_exact);
error_m = abs(kappa_m-kappa_exact);

plot(log2(m), log2(error_mn), '-o')
hold on
plot(log2(m), log2(error_m), '-o')
ylabel('log2(error)')
xlabel('log2(m)')

K_order_mn = (log(error_mn(4))-log(error_mn(2)))/(log(m(4))-log(m(2)))
K_order_m = (log(error_m(4))-log(error_m(2)))/(log(m(4))-log(m(2)))