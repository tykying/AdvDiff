%% Prescription of (Analytic) Diffusivity field

% Remark: Assumed domain = [0, 1]

% Input required:
% kappa_testcase: specify which kappa field to use
% kappa_scale: scale of kappa field

% Output:
% kappa(x, y): scalar diffusivity field
% Gradkappa(x, y): gradient of diffusivity field
% sqrt2kappa: sqrt(2*kappa(x))

% Accept x,y to be row vectors listing the positions of all particles
%% Computation
function kappa_fldStruct = Prescribed_DiffusivityTensorField(Profile, param)
kappa_scale = param.kappa_scale;


if contains(Profile, 'const')   % Constant diffusion
    sigma1 = 1;
    sigma2 = 1;
    phi = 0;
    if contains(Profile, 'skew')
        sigma1 = sqrt(3);
        sigma2 = sqrt(2);
        phi = 2*pi/5;
    end
    
    [Kxx, Kyy, Kxy]  = KCanon_to_KCarte(sigma1.^2, sigma2.^2, phi);
    
    BaseTensor = [Kxx, Kxy; Kxy, Kyy];

    kappa  = @(x, y) kappa_scale.* BaseTensor;
    Divkappa_r1 = @(x, y) 0.*x;
    Divkappa_r2 = @(x, y) 0.*y;
    Divkappa = @(x, y) [Divkappa_r1(x,y); Divkappa_r2(x,y)];
    
    param.kappa = kappa_scale.*BaseTensor;
    param.theta = [sqrt(kappa_scale)*sigma1, sqrt(kappa_scale)*sigma2, phi];
end

if strcmp(Profile, 'sinusoidal')
    Lx1 = pi/2; Ly1 = -pi;
    Lx2 = pi/3; Ly2 = pi/2;
    Lxp = pi/3; Lyp = pi/4;

    sigma1 = @(x, y) sin(Lx1*x + Ly1*y);
    sigma2 = @(x, y) cos(Lx2*x - Ly2*y);
    phi = @(x, y) (pi/2)*sin(Lxp*x).*sin(Lyp*y);
    
    % Formula from Kpolar_to_Kcart_vectorised
    Kxx = @(x, y) (cos(phi(x,y)) .*sigma1(x,y)).^2 + (sin(phi(x,y)) .* sigma2(x,y)).^2;
    Kyy = @(x, y) (sin(phi(x,y)) .*sigma1(x,y)).^2 + (cos(phi(x,y)) .* sigma2(x,y)).^2;
    Kxy = @(x, y)  cos(phi(x,y)) .* sin(phi(x,y)) .* (sigma1(x,y).^2 - sigma2(x,y).^2);

    BaseTensor = @(x, y) [Kxx(x,y), Kxy(x,y); Kxy(x,y), Kyy(x,y)];

    kappa  = @(x, y) kappa_scale.* BaseTensor(x, y);
    Divkappa_r1 = @(x, y) kappa_scale .* (-cos((1/2).*pi.*sin(Lxp.*x).*sin(Lyp.*y)) ...
        .*sin(Lx1.*x+Ly1.*y).^2.*pi.*Lxp.*cos(Lxp.*x).*sin(Lyp.*y).*sin((1/2).*pi ...
        .*sin(Lxp.*x).*sin(Lyp.*y)) +2.*cos((1/2).*pi.*sin(Lxp.*x).*sin(Lyp.*y)).^2 ...
        .*sin(Lx1.*x+Ly1.*y).*Lx1.*cos(Lx1.*x+Ly1.*y)+sin((1/2).*pi.*sin(Lxp.*x) ...
        .*sin(Lyp.*y)).*cos(Lx2.*x-Ly2.*y).^2.*pi.*Lxp.*cos(Lxp.*x).*sin(Lyp.*y) ...
        .*cos((1/2).*pi.*sin(Lxp.*x).*sin(Lyp.*y))-2.*sin((1/2).*pi.*sin(Lxp.*x) ...
        .*sin(Lyp.*y)).^2.*cos(Lx2.*x-Ly2.*y).*Lx2.*sin(Lx2.*x-Ly2.*y) ...
        -(1/2).*pi.*sin(Lxp.*x).*Lyp.*cos(Lyp.*y).*sin((1/2).*pi ...
        .*sin(Lxp.*x).*sin(Lyp.*y)).^2.*(sin(Lx1.*x+Ly1.*y).^2-cos(Lx2.*x-Ly2.*y).^2) ...
        +(1/2).*cos((1/2).*pi.*sin(Lxp.*x).*sin(Lyp.*y)).^2.*pi.*sin(Lxp.*x)...
        .*Lyp.*cos(Lyp.*y).*(sin(Lx1.*x+Ly1.*y).^2-cos(Lx2.*x-Ly2.*y).^2) ...
        +cos((1/2).*pi.*sin(Lxp.*x).*sin(Lyp.*y)).*sin((1/2).*pi.*sin(Lxp.*x) ...
        .*sin(Lyp.*y)).*(2.*sin(Lx1.*x+Ly1.*y).*Ly1.*cos(Lx1.*x+Ly1.*y) ...
        -2.*cos(Lx2.*x-Ly2.*y).*Ly2.*sin(Lx2.*x-Ly2.*y)) );
    Divkappa_r2 = @(x, y) kappa_scale .*(-(1/2).*pi.*Lxp.*cos(Lxp.*x).*sin(Lyp.*y) ...
        .*sin((1/2).*pi.*sin(Lxp.*x).*sin(Lyp.*y)).^2.*(sin(Lx1.*x+Ly1.*y).^2 ...
        -cos(Lx2.*x-Ly2.*y).^2)+(1/2).*cos((1/2).*pi.*sin(Lxp.*x).*sin(Lyp.*y)).^2 ...
        .*pi.*Lxp.*cos(Lxp.*x).*sin(Lyp.*y).*(sin(Lx1.*x+Ly1.*y).^2-cos(Lx2.*x-Ly2.*y).^2) ...
        +cos((1/2).*pi.*sin(Lxp.*x).*sin(Lyp.*y)).*sin((1/2).*pi.*sin(Lxp.*x).*sin(Lyp.*y)) ...
        .*(2.*sin(Lx1.*x+Ly1.*y).*Lx1.*cos(Lx1.*x+Ly1.*y)+2.*cos(Lx2.*x-Ly2.*y).*Lx2 ...
        .*sin(Lx2.*x-Ly2.*y))+sin((1/2).*pi.*sin(Lxp.*x).*sin(Lyp.*y)).*sin(Lx1.*x+Ly1.*y).^2 ...
        .*pi.*sin(Lxp.*x).*Lyp.*cos(Lyp.*y).*cos((1/2).*pi.*sin(Lxp.*x).*sin(Lyp.*y)) ...
        +2.*sin((1/2).*pi.*sin(Lxp.*x).*sin(Lyp.*y)).^2.*sin(Lx1.*x+Ly1.*y).*Ly1 ...
        .*cos(Lx1.*x+Ly1.*y)-cos((1/2).*pi.*sin(Lxp.*x).*sin(Lyp.*y)).*cos(Lx2.*x-Ly2.*y).^2. ...
        *pi.*sin(Lxp.*x).*Lyp.*cos(Lyp.*y).*sin((1/2).*pi.*sin(Lxp.*x).*sin(Lyp.*y)) ...
        +2.*cos((1/2).*pi.*sin(Lxp.*x).*sin(Lyp.*y)).^2.*cos(Lx2.*x-Ly2.*y).*Ly2.*sin(Lx2.*x-Ly2.*y));
    Divkappa = @(x, y) [Divkappa_r1(x,y); Divkappa_r2(x,y)];
    
    Divkappa_scale = kappa_scale;
    param.Lx1 = Lx1;
    param.Ly1 = Ly1;
    param.Lx2 = Lx2;
    param.Ly2 = Ly2;
    param.Lxp = Lxp;
    param.Lyp = Lyp;
end

if strcmp(Profile, 'TG_iso')
    Lx1 = pi; Ly1 = pi;
    Lx2 = pi; Ly2 = pi;
    
    % Formula from Kpolar_to_Kcart_vectorised
    Kxx = @(x, y) sin(Lx1*x).*sin(Ly1*y);
    Kyy = @(x, y) sin(Lx2*x).*sin(Ly2*y);
    Kxy = @(x, y) 0.0;
    
    BaseTensor = @(x, y) [Kxx(x,y), Kxy(x,y); Kxy(x,y), Kyy(x,y)];

    kappa  = @(x, y) kappa_scale.* BaseTensor(x, y);
    Divkappa_r1 = @(x, y) kappa_scale* Lx1*cos(Lx1*x).*sin(Ly1*y);
    Divkappa_r2 = @(x, y) kappa_scale* Ly2*sin(Lx2*x).*cos(Ly2*y);
    Divkappa = @(x, y) [Divkappa_r1(x,y); Divkappa_r2(x,y)];
    
    Divkappa_scale = kappa_scale;
    param.Lx1 = Lx1;
    param.Ly1 = Ly1;
    param.Lx2 = Lx2;
    param.Ly2 = Ly2;
end


kappa_fldStruct = struct('kappa', kappa, 'Divkappa', Divkappa,  'Profile', Profile, 'param', param);
end
