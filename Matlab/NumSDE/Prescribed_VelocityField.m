%% Prescription of Velocity Field

% Input required:
% veloc_testcase: specify which velocity field to use
% Pe: Peclet Number, to determine the scale of velocity field

% Output:
% u(x, y): velocity field
function veloc_fldStruct = Prescribed_VelocityField(Profile, param)
Pe = param.Pe;
kappa_scale = param.kappa_scale;
U_scale = param.U_scale;
GradU_scale = param.GradU_scale;

if contains(Profile, 'childress_soward')
    % psi_scale = u_scale;
    % psi =  @(x, y) sin(Lxc*x)*sin(Lyc*y) + lambda*cos(Lxs*x)*cos(Lys*y)
    
    lambda = 3/4;
    Lxs = 4/3;
    Lys = 1;
    Lxc = 1/2;
    Lyc = 3/2;
    
    u1 = @(x, y) (U_scale(1))*(-Lys*sin(Lxs*x).*cos(Lys*y) + lambda* Lyc*cos(Lxc*x).*sin(Lyc*y));
    u2 = @(x, y) (U_scale(2))*( Lxs*cos(Lxs*x).*sin(Lys*y) - lambda* Lxc*sin(Lxc*x).*cos(Lyc*y));
    
    ux = @(x, y) Lxc*cos(Lxc*x)*Lyc.*cos(Lyc*y) - lambda* Lxs*sin(Lxs*x).*Lys*sin(Lys*y);
    uy = @(x, y)   sin(Lxc*x)*Lyc^2.*sin(Lyc*y) + lambda* cos(Lxs*x)*Lys^2.*cos(Lys*y);
    vx = @(x, y)  -Lxc^2*sin(Lxc*x).*sin(Lyc*y) - lambda* Lxs^2*cos(Lxs*x).*cos(Lys*y);
    vy = @(x, y) Lxc*cos(Lxc*x)*Lyc.*cos(Lyc*y) + lambda* Lxs*sin(Lxs*x).*Lys*sin(Lys*y);
    
    Gradu1 = @(x, y) (U_scale(1))* [ ux(x,y) , uy(x,y) ];
    Gradu2 = @(x, y) (U_scale(2))* [ vx(x,y) , vy(x,y) ];
    
    
    % Append param struct
    param.Lxs = Lxs;
    param.Lys = Lys;
    param.Lxc = Lxc;
    param.Lyc = Lyc;
    param.lambda = lambda;
    
    u = @(x, y) [u1(x,y); u2(x,y)];    % Assume x,y is row vector [Python-friendly]
    Gradu = @(x,y) [Gradu1(x,y); Gradu2(x,y)];
end

if contains(Profile, 'tg_w_mean')
    % Guarantee homogenisation to work
    psi_scale = param.vortix_strength;
    vortix_length = param.vortix_length;
    
    Lx = vortix_length;
    Ly = vortix_length;
    
    b1 = param.b_mag*cos(param.b_angle);
    b2 = param.b_mag*sin(param.b_angle);
    
    psi = @(x, y) psi_scale*sin(2*pi/Lx*x).*sin(2*pi/Ly*y) - b1*y + b2*x;
    
    u1 = @(x, y) (psi_scale*2*pi/Ly)*(-cos(2*pi/Ly*y).*sin(2*pi/Lx*x)) + b1;
    u2 = @(x, y) (psi_scale*2*pi/Lx)*( sin(2*pi/Ly*y).*cos(2*pi/Lx*x)) + b2;
    
    Gradu1 = @(x, y) (2*pi*psi_scale/Ly)* [ -2*pi/Lx*cos(2*pi/Ly*y).*cos(2*pi/Lx*x) , 2*pi/Ly*sin(2*pi/Ly*y).*sin(2*pi/Lx*x) ];
    Gradu2 = @(x, y) (2*pi*psi_scale/Lx)* [ -2*pi/Lx*sin(2*pi/Ly*y).*sin(2*pi/Lx*x) , 2*pi/Ly*cos(2*pi/Ly*y).*cos(2*pi/Lx*x) ];
    
    u1_mean = @(x, y) b1;
    u2_mean = @(x, y) b2;
    
    % Append param struct
    param.Lx = Lx;
    param.Ly = Ly;
    param.b = [b1; b2];
    
    u = @(x, y) [u1(x,y); u2(x,y)];    % Assume x,y is row vector [Python-friendly]
    Gradu = @(x,y) [Gradu1(x,y); Gradu2(x,y)];
    
    u_mean = @(x, y) [u1_mean(x,y); u2_mean(x,y)];    % Assume x,y is row vector [Python-friendly]
end

if contains(Profile, 'TTG')
    % Defined on [0, 1]
    U_scale = param.U_scale;
    
    a = 3;
    b = -1;
    Lx = pi;
    Ly = 2*pi;
    k = -0.5;
        
    psi = @(x, y) U_scale*exp(k*(a*x+b*y)/(a+b))*sin(Lx*x)*sin(Ly*y);
    
    u1 = @(x, y)  -U_scale*k*b*exp(k*(a*x+b*y)/(a+b)).*sin(Lx*x).*sin(Ly*y)/(a+b)-U_scale*exp(k*(a*x+b*y)/(a+b)).*sin(Lx*x)*Ly.*cos(Ly*y);
    u2 = @(x, y)   U_scale*k*a*exp(k*(a*x+b*y)/(a+b)).*sin(Lx*x).*sin(Ly*y)/(a+b)+U_scale*exp(k*(a*x+b*y)/(a+b)).*Lx.*cos(Lx*x).*sin(Ly*y);
    
    Gradu1 = @(x, y) [ 0.*x , 0.*y ];
    Gradu2 = @(x, y) [ 0.*x , 0.*y ];
    
    % Append param struct
    param.Lx = Lx;
    param.Ly = Ly;
    param.a = a;
    param.b = b;
    
    u = @(x, y) [u1(x,y); u2(x,y)];    % Assume x,y is row vector [Python-friendly]
    Gradu = @(x,y) [Gradu1(x,y); Gradu2(x,y)];
    
    u_mean = @(x, y) [u1(x,y); u2(x,y)];    % Assume x,y is row vector [Python-friendly]

end


if contains(Profile, 'cosine')
    vortix_length = param.vortix_length;
    
    L = param.L;
    Lx = 0;
    Ly = 5*vortix_length;  % vortix_length = 100km
    psi_scale = 3*U_scale *Ly/(2*pi);  % U_scale = 0.1 m/s
    %psi_scale = 20 * param.U_scale *Ly/(2*pi);
        
    psi = @(x, y) -psi_scale.*0.5*(sin(2*pi/Ly*y) + 2*pi/Ly*y);
    
    u1 = @(x, y) (psi_scale*2*pi/Ly).*0.5*(cos(2*pi/Ly*y) + 1);
    u2 = @(x, y) psi_scale*0.*y;
    
    Gradu1 = @(x, y) (psi_scale*0.5*(2*pi/Ly)^2)* [ 0 , -sin(2*pi/Ly*y) ];
    Gradu2 = @(x, y) (psi_scale)* [ 0 , 0 ];
    
    % Append param struct
    param.psi_scale = psi_scale;
    param.Lx = Lx;
    param.Ly = Ly;
    
    u = @(x, y) [u1(x,y); u2(x,y)];    % Assume x,y is row vector [Python-friendly]
    Gradu = @(x,y) [Gradu1(x,y); Gradu2(x,y)];
    
    u_mean = @(x, y) [u1(x,y); u2(x,y)];    % Assume x,y is row vector [Python-friendly]
end


if contains(Profile, 'linear_jet')
    nu_1 = 4;
    phi_A = -2*pi/5;

    if contains(Profile, '_CplxEigen')
        nu_2 = -sqrt(3);
    elseif contains(Profile, '_RealEigen')
        nu_2 = sqrt(3);
    end
    if contains(Profile, '_ISO')
        nu_2 = nu_2/10;
    end
    
    if contains(Profile, '_INC')
        TrA = 0;
    else
        TrA = 2;
    end
    
    s2A = sin(2*phi_A); 
    c2A = cos(2*phi_A);
    
    theta_nus =  [0.5*(nu_1+nu_2), 0.5*(nu_1-nu_2)];
    theta_trAh = 0.5*TrA;
    
    theta_A = [theta_nus, phi_A, theta_trAh];
    
    A_base = theta_A(1)*[-s2A, c2A; c2A, s2A] + theta_A(2)*[0 1; -1 0] + theta_A(4)*eye(2)
    
    A = GradU_scale * A_base;
 
    if nu_1*nu_2 > 0
        assert(isreal(eig(A)))
    else
        assert(~isreal(eig(A)))
    end
    
    b = param.b_mag .* [cos(param.b_angle); sin(param.b_angle)];
    
    % Append param struct

    theta = [param.b_mag, param.b_angle, GradU_scale*theta_nus, phi_A, GradU_scale*theta_trAh];
    param.A_base = A_base;
    param.A = A;
    param.b = b;
    param.nu_1 = nu_1;
    param.nu_2 = nu_2;
    param.phi_A = phi_A;
    param.TrA = TrA;
    param.theta = theta;
    
    u = @(x, y) A*[x; y] + b;    % Assume x,y is row vector [Python-friendly]
    Gradu = @(x,y) A;
    
    u_mean = @(x, y) A*[x; y] + b;    % Assume x,y is row vector [Python-friendly]
end


if contains(Profile, 'zero')
    psi = @(x, y) zeros(size(x));
    
    u1 = @(x, y) zeros(size(x));
    u2 = @(x, y) zeros(size(x));
    
    Gradu1 = @(x, y) [ 0 , 0 ];
    Gradu2 = @(x, y) [ 0 , 0 ];
        
    u = @(x, y) [u1(x,y); u2(x,y)];    % Assume x,y is row vector [Python-friendly]
    Gradu = @(x,y) [Gradu1(x,y); Gradu2(x,y)];
    
    u_mean = @(x, y) [u1(x,y); u2(x,y)];    % Assume x,y is row vector [Python-friendly]
end



% Wrap up for output
veloc_fldStruct = struct('u', u, 'Gradu', Gradu, 'u_mean', u_mean, 'psi', psi, 'Profile', Profile, 'param', param);


end
