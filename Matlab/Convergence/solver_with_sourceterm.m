m = [];
error = [];

% % Spatial: by running to steady state
% % Lax-Wendoff (Without source and decay term; initialise with steaty)
%  m(           1 ) =          16 ; error(           1 ) =    1.6233375109468717E-002
%  m(           2 ) =          32 ; error(           2 ) =    3.5645582156477048E-003
%  m(           3 ) =          64 ; error(           3 ) =    1.0249151492071390E-003
%  m(           4 ) =         128 ; error(           4 ) =    3.1706854097251309E-004
%  m(           5 ) =         256 ; error(           5 ) =    1.2111513691564434E-004 
% 
% 
% 
%  % MC limiter
%  m(           1 ) =          16 ; error(           1 ) =    7.5472823218433668E-002
%  m(           2 ) =          32 ; error(           2 ) =    2.1884791443304644E-002
%  m(           3 ) =          64 ; error(           3 ) =    7.2307637677371207E-003
%  m(           4 ) =         128 ; error(           4 ) =    2.4975713827545504E-003
%  m(           5 ) =         256 ; error(           5 ) =    8.8271812906080902E-004
% 
% 
% % Upwind
%  m(           1 ) =          16 ; error(           1 ) =   0.11036916877298189     
%  m(           2 ) =          32 ; error(           2 ) =    9.2844197535228670E-002
%  m(           3 ) =          64 ; error(           3 ) =    6.4660942979935257E-002
%  m(           4 ) =         128 ; error(           4 ) =    3.9549298784285938E-002
%  m(           5 ) =         256 ; error(           5 ) =    2.2226556727427956E-002
% 
% % Spatial: with decay term; no diffusion
%  m(           1 ) =          16 ; error(           1 ) =    1.8299549097049810E-003
%  m(           2 ) =          32 ; error(           2 ) =    4.9596145787407487E-004
%  m(           3 ) =          64 ; error(           3 ) =    1.2692588355876766E-004
%  m(           4 ) =         128 ; error(           4 ) =    3.2231685049711650E-005
%  m(           5 ) =         256 ; error(           5 ) =    8.3737804688066507E-006

% Spatial: advection + diffusion + source
% nts = 512*m
% MC Limiter
 m(           1 ) =          16 ; error(           1 ) =    4.4328073570370734E-006
 m(           2 ) =          32 ; error(           2 ) =    6.5107311067579653E-007
 m(           3 ) =          64 ; error(           3 ) =    1.1730684828896538E-007
 m(           4 ) =         128 ; error(           4 ) =    2.9795349086797215E-008
 m(           5 ) =         256 ; error(           5 ) =    8.6683011857651947E-009

% Lax-Wendoff
 m(           1 ) =          16 ; error(           1 ) =    4.4693281090053769E-006
 m(           2 ) =          32 ; error(           2 ) =    1.1317049076188292E-006
 m(           3 ) =          64 ; error(           3 ) =    2.8291235286278855E-007
 m(           4 ) =         128 ; error(           4 ) =    7.0592857764786114E-008
 m(           5 ) =         256 ; error(           5 ) =    1.7885161078856209E-008
 
 m = [];
 error = [];
%   Timestepping to periodic state
 nts(           1 ) =       16384 ; error(           1 ) =    7.9150723590806441E-003
 nts(           2 ) =       32768 ; error(           2 ) =    2.2028023524653930E-003
 nts(           3 ) =       65536 ; error(           3 ) =    5.6600358679491250E-004
 nts(           4 ) =      131072 ; error(           4 ) =    1.4297735330849590E-004

 
h = 1./m;
h = 1./nts;

e_vis = error(2:end)./error(1:end-1);
h_vis = h(2:end)./h(1:end-1);

subplot(1,2,1)
plot(1:length(e_vis),log(e_vis)./log(h_vis))
subplot(1,2,2)
loglog(h,error)

p = polyfit(log(h),log(error),1)


%N = kg * m^2/2


% Dark Magic
88/5.5
64/4