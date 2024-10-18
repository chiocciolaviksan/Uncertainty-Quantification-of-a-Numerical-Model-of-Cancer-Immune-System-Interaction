%% ------------------------------------------------------------------- %%
%                   script that provides the 
%                        INITIAL SETTING
% --------------------------------------------------------------------


% -------------------------------- % 
%       TEMPORAL DOMAIN 
% -------------------------------- % 
T = 1800;       
I_time  = [0,T];                              

Mtime = 100;
dt = T/Mtime;
time = (dt:dt:T);

punti = 20;
ind = round(linspace(1,length(time),punti));

% -------------------------------- % 
%     SPACIAL DOMAIN = [A,A+L]
% -------------------------------- % 
A = 0;
L = 1000;

I_space = [A,A+L];                              
xx = linspace(I_space(1),I_space(2),punti);   

% -------------------------------- % 
%           PARAMETERS
% -------------------------------- % 
b = 0;  
nu_nominale = 200;
chi = 10/nu_nominale;
a = 1.e-04;

kappa_init = 50;            
sigma_init = 10.;           
c_init = L/4;


% -----------------
% range

mu = [700 1100];
nu = [100 300];   
kappa1 = [30 70];
rho = [5 50];
c = [650 850];
sigma1 = [5 15];


% -----------------
% nominal value

% mu = 900;             % visc. phi
% nu = 200;             % visc. u
% kappa1 = 1;         
% rho = 27.5            % (e-4)
% c = 3*1000/4; = 750         
% sigma1 = 10;      


% -------------------------------- % 
%    Initial condition for u
% -------------------------------- % 
u0 = @(x,t) kappa_init * exp( -((x-c_init).^2) / (2*sigma_init^2) );


% ------------------------------------------------ % 
%  Initial data for phi and its derivative in x
% ------------------------------------------------ % 
phi0 = @(x,t)  0 + 0*x;
phix0 = @(x,t)  0 + 0*x;


% ---------------------------------------------------- % 
%  Neumann boundary conditions for the two variables
% ---------------------------------------------------- % 
fluxu0 = @(t)  0 + 0*t;
fluxu1 = @(t)  0 + 0*t;
fluxp0 = @(t)  0 + 0*t;
fluxp1 = @(t)  0 + 0*t;


% -------------------------------- % 
%          FORCING TERMS
% -------------------------------- % 
force_u = @(x,t) zeros(size(x));
force_phi = @(x,t,z) z(1) * exp(-z(2)*1.e-4*t) * exp(- (x - z(3)).^2/(2*z(4)^2));


