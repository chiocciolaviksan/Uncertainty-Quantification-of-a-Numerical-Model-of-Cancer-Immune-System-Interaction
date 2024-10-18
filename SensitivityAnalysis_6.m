%% ================================================================== %%
%                     SENSITIVITY ANALYSIS 
%                   QoI Median - 6 parameters
% =================================================================== %%

close all
clear all
clc

% -----------------------------
% path and changing folders 
% cd ..
addpath(genpath(pwd))

% activation of the Parallel Toolbox
activate_parallel


%% ------------------------------------------------------------------- %%
%  Initial data and parameters
% ---------------------------------------
run InitialDataSetting.m


% --------------------------------------------------------------------
% creation of the 'problem' structure with 6 parameters
problem_6 = @(y) pbParameters(a,b,chi,y(1),y(2),A,L,T,...
    u0,phi0,phix0,fluxu0,fluxu1,fluxp0,fluxp1,force_u, @(x,t) force_phi(x,t,y([3,4,5,6])));


%% -------------------------------------------------------------------- % 
%  Sparse Grid for FORWARD analysis
% ------------------------------------

w = 1;
lev2knots = @lev2knots_doubling;         
rule = @(ii) sum(ii-1);    

% creation of the nodes
knots1 = @(n) knots_CC(n,mu(1),mu(2));  
knots2 = @(n) knots_CC(n,nu(1),nu(2));  
knots3 = @(n) knots_CC(n,kappa1(1),kappa1(2));  
knots4 = @(n) knots_CC(n,rho(1),rho(2));  
knots5 = @(n) knots_CC(n,c(1),c(2));  
knots6 = @(n) knots_CC(n,sigma1(1),sigma1(2)); 


% number of parameters selected for the forward analysis
N_fwd = 6;
param_6 = {'\mu','\nu','k_1','\rho','c','\sigma_{\phi}'}; 

% domain of the different parameters
domain_6 = [mu(1) nu(1) kappa1(1) rho(1) c(1) sigma1(1); 
            mu(2) nu(2) kappa1(2) rho(2) c(2) sigma1(2)]; 

% labels for the plots
legend_labels = cell(1,length(ind));
leg_v1 = {'\mu','k_1','c'};                      
leg_v2 = {'\nu','\rho', '\sigma_{\phi}'};      


% number of samples for the Morris analysis
K_Morris = 10000;   

% number of random samples for the Monte Carlo method
M = 50000;



%% ------------------------------------------- %
%    SPARSE GRID GENERATION and EVALUATION  
% -------------------------------------------- %

% Median as QoI
Med =  @(y) QoI_M_wrapper(problem_6(y));   

% Integral as QoI
I = @(y) QoI_I_wrapper(problem_6(y));    


S_6 = create_sparse_grid(N_fwd,w,{knots1, knots2, knots3, knots4, knots5, knots6},lev2knots,rule,[]);
Sr_6 = reduce_sparse_grid(S_6);

M_on_Sr_6 = evaluate_on_sparse_grid(Med,S_6,Sr_6,[],[],[],10);
I_on_Sr_6 = evaluate_on_sparse_grid(I,S_6,Sr_6,[],[],[],10);


%% -------------------------------------------------------------------- %%
%     RESPONSE SURFACES 
% ----------------------------- %
run ResponseSurf_graph.m


%% ------------------------------------------------------------------- %%
%    Monte Carlo Method
% ----------------------------- %

map = get_interval_map(domain_6(1,:),domain_6(2,:),'uniform');
y_samples = map( rand(N_fwd,M)*2-1 );

M_vals_rand_6 = interpolate_on_sparse_grid(S_6,Sr_6,M_on_Sr_6,y_samples);
I_vals_rand_6 = interpolate_on_sparse_grid(S_6,Sr_6,I_on_Sr_6,y_samples);


for i = 2:length(ind)
     [y_pdf_M(i,:),x_pdf_M(i,:)] = ksdensity(M_vals_rand_6(ind(i),:));
     [y_pdf_I(i,:),x_pdf_I(i,:)] = ksdensity(I_vals_rand_6(ind(i),:));
  
end 


%% -------------------------------------------------------------------- %% 
%    HISTOGRAMS and PDFs
% ----------------------------- %
run HisPdf_graph.m


%% ------------------------------------------ %%
%       SOBOL' and MORRIS' INDICES 
% -------------------------------------------

tempo = linspace(I_time(1),I_time(2),size(M_on_Sr_6,1));

% --------------------------
%  Sobol for MEDIAN QoI
for i = 1:size(M_on_Sr_6,1)
    [Sob_idx_M(:,i), Tot_Sob_idx_M(:,i)] = compute_sobol_indices_from_sparse_grid(S_6,Sr_6,M_on_Sr_6(i,:),domain_6,'legendre');

end


% --------------------------
%  Sobol for INTEGRAL QoI
for i = 1:size(I_on_Sr_6,1)
    [Sob_idx_I(:,i), Tot_Sob_idx_I(:,i)] = compute_sobol_indices_from_sparse_grid(S_6,Sr_6,I_on_Sr_6(i,:),domain_6,'legendre');
end


% --------------------------
%  Morris for MEDIAN QoI
[emme_Morris_M,mi_Morris_M,sigma_Morris_M] = Morris(S_6,Sr_6,M_on_Sr_6,domain_6,K_Morris);
emme_tot_M = sum(emme_Morris_M,1);
mi_tot_M = sum(mi_Morris_M,1);


% --------------------------
%  Morris for INTEGRAL QoI
[emme_Morris_I,mi_Morris_I,sigma_Morris_I] = Morris(S_6,Sr_6,I_on_Sr_6,domain_6,K_Morris);
emme_tot_I = sum(emme_Morris_I,1);
mi_tot_I = sum(mi_Morris_I,1);


%% ---------------------------------------- %
%   Graphs of Sobol and Morris indices
% ----------------------------------------- %
run SobolMorris_graph.m



%% ------------------------------------------------------------------- %%
%             ================ FUNCTIONS ================
% -------------------------------------------------------------------- %

% function for creating the 'problem' structure
function problem = pbParameters(a,b,chi,mu,nu,A,L,T,...
    u0,phi0, phix0,fluxu0,fluxu1,fluxp0,fluxp1,force_u,force_phi)
    
    problem.a = a;
    problem.b = b;
    problem.chi = chi;
    problem.mu = mu;
    problem.nu = nu;
    problem.A = A;
    problem.L = L;
    problem.T = T;
    
    problem.u0 = u0;            
    problem.phi0 = phi0;
    problem.phix0 = phix0;
    
    problem.fluxu0 = fluxu0;    
    problem.fluxu1 = fluxu1;
    problem.fluxp0 = fluxp0;
    problem.fluxp1 = fluxp1;
    
    problem.force_u = force_u;
    problem.force_phi = force_phi;

end

% --------------------------------------------------------------- %
% WRAPPER INTEGRAL
function I_only = QoI_I_wrapper(problem)

    [~,~,I_only,~,~,~,~,~,~] = Solutore_HDG_KS1D(problem);

end 

% --------------------------------------------------------------- %
% WRAPPER MEDIAN
function med_u = QoI_M_wrapper(problem)

    [med_u,~,~,~,~,~,~,~,~] = Solutore_HDG_KS1D(problem);

end 
