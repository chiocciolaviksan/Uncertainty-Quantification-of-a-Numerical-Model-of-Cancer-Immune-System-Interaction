%% ================================================================== %%
%          Bayesiana Inversione and MCMC Slice Sampling
%                 QoI Integral - 5 parameters
% =================================================================== %%

close all
clear all
clc

% -----------------------------
% path and changing folders 
% cd ..
addpath(genpath(pwd))

% activation of the Parallel Toolbox
% activate_parallel


%% ------------------------------------------------------------------- %%
%  Initial data and parameters
% ---------------------------------------
run InitialDataSetting.m

mu = 900;

% --------------------------------------------------------------------
% creation of the 'problem' structure with 5 parameters
problem_5 = @(y) pbParameters(a,b,chi,mu,y(1),A,L,T,...
    u0,phi0,phix0,fluxu0,fluxu1,fluxp0,fluxp1,force_u, @(x,t) force_phi(x,t,y([2,3,4,5])));


%% -------------------------------------------------------------------- % 
%  Sparse Grid for INVERSE analysis
% ------------------------------------
w = 1;

% % creation of the nodes
knots1_inv = @(n) knots_CC(n,nu(1),nu(2));  
knots2_inv = @(n) knots_CC(n,kappa1(1),kappa1(2));  
knots3_inv = @(n) knots_CC(n,rho(1),rho(2));  
knots4_inv = @(n) knots_CC(n,c(1),c(2));  
knots5_inv = @(n) knots_CC(n,sigma1(1),sigma1(2));


% Number of parameters selected for the inverse analysis
N_inv = 5;
param_5 = {'\nu','k_1','\rho','c','\sigma_{\phi}'};

% Domain of the different parameters
domain_inv = [nu(1) kappa1(1) rho(1) c(1) sigma1(1); 
              nu(2) kappa1(2) rho(2) c(2) sigma1(2)]; 

% Reasonable square for the response surfaces
domainPost = [-5 -5 -5 -5 -5; 
               5  5  5  5  5];  

% labels for the plots
legend_labels = cell(1, length(ind));
leg_v1 = {'\nu','\rho','c'};            % 1 3 5
leg_v2 = {'k_1','c','\sigma_{\phi}'};     % 2 4 


%% -------------------------------------------------------------------- % 
%  Gauss Approximation parameters
% ------------------------------------
% additive Gaussian noise
sigma_eps = 2;  

% true value of the parameters
y_star = [200; 40; 27.5; 3*1000/4; 10];

% Number of synthetic noisy data points
N_dati = 10;
% N_dati = 20;
% N_dati = 50;
i_t = round(linspace(3,length(time),N_dati));

% discretization of the time domain
t_k = time(i_t);
J = length(i_t);


% ------------------------------------------- %
%  MCMC Slice Sampling parameters
% ------------------------------------
th = 3;
bur = 100;

% number of samples
N_MCMC_samples = 10000;  

sampler = 'slice';

% width
peso_w = [30 20 10 20 5];


%% =================================================================== %%
%                       BAYESIAN Inversion
% ==================================================================== %%

% 'true' value
u_star = M_reshape(problem_5(y_star),i_t);


% noisy data 
error = sigma_eps*randn(J,1);
u_tilde = u_star + error;


[u_graph,~,~,~,~,~,~,~] = Solutore_HDG_KS1D(problem_5(y_star));
t_graph = linspace(I_time(1),I_time(2),length(u_graph));



%% ---------------------------------------------------------------------- %
%     plot campionamento i_t
% ------------------------------

figure('Name','Sampling')
set(figure, 'Color', 'white');
plot(u_graph,'-o','Color',"#77AC30",'LineWidth',1.5)
hold on
plot(i_t, u_star,'ro', 'MarkerSize',7,'LineWidth',1.5)
hold on
plot(i_t, u_tilde,'*','Color',"#0072BD",'MarkerSize',7,'LineWidth',1.5)

ylabel('QoI Median','FontSize', 14)
xlabel('Time','FontSize', 14)
title(['We consider ', num2str(N_dati), ' data points'],'FontSize', 14,'FontWeight', 'bold')
legend('q^{star}','q^{star}_{k} (sampled)','q^{tilde}_{k} (noise)')
legend('FontSize', 14)



%% --------------------------------------------------------------------- %
%    SPARSE GRID GENERATION and EVALUATION  
% -------------------------------------------- %

lev2knots = @lev2knots_doubling;         
rule = @(ii) sum(ii-1);

S_inv = create_sparse_grid(N_inv,w,{knots1_inv, knots2_inv, knots3_inv, knots4_inv, knots5_inv},lev2knots,rule,[]);
Sr_inv = reduce_sparse_grid(S_inv);

M_5 = @(y) M_reshape(problem_5(y), i_t);

M_on_Sr_5 = evaluate_on_sparse_grid(M_5,S_inv,Sr_inv,[],[],[],2);  
M_approx_5 = @(y) interpolate_on_sparse_grid(S_inv,Sr_inv,M_on_Sr_5,y);



%% -------------------------------- %
%     GAUSSIAN APPROXIMATION
% --------------------------------- %

% definition of misfits
misfits = @(y) u_tilde - M_approx_5(y);  

% definition of the Negative Log-Likelyhood
NLL = @(y,sigma_noise) sum( misfits(y).^2  / (2*sigma_noise^2)) + J*log(sigma_noise) + J*0.5*log(2*pi);

% The mean can be taken as the mode of the posterior pdf or with the  
% minimization of the sum of the squared errors:
LS = @(y) sum(misfits(y).^2);


%% -------------------------------- %
%     Mean - y-MAP
% --------------------------------- %
% Minimizations of the negative-log-likelyhood
n_random = 100;

pti_start = zeros(n_random,N_inv);
for i=1:N_inv
    pti_start(:,i) = domain_inv(1,i)+(domain_inv(2,i)-domain_inv(1,i))*rand(n_random,1);
end 

y_MAP_check = zeros(N_inv,n_random);
LS_finalvalues = zeros(1,n_random);
out_flags = zeros(1,n_random);

tic
for i =1:n_random 
    [y_MAP_check(:,i),LS_finalvalues(i),out_flags(i)] = fminsearch(LS,pti_start(i,:)');
end 
toc


%% -------------------------------------------------------------
%   Removal of those that do not converge

flag_toremove = find(out_flags ~= 1);
rimossi = flag_toremove;
out_flags_backup = out_flags;
flag_toremove_backup = flag_toremove;
y_MAP_check_backup = y_MAP_check;

number_removed = length(flag_toremove_backup);
n_random = n_random - number_removed;

pti_start(flag_toremove,:) = [];
y_MAP_check(:,flag_toremove) = [];
LS_finalvalues(flag_toremove) = [];
out_flags(flag_toremove) = [];
flag_toremove = [];


%% -------------------------------------------------------------
%  Diagnostic plots with Trajectories
% --------------------------------------
run Minimizations_graph.m


%% --------------------------------------------------------------------- %
%   y_map and sigma_eps approximation
% ---------------------------------------

[~,ind_MAP_choose] = min(LS_finalvalues);
y_MAP_min = y_MAP_check(:,ind_MAP_choose);

y_MAP = y_MAP_min;
sigma_eps_approx = sqrt(mean(misfits(y_MAP).^2));
sigma_eps_approx_ideal = sqrt(mean(misfits(y_star).^2));


%% -------------------------------------------------------------
%  plot of the Negative-Log-likelyhood
% -----------------------------------------
run NLL_graph.m


%% -------------------------------------------------------------------- %
%  covariance matrix
% ---------------------
% The covariance matrix Sigma_post can be taken as the inverse of the Hessian 
% of the NLL at yMAP. A convenient approximation for Sigma_post is
%           Sigma_post ~ sigma_epsilon^2(Jy(u)'Jy(u))^-1

Jac_at_MAP = zeros(J,N_inv);

for i = 1:J
 Jac_at_MAP(i,:) = derive_sparse_grid(S_inv,Sr_inv,M_on_Sr_5(i,:),domain_inv,y_MAP)';
end

Sigma_post = sigma_eps_approx^2*inv(Jac_at_MAP'*Jac_at_MAP);


%% ----------------------- %
%   correlation matrix
% ----------------------- %

esse = sqrt(diag(Sigma_post));    

% Formation of the inverse standard deviation matrix
Ddi = inv(diag(esse));

% Calculation of the correlation matrix
R = Ddi * Sigma_post * Ddi;

for i=1:N_inv
    left(i) = (y_MAP(i)-(10*esse(i)));
    right(i) = (y_MAP(i)+(10*esse(i)));
    ki(i,:) = linspace(left(i),right(i),100);
end 

for i=1:N_inv
    enne(i,:) = normpdf(ki(i,:), y_MAP(i), esse(i));
end


%% --------------------------------------------
%           Plot Approx Gaussiana
% ---------------------------------------------
run GaussApprox_graph.m


%% --------------------------------------------------------------------- %
%  forward UQ analysis based on the posterior
% ------------------------------------------------

knotsPost_5 = @(n) knots_normal(n,0,1); 

rule = @(ii) sum(ii-1);
lev2knotsPost_5 = @lev2knots_lin;

SPost_5 = create_sparse_grid(N_inv,w,knotsPost_5,lev2knotsPost_5,rule);
SrPost_5 = reduce_sparse_grid(SPost_5);

H = chol(Sigma_post);


%% ---------------------------- %
%    Posterior I_post
% ------------------------
M = 50000;

IPost_5 = @(z) QoI_I_wrapper(problem_5(H'*z+y_MAP));
I_on_SrPost_5 = evaluate_on_sparse_grid(IPost_5,SPost_5,SrPost_5,[],[],[],10);

z_samples = randn(N_inv,M);
y_post_samples = H'*z_samples+y_MAP;

I_vals_randnPost_5 = interpolate_on_sparse_grid(SPost_5,SrPost_5,I_on_SrPost_5,z_samples);

for i=1:length(ind)
    [y_pdf_Post(i,:),x_pdf_Post(i,:)] = ksdensity(I_vals_randnPost_5(ind(i),:));
end 

%% ---------------------------- %
%    Prior I_prior
% --------------------

IPrior_5 = @(y) QoI_I_wrapper(problem_5(y));

I_on_SrPrior_5 = evaluate_on_sparse_grid(IPrior_5,S_inv,Sr_inv,[],[],[],10);

y_samples_inv = zeros(N_inv,M);
for i=1:N_inv
    y_samples_inv(i,:) = domain_inv(1,i)+(domain_inv(2,i)-domain_inv(1,i))*rand(1,M);
end 

I_vals_randPrior_5 = interpolate_on_sparse_grid(S_inv,Sr_inv,I_on_SrPrior_5,y_samples_inv);

% ks-density
y_pdf_Prior = zeros(length(ind),100);
x_pdf_Prior = zeros(length(ind),100);
for i=1:length(ind)
    [y_pdf_Prior(i,:),x_pdf_Prior(i,:)] = ksdensity(I_vals_randPrior_5(ind(i),:));
end 

Qoi_star = QoI_I_wrapper(problem_5(y_star)); 



%% ---------------------------------------------------------------------
%           Prior-based and Posterior-based pdf over time 
%                   Considering INTEGRAL QoI
% ----------------------------------------------------------------------

quali = [5 10 15];
ind_local = ind(quali);


figure('Name','I_Prior and I_Post')
for i=1:length(ind_local)

    subplot(1,length(ind_local),i)

    % prior - orange
    h = histogram(I_vals_randPrior_5(ind_local(i),:),'Normalization','pdf','FaceColor',"#D95319",'DisplayName','Histogram prior');  
    hold on
    plot(x_pdf_Prior(quali(i),:),y_pdf_Prior(quali(i),:),'Color',"#C14410",'LineWidth',2.5,'DisplayName','Pdf prior')
    
    % posterior - blue 
    l = histogram(I_vals_randnPost_5(ind_local(i),:),'Normalization','pdf','FaceColor',"#0072BD",'DisplayName','Histogram posterior');
    hold on
    plot(x_pdf_Post(quali(i),:),y_pdf_Post(quali(i),:),'Color',"#4DBEEE",'LineWidth',2.5,'DisplayName','Pdf posterior')
    
    
    % xline(Qoi_star(quali(i)),'color',"#EDB120",'LineWidth',2.5,'DisplayName','Qoi star')
   
    xlabel('Values')
    ylabel('Frequency')
    title(['T = ', num2str(time(ind_local(i))), ' [s]'])
    
    
end 

sgtitle('Prior-based and Posterior-based pdf over time')
legend show



%% ========================================================
%     Evaluation Metrics CV and CF for I_prior e I_post
% ---------------------------------------------------------

media_prior = mean(I_vals_randPrior_5(end,:));
stdev_prior = std(I_vals_randPrior_5(end,:));

media_post_ga = mean(I_vals_randnPost_5(end,:));
stdev_post_ga = std(I_vals_randnPost_5(end,:));

CV_prior = stdev_prior/abs(media_prior);
CV_post_ga = stdev_post_ga/abs(media_post_ga);

CF_gauss = stdev_prior/stdev_post_ga;


%% =================================================================== %%
%                      Inverse UQ using MCMC 
%                         Slice Sampling 
% ==================================================================== %%

s0 = y_MAP'; 
tic
[MCMC_samples,neval] = slicesample(s0,N_MCMC_samples,'logpdf',@(y) -NLL(y',sigma_eps_approx),'thin',th,'burnin',bur,'width',peso_w); 
toc

maxlag = 50;
autocorr = zeros(2*maxlag+1,N_inv);
for d = 1:N_inv
    autocorr(:,d) = xcorr(MCMC_samples(:,d)-mean(MCMC_samples(:,d)),maxlag,'normalized');
end


%% -----------------------------------
%     plot for MCMC Slice Sampling
% ------------------------------------
run MCMC_graph.m


%% -----------------------------------
%     plot of Global Results
% ------------------------------------
figure('Name', 'Results');

for i = 1:N_inv
    for j = 1:(i-1)

        %----- contour i,j
        subplot(N_inv, N_inv, (i-1)*N_inv + j)

        % MCMC samples
        plot(MCMC_samples(:,j), MCMC_samples(:,i), 'x', 'Color', "#77AC30", 'LineWidth', 1.5, 'DisplayName', 'MCMC samples');
        hold on

        % MCMC contour
        [Join, Xi] = ksdensity([MCMC_samples(:,j), MCMC_samples(:,i)]);
        contour(reshape(Xi(:,1), 30, 30), reshape(Xi(:,2), 30, 30), reshape(Join, 30, 30), 20, 'color', 'k', 'DisplayName', 'MCMC contour');
        
        % Gauss contour
        y_M = [y_MAP(j) y_MAP(i)];
        Sigma_sub = Sigma_post([j, i], [j, i]);
        contour(reshape(Xi(:,1),30,30),reshape(Xi(:,2),30,30),reshape(mvnpdf(Xi,y_M,Sigma_sub),30,30),5,'LineColor',"#0072BD",'LineWidth',1.5,'DisplayName','Gauss contour')

        % y-MAP-check
        plot(y_MAP_check(j,:), y_MAP_check(i,:), 'o', 'Color', '#C87828', 'LineWidth', 1.5, 'DisplayName', 'y-MAP-check');
        hold on
        % plot(y_MAP_check(j,rimossi), y_MAP_check(i,rimossi), 'mo', 'LineWidth', 1.5, 'DisplayName', 'y-MAP removed');
        hold on
        
        % minimum of the NLL
        plot(y_MAP_min(j), y_MAP_min(i), 'xc', 'LineWidth', 3, 'DisplayName', 'y-MAP minimum');
        

        % setting of the prior region limits
        xline(domain_inv(1,j),'k', 'HandleVisibility', 'off')
        hold on
        xline(domain_inv(2,j),'k', 'HandleVisibility', 'off')
        hold on
        yline(domain_inv(1,i),'k', 'HandleVisibility', 'off')
        hold on
        yline(domain_inv(2,i),'k', 'HandleVisibility', 'off')
        
        xlim([domain_inv(1,j) domain_inv(2,j)])
       
        xlabel(param_5(j))
        ylabel(param_5(i))
        title(sprintf('%s - %s', param_5{j}, param_5{i}))
        
        if ((i-1)*N_inv + j)==((N_inv*N_inv)-1)

            % first legend
            legend1 = legend('show');
            set(legend1, 'Position', [0.8, 0.6, 0.1, 0.1])
        end 

        
    end 

    %----- histogram i - diagonal
    subplot(N_inv, N_inv, (i-1)*(N_inv+1) + 1)
    
    % histogram
    histogram(MCMC_samples(:,i), 'Normalization', 'pdf', 'FaceColor', "#77AC30" ,'DisplayName', 'Histogram MCMC');
    hold on
    
    % ksdensity MCMC
    [y_ksMCMC, x_ksMCMC] = ksdensity(MCMC_samples(:,i));
    plot(x_ksMCMC, y_ksMCMC, 'color', '#4C6E1E','LineWidth', 1.5, 'DisplayName', 'ks-density MCMC');

    % gaussian approx marginals
    plot(ki(i,:), enne(i,:), 'color', "#0072BD", 'LineWidth', 1.5, 'DisplayName', 'Gaussian Approximation');
    hold on

    % y-star
    xline(y_star(i), 'color','#C83C4B', 'LineWidth', 2.5, 'DisplayName', 'y-star');

    % y-minimum
    xline(y_MAP_min(i), 'c', 'LineWidth', 2.5, 'DisplayName', 'y-minimum');

    % xline(domain_inv(1,i),'k')
    % xline(domain_inv(2,i),'k')
    
   
    xlim([domain_inv(1,i) domain_inv(2,i)])

    title(parametri(i))
end 

sgtitle('Comparison: Gauss Approx vs MCMC slice sampling')
legend2 = legend('show');
set(legend2, 'Position', [0.8, 0.8, 0, 0.1])



%% ----------------------------------------------------------------------
%    Posterior MCMC
% ------------------------

Q_surrMCMC = interpolate_on_sparse_grid(S_inv,Sr_inv,I_on_SrPrior_5,MCMC_samples');

pdf_post_MCMC_surr = zeros(length(ind),100);
xx_post_MCMC_surr = zeros(length(ind),100);
for i=1:length(ind)
    [pdf_post_MCMC_surr(i,:),xx_post_MCMC_surr(i,:)] = ksdensity(Q_surrMCMC(ind(i),:));
end 


%% ------------------------------------------------------------------
%   Prior-based and Posterior-based pdf over time MCMC
% -------------------------------------------------------------------

figure('Name','Plot of the pdf')
for i=1:length(ind_local)
    subplot(1,length(ind_local),i)
    
    % ---- prior orange -----
    % h = histogram(Q_surr(ind_local(i),:),'Normalization','pdf',FaceColor='m');
    % hold on
    % plot(xx_prior(quali(i),:),pdf_prior(quali(i),:),'-r','LineWidth',1.5,'DisplayName','prior')
    % hold on
    h = histogram(I_vals_randPrior_5(ind_local(i),:),'Normalization','pdf','FaceColor',"#D95319",'DisplayName','Histogram Prior');  
    hold on
    plot(x_pdf_Prior(quali(i),:),y_pdf_Prior(quali(i),:),'Color',"#C14410",'LineWidth',2.5,'DisplayName','Pdf Prior')
    
    % ---- gaussian approssimazion blue -----
    p = histogram(I_vals_randnPost_5(ind_local(i),:),'Normalization','pdf','FaceColor',"#0072BD",'DisplayName','Hist. Gauss Approx.');
    hold on
    plot(x_pdf_Post(quali(i),:),y_pdf_Post(quali(i),:),'Color',"#4DBEEE",'LineWidth',2.5,'DisplayName','Pdf Gauss Approx.')
    % azzurro "#4DBEEE"
    
    % ---- MCMC green -----
    l = histogram(Q_surrMCMC(ind_local(i),:),'Normalization','pdf','FaceColor', "#77AC30",'DisplayName','Hist. Posterior MCMC');
    hold on
    plot(xx_post_MCMC_surr(quali(i),:),pdf_post_MCMC_surr(quali(i),:),'Color', '#4C6E1E','LineWidth',2.5,'DisplayName','Pdf Posterior MCMC')
    
    % ---- posterior pde -----
    % hold on
    % plot(xx_post_MCMC,pdf_post_MCMC,'LineWidth',2,'DisplayName','posterior MCMC pde')
 
    xlabel('Values')
    ylabel('Frequency')
    title(['T = ', num2str(time(ind_local(i))), ' [s]'])
end 

sgtitle('Prior-based and Posterior-based pdf over time')
legend show


%% ========================================================
%     Evaluation Metrics CV and CF for I_prior e I_post
% ---------------------------------------------------------

media_prior = mean(I_vals_randPrior_5(end,:));
stdev_prior = std(I_vals_randPrior_5(end,:));

media_post_MC = mean(Q_surrMCMC(end,:));
stdev_post_MC = std(Q_surrMCMC(end,:));

CV_prior = stdev_prior/abs(media_prior);
CV_post_MC = stdev_post_MC/abs(media_post_MC);

CF_MC = stdev_prior/stdev_post_MC;


%% ======================================================================
%                   Evaluation Metrics CV and CF 
%   for the PARAMETERS' distributions rho_prior and rho_post
% ======================================================================

for i=1:N_inv

    mu_u(i) = (domain_inv(2,i)+domain_inv(1,i))/2;
    std_u(i) = sqrt(((domain_inv(2,i)-domain_inv(1,i))^2)/12);

    mu_mc(i) = mean(MCMC_samples(:,i));
    std_mc(i) = std(MCMC_samples(:,i));
end 

mu_ga = y_MAP';
std_ga = esse';


CV_u = std_u./abs(mu_u);
CV_ga = std_ga./abs(mu_ga);
CV_mc = std_mc./abs(mu_mc);

CF_ga = std_u./std_ga;
CF_mc = std_u./std_mc;


%% ------------------------------------------------------------------- %%
%              ================ FUNZIONI ================
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

% --------------------------------------------------------------- %
% u sampled with reshape
function u_k = M_reshape(problem,i_t)
    
    [u,~,~,~,~,~,~,~] = Solutore_HDG_KS1D(problem);

    u_k = reshape(u(i_t),length(i_t),1);

    
end 
