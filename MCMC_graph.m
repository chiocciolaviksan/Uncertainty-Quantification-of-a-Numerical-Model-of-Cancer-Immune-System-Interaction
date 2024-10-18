%% ------------------------------------------------------------------- %%
%                  Script that implements the graphs for 
%                              MCMC ANALYSIS
% --------------------------------------------------------------------

% sistemo incompatibilità di nomenclatrura tra i due script
if ~exist('parametri', 'var')
    parametri = param_5;
end

%% -----------------------------------
%   samples nella regione surrogata
% ------------------------------------

% verifica che i samples MCMC siano tutti all'interno della regione 
% di definizione del modello surrogato
figure('Name','Regione Surrogata')
set(figure, 'Color', 'white');
                                    
% nu - kappa1
subplot 131
plot(MCMC_samples(:,1),MCMC_samples(:,2),'x','Color',"#77AC30",'LineWidth',2)
xline(domain_inv(1,1),LineWidth=2)
xline(domain_inv(2,1),LineWidth=2)
yline(domain_inv(1,2),LineWidth=2)
yline(domain_inv(2,2),LineWidth=2)
hold on
plot(y_star(1),y_star(2),'o','color','#C83C4B','LineWidth',3)

xlabel(parametri(1),'FontSize', 14, 'FontWeight', 'bold')
ylabel(parametri(2),'FontSize', 14, 'FontWeight', 'bold')

% rho - c
subplot 132
plot(MCMC_samples(:,3),MCMC_samples(:,4),'x','Color',"#77AC30",'LineWidth',2)
xline(domain_inv(1,3),LineWidth=2)
xline(domain_inv(2,3),LineWidth=2)
yline(domain_inv(1,4),LineWidth=2)
yline(domain_inv(2,4),LineWidth=2)
hold on
plot(y_star(3),y_star(4),'o','color','#C83C4B','LineWidth',3)

xlabel(parametri(3),'FontSize', 14, 'FontWeight', 'bold')
ylabel(parametri(4),'FontSize', 14, 'FontWeight', 'bold')

% c - sigma
subplot 133
plot(MCMC_samples(:,4),MCMC_samples(:,5),'x','Color',"#77AC30",'LineWidth',2)
xline(domain_inv(1,4),LineWidth=2)
xline(domain_inv(2,4),LineWidth=2)
yline(domain_inv(1,5),LineWidth=2)
yline(domain_inv(2,5),LineWidth=2)
hold on
plot(y_star(4),y_star(5),'o','color','#C83C4B','LineWidth',3)

xlabel(parametri(4),'FontSize', 14, 'FontWeight', 'bold')
ylabel(parametri(5),'FontSize', 14, 'FontWeight', 'bold')

sgtitle(['Surrogate Region'],'FontWeight', 'bold')
legend('MCMC Samples','Surrogate region limits','','','','y_{star}')


%% ----------------------------
%     diagnostic plots 
% -----------------------------

figure('Name','Samples')
set(figure, 'Color', 'white');
for i=1:N_inv
    subplot(2,5,i)
    plot(MCMC_samples(:,i),'x','Color',"#77AC30",'LineWidth', 2)
    title(parametri(i),'FontSize', 14, 'FontWeight', 'bold')
end 
sgtitle('MCMC Samples','FontWeight', 'bold')

%%
% autocorrelation
figure('Name','Autocorrelation')
set(figure, 'Color', 'white');
for i=1:N_inv
    subplot(2,5,i)
    plot(-maxlag:maxlag,autocorr(:,i),'Color',"#77AC30",'LineWidth', 2)
    yline(0, 'k')
    ylim([-0.05 1])
    title(parametri(i),'FontSize', 14, 'FontWeight', 'bold')
end 
sgtitle('MCMC Autocorrelation','FontWeight', 'bold')


%% ---------------------------------------------------------------
%     samples over pdf
[X1,X2] = meshgrid(ki(1,:),ki(2,:));      % nu - kappa1

[X3,X4] = meshgrid(ki(3,:),ki(4,:));      % rho - centro

[X5,X6] = meshgrid(ki(4,:),ki(5,:));      % centro - sigma

% prendo gli X calcolati in precedenza (~ line 700)
X = [X1(:)'; X2(:)'; X3(:)'; X4(:)'; X6(:)'];

figure('Name','samples over pdf')
set(figure, 'Color', 'white');
subplot 131
% NLL_isolines = exp(-NLL(X,sigma_eps_approx));
% contour(X1,X2,reshape(NLL_isolines,100,100),50,'DisplayName','isolines exp(-NLL)');
hold on
grid on
plot(MCMC_samples(:,1),MCMC_samples(:,2),'x','color',"#77AC30",'LineWidth',1.5,'DisplayName','samples')
contour(X1,X2,reshape(mvnpdf([X1(:) X2(:)],y_MAP([1:2])',Sigma_post(1:2,1:2)),100,100),50,'LineColor',"#0072BD",'LineWidth',1.5,'DisplayName','gauss approx')
plot(y_MAP_check(1,:), y_MAP_check(2,:), 'o', 'Color', '#C87828', 'LineWidth', 1.5, 'DisplayName', 'y-MAP-check');
hold on
% plot(y_MAP_check(1,rimossi), y_MAP_check(2,rimossi), 'mo', 'LineWidth', 1.5, 'DisplayName', 'y-MAP removed');
% hold on
xlabel(parametri(1),'FontSize', 14, 'FontWeight', 'bold')
ylabel(parametri(2),'FontSize', 14, 'FontWeight', 'bold')
% legend show

subplot 132
% NLL_isolines = exp(-NLL(X,sigma_eps_approx));
% contour(X3,X4,reshape(NLL_isolines,100,100),50,'DisplayName','isolines exp(-NLL)');
hold on
grid on
plot(MCMC_samples(:,3),MCMC_samples(:,4),'x','color',"#77AC30",'LineWidth',1.5,'DisplayName','samples')
contour(X3,X4,reshape(mvnpdf([X3(:) X4(:)],y_MAP([3:4])',Sigma_post(3:4,3:4)),100,100),50,'LineColor',"#0072BD",'LineWidth',1.5,'DisplayName','gauss approx')
plot(y_MAP_check(3,:), y_MAP_check(4,:), 'o', 'Color', '#C87828', 'LineWidth', 1.5, 'DisplayName', 'y-MAP-check');
hold on
% plot(y_MAP_check(3,rimossi), y_MAP_check(4,rimossi), 'mo', 'LineWidth', 1.5, 'DisplayName', 'y-MAP removed');
xlabel(parametri(3),'FontSize', 14, 'FontWeight', 'bold')
ylabel(parametri(4),'FontSize', 14, 'FontWeight', 'bold')
% legend show

subplot 133
% NLL_isolines = exp(-NLL(X,sigma_eps_approx));
% contour(X5,X6,reshape(NLL_isolines,100,100),50,'DisplayName','isolines exp(-NLL)');
hold on
grid on
plot(MCMC_samples(:,4),MCMC_samples(:,5),'x','color',"#77AC30",'LineWidth',1.5,'DisplayName','samples')
contour(X5,X6,reshape(mvnpdf([X5(:) X6(:)],y_MAP([4:5])',Sigma_post(4:5,4:5)),100,100),50,'LineColor',"#0072BD",'LineWidth',1.5,'DisplayName','gauss approx')
plot(y_MAP_check(4,:), y_MAP_check(5,:), 'o', 'Color', '#C87828', 'LineWidth', 1.5, 'DisplayName', 'y-MAP-check');
hold on
% plot(y_MAP_check(4,rimossi), y_MAP_check(5,rimossi), 'mo', 'LineWidth', 1.5, 'DisplayName', 'y-MAP removed');
xlabel(parametri(4),'FontSize', 14, 'FontWeight', 'bold')
ylabel(parametri(5),'FontSize', 14, 'FontWeight', 'bold')
legend show

sgtitle('Analysis of non-convergence points','FontSize', 14, 'FontWeight', 'bold')


%% ---------------------------------------
%  histogram & ksdensities of marginals
% ----------------------------------------

figure('Name','Istogrammi')
for i=1:N_inv

    [mar_MCMC,xi_MCMC] = ksdensity(MCMC_samples(:,i));
    % m = mean(MCMC_samples(:,i));

    subplot(2,3,i)
    hh = histogram(MCMC_samples(:,i),'Normalization','pdf','FaceColor', "#77AC30" ,'DisplayName', 'histogram MCMC');
    hold on
    plot(xi_MCMC,mar_MCMC,'Color', '#4C6E1E','LineWidth', 2, 'DisplayName', 'ksdensity MCMC')
    xline(y_star(i),'Color','#C83C4B', 'LineWidth', 2.5, 'DisplayName', 'y-star')
    plot(ki(i,:),enne(i,:),'color', "#0072BD", 'LineWidth', 1.5, 'DisplayName', 'Gaussian approximation')
    title(parametri(i))
   
end 
legend show



%% -------------------------------------------------------------------
% joint ksdensity
[Join_a,Xi_a] = ksdensity([MCMC_samples(:,1) MCMC_samples(:,2)]);
[Join_b,Xi_b] = ksdensity([MCMC_samples(:,3) MCMC_samples(:,4)]);
[Join_c,Xi_c] = ksdensity([MCMC_samples(:,4) MCMC_samples(:,5)]);


figure('Name','Joint ksdensity')
set(figure, 'Color', 'white');

subplot(131)
plot(MCMC_samples(:,1), MCMC_samples(:,2), 'x', 'Color', "#77AC30", 'LineWidth', 2, 'DisplayName', 'Samples')
hold on
contour(reshape(Xi_a(:,1), 30, 30), reshape(Xi_a(:,2), 30, 30), reshape(Join_a, 30, 30), 20, 'color', 'k', 'DisplayName', 'Ksdensity MCMC');
contour(reshape(Xi_a(:,1),30,30),reshape(Xi_a(:,2),30,30),reshape(mvnpdf(Xi_a,y_MAP(1:2)',Sigma_post([1:2],[1:2])),30,30),5,'LineColor',"#0072BD",'LineWidth',1.5,'DisplayName','Gauss approx')
xlabel(parametri(1),'FontSize', 14, 'FontWeight', 'bold')
ylabel(parametri(2),'FontSize', 14, 'FontWeight', 'bold')
title(sprintf('%s - %s', parametri{1}, parametri{2}),'FontSize', 14, 'FontWeight', 'bold')

subplot(132)
plot(MCMC_samples(:,3), MCMC_samples(:,4), 'x', 'Color', "#77AC30", 'LineWidth', 2, 'DisplayName', 'Samples')
hold on
contour(reshape(Xi_b(:,1), 30, 30), reshape(Xi_b(:,2), 30, 30), reshape(Join_b, 30, 30), 20, 'color', 'k', 'DisplayName', 'Ksdensity MCMC');
contour(reshape(Xi_b(:,1),30,30),reshape(Xi_b(:,2),30,30),reshape(mvnpdf(Xi_b,y_MAP(3:4)',Sigma_post([3:4],[3:4])),30,30),5,'LineColor',"#0072BD",'LineWidth',1.5,'DisplayName','Gauss approx')
xlabel(parametri(3),'FontSize', 14, 'FontWeight', 'bold')
ylabel(parametri(4),'FontSize', 14, 'FontWeight', 'bold')
title(sprintf('%s - %s', parametri{3}, parametri{4}),'FontSize', 14, 'FontWeight', 'bold')

subplot(133)
plot(MCMC_samples(:,4), MCMC_samples(:,5), 'x', 'Color', "#77AC30", 'LineWidth', 2, 'DisplayName', 'Samples')
hold on
contour(reshape(Xi_c(:,1), 30, 30), reshape(Xi_c(:,2), 30, 30), reshape(Join_c, 30, 30), 20, 'color', 'k', 'DisplayName', 'Ksdensity MCMC');
contour(reshape(Xi_c(:,1),30,30),reshape(Xi_c(:,2),30,30),reshape(mvnpdf(Xi_c,y_MAP(4:5)',Sigma_post([4:5],[4:5])),30,30),5,'LineColor',"#0072BD",'LineWidth',1.5,'DisplayName','Gauss approx')
xlabel(parametri(4),'FontSize', 14, 'FontWeight', 'bold')
ylabel(parametri(5),'FontSize', 14, 'FontWeight', 'bold')
title(sprintf('%s - %s', parametri{4}, parametri{5}),'FontSize', 14, 'FontWeight', 'bold')
legend show 

sgtitle('Joint ksdensity MCMC with Gauss Approx. Overlapped','FontWeight', 'bold')

