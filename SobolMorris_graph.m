%% ------------------------------------------------------------------- %%
%              Script that implements the graph of 
%                    SOBOL and MORRIS INDICES 
% --------------------------------------------------------------------

% ------------------------------------------------------------------- %%
%   Sobol - Morris for MEDIAN QoI
% ----------------------------------------

figure('Name','Paragone Sobol vs Morris - Mediana')
for j = 1:N_fwd

    subplot(N_fwd,1,j)                                          
    
    % Total Sobol - light blue
    plot(tempo,Sob_idx_M(j,:),'-o','Color',"#4DBEEE",'LineWidth',1.5,'DisplayName','Sobol total indices')                     
    hold on
    grid on

    % Sobol - blue
    plot(tempo,Tot_Sob_idx_M(j,:),'-o','Color',"#0072BD",'LineWidth',1.5,'DisplayName','Sobol indices')                  
    hold on
   
    % Morris - green
    plot(tempo,emme_Morris_M(j,:)./emme_tot_M,'-o','Color',"#77AC30",'LineWidth',1.5,'DisplayName','Morris indices with modulus')     
    

    ylim([0 1])
    ylabel('Indices')
    title([param_6{j}])

end 
xlabel('Time [s]')
legend('Position', [0.33, 0.9, 0.1, 0.1]); 
legend show
sgtitle('Comparison Sobol vs Morris - MEDIAN')


%% ------------------------------------------------------------------- %%
%   Sobol - Morris for INTEGRAL QoI
% ----------------------------------------

figure('Name','Paragone Sobol vs Morris - Integrale')
for j = 1:N_fwd

    subplot(N_fwd,1,j)                                          
    
    % Totali Sobol - light blue
    plot(tempo,Sob_idx_I(j,:),'-o','Color',"#4DBEEE",'LineWidth',1.5,'DisplayName','Sobol total indices')                     
    hold on
    grid on

    % Sobol - blue
    plot(tempo,Tot_Sob_idx_I(j,:),'-o','Color',"#0072BD",'LineWidth',1.5,'DisplayName','Sobol indices')                  
    hold on
    
    % Morris - green
    plot(tempo,emme_Morris_I(j,:)./emme_tot_I,'-o','Color',"#77AC30",'LineWidth',1.5,'DisplayName','Morris indices with modulus')     
    

    ylim([0 1])
    ylabel('Indices')
    title([param_6{j}])

end 
xlabel('Time [s]')
legend('Position', [0.33, 0.9, 0.1, 0.1]); 
legend show
sgtitle('Comparison Sobol vs Morris - INTEGRAL')