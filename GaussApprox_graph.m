%% ------------------------------------------------------------------- %%
%              Script that implements the graphs for 
%                      GAUSSIAN APPROXIMATION
% -------------------------------------------------------------------- %%

% sistemo incompatibilità di nomenclatrura tra i due script
if ~exist('parametri', 'var')
    parametri = param_5;
end 


%% -----------------------------------------%
%     Plot  MARGINALI APPROX GAUSSIANA
% -------------------------------------------

figure('Name','Marginals - Approx Gaussiana')
set(figure, 'Color', 'white');
for i=1:N_inv

    subplot (2,5,i)
    plot(ki(i,:),enne(i,:),'LineWidth',2,'Color',"#0072BD",'DisplayName','Gauss approx')     % blu
    hold on
    xline(y_star(i),'LineWidth',1.5,'Color',"#D95319",'DisplayName','y-star')                  % rosso
    xline(y_MAP_min(i),'c','LineWidth',2,'DisplayName','y-MAP')                % azzurro
    title([parametri(i)], 'FontSize', 14, 'FontWeight', 'bold')

end 
legend show
legend('Position', [0.7, 0.1, 0.1, 0.05]);
legend('FontSize', 12)
sgtitle('Marginals - Gaussian Approximation','FontWeight', 'bold')

