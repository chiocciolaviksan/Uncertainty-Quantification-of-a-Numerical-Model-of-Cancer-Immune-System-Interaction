%% ------------------------------------------------------------------- %%
%                  Script that implements the graphs for 
%              the control and the TRAJECTORIES of y_MAP_check
% --------------------------------------------------------------------

% sistemo incompatibilità di nomenclatrura tra i due script
if ~exist('parametri', 'var')
    parametri = param_5;
end 


%% -------------------------------------------------------------
%   Diagnostic plot with Trajectories
% --------------------------------------

figure('Name','Trajectories')
subplot 131
plot(pti_start(:,1),pti_start(:,2),'bx')
hold on
plot(y_star(1),y_star(2),'ro','MarkerSize', 5, 'LineWidth', 2)
plot(y_MAP_check(1,:),y_MAP_check(2,:),'k^')
for i=1:n_random
    plot([pti_start(i,1) y_MAP_check(1,i)],[pti_start(i,2) y_MAP_check(2,i)]);
end
xlabel(parametri(1))
ylabel(parametri(2))

subplot 132
plot(pti_start(:,3),pti_start(:,4),'bx')
hold on
plot(y_star(3),y_star(4),'ro', 'MarkerSize', 5, 'LineWidth', 2)
plot(y_MAP_check(3,:),y_MAP_check(4,:),'k^')
for i=1:n_random
    plot([pti_start(i,3) y_MAP_check(3,i)],[pti_start(i,4) y_MAP_check(4,i)]);
end
xlabel(parametri(3))
ylabel(parametri(4))

subplot 133
plot(pti_start(:,4),pti_start(:,5),'bx')
hold on
plot(y_star(4),y_star(5),'ro', 'MarkerSize', 5, 'LineWidth', 2)
plot(y_MAP_check(4,:),y_MAP_check(5,:),'k^')
for i=1:n_random
    plot([pti_start(i,4) y_MAP_check(4,i)],[pti_start(i,5) y_MAP_check(5,i)]);
end
xlabel(parametri(4))
ylabel(parametri(5))
legend('starting points','y-STAR','y-MAP')
sgtitle('Trajectories')


%% ---------------------------------------------------------------------
%    Diagnostic plot of minimizations
% ---------------------------------------
set(groot, 'defaultAxesColorOrder', 'remove');

figure('Name','Diagnostic of minimizations')
set(figure, 'Color', 'white');
for i=1:N_inv
    subplot(2,5,i)
    plot(y_MAP_check(i,:)/y_star(i),'o','LineWidth',1.5)
    hold on
    ylim([0 2])
    plot([1 n_random],[1 1],'-x','LineWidth',1.5)
    title([parametri{i}], 'FontSize', 14, 'FontWeight', 'bold')
end 
sgtitle('Minimizations of the NLL','FontWeight', 'bold')
legend('y-MAP-check','y-STAR')
legend('FontSize', 12);