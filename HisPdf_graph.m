%% ------------------------------------------------------------------- %% 
%               Script that implements the graph of 
%                       HISTOGRAMS and PDF
% --------------------------------------------------------------------


% ----------------------------------------- %
%   Histograms in SUBPLOT for MEDIAN QoI
% ----------------------------------------- %

set(groot, 'defaultAxesColorOrder', 'remove');
valori = cell(1, length(ind));   

figure('Name','Histogram and Pdf - subplot')
set(figure, 'Color', 'white');

for i = 2:length(ind)
    
    subplot(4,5,i)
    h = histogram(M_vals_rand_6(ind(i),:),'Normalization','pdf');
    valori{i} = h.Values;

    hold on
    plot(x_pdf_M(i,:),y_pdf_M(i,:),'r','LineWidth',1.5)
    
   
    xlabel('Values')
    ylabel('Frequency')
    title(['T ',num2str(time(ind(i))), ' [s]'])
   
end 

legend('Position', [0.13, 0.75, 0.1, 0.15]); 
legend('Histogram','Pdf')
sgtitle('Histogram and Pdf over time - MEDIAN')


%% -------------------------------------------------------------------- %
%   Histograms in SUBPLOT for INTEGRAL QoI
% ------------------------------------------ %

set(groot, 'defaultAxesColorOrder', 'remove');

figure('Name','Histogram and Pdf - subplot')
set(figure, 'Color', 'white');

for i = 2:length(ind)
    
    subplot(4,5,i)
    h = histogram(I_vals_rand_6(ind(i),:),'Normalization','pdf');
    valori{i} = h.Values;

    hold on
    plot(x_pdf_I(i,:),y_pdf_I(i,:),'r','LineWidth',1.5)
    
    xlabel('Values')
    ylabel('Frequency')
     title(['T ',num2str(time(ind(i))), ' [s]'])
    
end 

legend('Position', [0.13, 0.75, 0.1, 0.15]); 
legend('Histogram','Pdf')
sgtitle('Histogram and Pdf over time - INTEGRAL')


%% ------------------------------------------------------------------- %%
%    PDF hold on MEDIAN QoI
% -----------------------------

rainbowColors = summer(20);
set(groot, 'defaultAxesColorOrder', rainbowColors);
legend_labels = cell(1, (length(ind)-1));

figure('Name','Pdf - hold on - Median')
set(figure, 'Color', 'white');

for i = 2:length(ind)
    legend_labels{i-1} = [num2str(time(ind(i))) ' [s]'];
    plot(x_pdf_M(i,:),y_pdf_M(i,:),'LineWidth',1.5)
    hold on
   
end  
xlabel('Values','FontSize', 14)
ylabel('Pdf','FontSize', 14)
legend(legend_labels)
title('Pdf over time - MEDIAN','FontSize', 14,'FontWeight', 'bold')

% desktopPath = fullfile(getenv('HOME'), 'Desktop');
% saveas(gcf, fullfile(desktopPath, 'pdf_med.pdf'));
% set(groot, 'defaultAxesColorOrder', 'remove');


%% ------------------------------------------------------------------- %%
%    PDF hold on INTEGRAL QoI
% ------------------------------

rainbowColors = summer(20);
set(groot, 'defaultAxesColorOrder', rainbowColors);
legend_labels = cell(1, (length(ind)-1));

figure('Name','Pdf - hold on - Integral')
set(figure, 'Color', 'white');

for i = 2:length(ind)
    legend_labels{i-1} = [num2str(time(ind(i))) ' [s]'];
    plot(x_pdf_I(i,:),y_pdf_I(i,:),'LineWidth',1.5)
    hold on
   
end  
xlabel('Values','FontSize', 14)
ylabel('Pdf','FontSize', 14)
legend(legend_labels)
title('Pdf over time - INTEGRAL','FontSize', 14,'FontWeight', 'bold')

% desktopPath = fullfile(getenv('HOME'), 'Desktop');
% saveas(gcf, fullfile(desktopPath, 'pdf_int.pdf'));
% set(groot, 'defaultAxesColorOrder', 'remove');
