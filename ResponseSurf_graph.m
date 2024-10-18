%% ------------------------------------------------------------------- %%
%              Script that implements the graphs of the 
%                        RESPONSE SURFACES
% -------------------------------------------------------------------- %%

% ------------------------------------------ %
%    Response Surfaces for MEDIAN QoI
% ------------------------------------------ %

if ~exist('couples', 'var')
    couples = 1:N_fwd;
    num = floor(length(couples)/2);
end 

colors = colormap('default');
v = ceil(linspace(1,size(colors,1),length(ind)));
colors = colors(v,:);

offset = 0.2;               
lineWidth = 2;              

for j = 1:num
    figure(j)
    set(figure, 'Color', 'white');
    legend_labels = cell(length(ind), 1);  

    for i = 1:length(ind)
        legend_labels{i} = ['Time: ' num2str(time(ind(i)))];

        [V1,V2,out_x,out_y,out_FIP] = Myplot(S_6,Sr_6,domain_6,M_on_Sr_6(ind(i),:),[]);

        % Normalize the z values to create the gradient
        z_data = out_FIP{1, j};
        z_min = min(z_data(:));
        z_max = max(z_data(:));
        z_norm = (z_data - z_min) / (z_max - z_min);  % Normalizza i valori z tra 0 e 1

        % Offset to avoid overly dark colors
        z_norm = z_norm * (1 - offset) + offset;

        % Create an RGB matrix for each point on the surface
        c_data = zeros(size(z_data, 1), size(z_data, 2), 3);
        for k = 1:3
            c_data(:, :, k) = colors(i, k) * z_norm;
        end

        % Plot of the surface with color gradients based on the z values and a grid
        surf(out_x{1, j}, out_y{1, j}, z_data, c_data, 'FaceColor', 'interp', 'EdgeColor', 'k');  
     
        % colorbar
        hold on
    end

    % Creation of the legend with the lightest color (maximum value)
    legend_handles = gobjects(length(ind), 1);
    for i = 1:length(ind)
        legend_handles(i) = plot(NaN, NaN, 'Color', colors(i, :) * (1 - offset) + offset, ...
                                 'DisplayName', [legend_labels{i} ' [s]'], 'LineWidth', lineWidth);
    end
    legend(legend_handles);

    xlabel([leg_v1{j}], 'FontSize', 14, 'FontWeight', 'bold')
    ylabel([leg_v2{j}], 'FontSize', 14, 'FontWeight', 'bold')
    zlabel('QoI Median','FontSize', 14)
    sgtitle(['Response Surfaces - MEDIAN'],'FontWeight', 'bold')
end


%%
% ------------------------------------------ %
%    Response Surfaces for MEDIAN QoI
% ------------------------------------------ %


if ~exist('couples', 'var')
    couples = 1:N_fwd;
    num = floor(length(couples)/2);
end 

if ~exist('colors', 'var')
    colors = colormap('default');
    v = ceil(linspace(1,size(colors,1),length(ind)));
    colors = colors(v,:);

    offset = 0.2;               
    lineWidth = 2;              

end 

for j = 1:num
    figure(j)
    set(figure, 'Color', 'white');
    legend_labels = cell(length(ind), 1);  

    for i = 1:length(ind)
        legend_labels{i} = ['Time: ' num2str(time(ind(i)))];

        [V1,V2,out_x,out_y,out_FIP] = Myplot(S_6,Sr_6,domain_6,I_on_Sr_6(ind(i),:),[]);
        
        % Normalize the z values to create the gradient
        z_data = out_FIP{1, j};
        z_min = min(z_data(:));
        z_max = max(z_data(:));
        z_norm = (z_data - z_min) / (z_max - z_min);  % Normalizza i valori z tra 0 e 1

        % Offset to avoid overly dark colors
        z_norm = z_norm * (1 - offset) + offset;

        % Create an RGB matrix for each point on the surface
        c_data = zeros(size(z_data, 1), size(z_data, 2), 3);
        for k = 1:3
            c_data(:, :, k) = colors(i, k) * z_norm;
        end

        % Plot of the surface with color gradients based on the z values and a grid
        surf(out_x{1, j}, out_y{1, j}, z_data, c_data, 'FaceColor', 'interp', 'EdgeColor', 'k');  
     
        % colorbar
        hold on
    end

    % Creation of the legend with the lightest color (maximum value)
    legend_handles = gobjects(length(ind), 1);
    for i = 1:length(ind)
        legend_handles(i) = plot(NaN, NaN, 'Color', colors(i, :) * (1 - offset) + offset, ...
                                 'DisplayName', [legend_labels{i} ' [s]'], 'LineWidth', lineWidth);
    end
    legend(legend_handles);

    xlabel([leg_v1{j}], 'FontSize', 14, 'FontWeight', 'bold')
    ylabel([leg_v2{j}], 'FontSize', 14, 'FontWeight', 'bold')
    zlabel('QoI Median','FontSize', 14)
    sgtitle(['Response Surfaces - INTEGRAL'],'FontWeight', 'bold')
end