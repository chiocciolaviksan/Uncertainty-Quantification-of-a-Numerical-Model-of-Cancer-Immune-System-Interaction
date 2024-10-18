%% ------------------------------------------------------------------- %%
%               Script che implementa i grafici per 
%                 la Negative Log-Likelyhood NLL
% -------------------------------------------------------------------- %%

for i = 1:N_inv
    yi_temp(i,:) = linspace(domain_inv(1,i),domain_inv(2,i),100);
end

% Creazione delle griglie
[Y2, Y6] = meshgrid(yi_temp(2,:), yi_temp(5,:));
Y1 = y_star(1)*ones(size(Y2));
Y3 = y_star(3)*ones(size(Y2));
Y4 = y_star(4)*ones(size(Y2));

% each column is a evaluation point
yy_plot = [Y1(:)'; Y2(:)'; Y3(:)'; Y4(:)'; Y6(:)']; 

NLL_plot = zeros(size(yy_plot,2),1);
for i = 1:size(yy_plot,2)
    NLL_plot(i) = NLL(yy_plot(:,i),sigma_eps_approx);
end


%%
figure
set(figure, 'Color', 'white');
surf(Y2,Y6,reshape(NLL_plot,100,100))
c = colorbar;
xlabel(parametri(2), "FontSize", 14, "FontWeight", 'bold')
ylabel(parametri(5), "FontSize", 14, "FontWeight", 'bold')
title('Negative-Log-Likelyhood')
