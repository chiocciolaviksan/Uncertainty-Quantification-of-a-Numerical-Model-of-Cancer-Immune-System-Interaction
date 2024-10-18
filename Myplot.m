%----------------------------------------------------
%  Implementazione solo per il mio grafico
%----------------------------------------------------
function [V1,V2,out_x,out_y,out_FIP] = Myplot(S,Sr,domain,f_values,NP) 

    % PLOTS_SPARSE_GRIDS_INTERPOLANT plots the sparse grid interpolant of a function. Different plots
    % are produced depending on the number of dimensions of the sparse grid:
    
    % if N>3, a number of bidimensional cuts will be considered, and for each of them a surf will be generated.
    %       In other words, all variables but two will be frozen to their average value and the resulting
    %       two-dimensional plot will be produced
    
    N = size(domain,2);
    
    % ===================================== %
    %      dipende dal valore di N
    % ===================================== %
    if (N < 3)
        error('Usare la funzione apposita plot_sparse_grids_interpolant')
    elseif (N == 3)
        
        % ===================================== %
        %      caso N = 3
        % ===================================== %
        % value of 'plot_grid_size'
        NP = 20;

        % extract info on lower and upper ends of each direction
        aa_vec = domain(1,:);
        bb_vec = domain(2,:);
        avg_vec = (aa_vec + bb_vec) / 2;
        
        % wrap interpolate on sparse grid into a @-function for ease of plotting
        f_interp = @(x) interpolate_on_sparse_grid(S, Sr, f_values, x);

        % initialization
        CUTS = N;
        V1 = zeros(1, CUTS);
        V2 = zeros(1, CUTS);
        out_x = cell(1, CUTS);
        out_y = cell(1, CUTS);
        out_FIP = cell(1, CUTS);

        % per le direzioni uso una lookup table
        couple_loc = [1 2; 1 3; 2 3];

        for ii = 1:CUTS
            v1 = couple_loc(ii, 1);
            v2 = couple_loc(ii, 2);
            V1(ii) = v1;
            V2(ii) = v2;

            % generate a mesh grid over the cut
            xp = linspace(aa_vec(v1), bb_vec(v1), NP);
            yp = linspace(aa_vec(v2), bb_vec(v2), NP);
            
            [XP, YP] = meshgrid(xp, yp);
            out_x{ii} = XP;
            out_y{ii} = YP;
            
            nb_pts = length(xp) * length(yp);
           
            PTS = avg_vec' * ones(1, nb_pts);
            
            % replace lines of non-constant directions
            PTS(v1, :) = XP(:)';
            PTS(v2, :) = YP(:)';
            
            % interpolate on sparse grid
            f_interp_eval = f_interp(PTS);
            
            % reshape to use surf
            FIP = reshape(f_interp_eval, size(XP));
            out_FIP{ii} = FIP; 
        end 

    % ===================================== %
    %      caso N > 3 DISPARI
    % ===================================== %    
    elseif (N > 3 && mod(N, 2) ~= 0)
        
        % value of 'plot_grid_size'
        NP = 20;

        % extract info on lower and upper ends of each direction
        aa_vec = domain(1,:);
        bb_vec = domain(2,:);
        avg_vec = (aa_vec + bb_vec) / 2;
        
        % wrap interpolate on sparse grid into a @-function for ease of plotting
        f_interp = @(x) interpolate_on_sparse_grid(S, Sr, f_values, x);

        % initialization
        CUTS = floor(N / 2) + 1;
        V1 = zeros(1, CUTS);
        V2 = zeros(1, CUTS);
        out_x = cell(1, CUTS);
        out_y = cell(1, CUTS);
        out_FIP = cell(1, CUTS);
        
        % per le direzioni uso una lookup table
        couple_loc = zeros(CUTS, 2);

        % Prima coppie consecutive
        for i = 1:(CUTS - 1)
            couple_loc(i, :) = [2 * i - 1, 2 * i];
        end
        
        % Ultima coppia: penultimo e ultimo
        couple_loc(CUTS, :) = [N - 1, N];

        for ii = 1:CUTS
            v1 = couple_loc(ii, 1);
            v2 = couple_loc(ii, 2);
            V1(ii) = v1;
            V2(ii) = v2;

            % generate a mesh grid over the cut
            xp = linspace(aa_vec(v1), bb_vec(v1), NP);
            yp = linspace(aa_vec(v2), bb_vec(v2), NP);
            
            [XP, YP] = meshgrid(xp, yp);
            out_x{ii} = XP;
            out_y{ii} = YP;
            
            nb_pts = length(xp) * length(yp);
           
            PTS = avg_vec' * ones(1, nb_pts);
            
            % replace lines of non-constant directions
            PTS(v1, :) = XP(:)';
            PTS(v2, :) = YP(:)';
            
            % interpolate on sparse grid
            f_interp_eval = f_interp(PTS);
            
            % reshape to use surf
            FIP = reshape(f_interp_eval, size(XP));
            out_FIP{ii} = FIP; 
        end 

    else 
        % ===================================== %
        %      caso N > 3 PARI
        % ===================================== %
        % value of 'plot_grid_size'
        NP = 20;
        
        % extract info on lower and upper ends of each direction
        aa_vec = domain(1,:);
        bb_vec = domain(2,:);
        avg_vec = (aa_vec + bb_vec) / 2;
        
        % wrap interpolate on sparse grid into a @-function for ease of plotting
        f_interp = @(x) interpolate_on_sparse_grid(S, Sr, f_values, x);
        
        % initialization
        CUTS = floor(N / 2);
        V1 = zeros(1, CUTS);
        V2 = zeros(1, CUTS);
        out_x = cell(1, CUTS);
        out_y = cell(1, CUTS);
        out_FIP = cell(1, CUTS);

        for ii = 1:CUTS
            
            couple_loc = [2 * ii - 1, 2 * ii];
            v1 = couple_loc(1);
            v2 = couple_loc(2);
            V1(ii) = v1;
            V2(ii) = v2;
            
            % generate a mesh grid over the cut
            xp = linspace(aa_vec(v1), bb_vec(v1), NP);
            yp = linspace(aa_vec(v2), bb_vec(v2), NP);
            
            [XP, YP] = meshgrid(xp, yp);
            out_x{ii} = XP;
            out_y{ii} = YP;
            
            nb_pts = length(xp) * length(yp);
           
            PTS = avg_vec' * ones(1, nb_pts);
            
            % replace lines of non-constant directions
            PTS(v1, :) = XP(:)';
            PTS(v2, :) = YP(:)';
            
            % interpolate on sparse grid
            f_interp_eval = f_interp(PTS);
            
            % reshape to use surf
            FIP = reshape(f_interp_eval, size(XP));
            out_FIP{ii} = FIP; 
            
        end

    % chiusura if    
    end              

% chiusura della funzione principale
end