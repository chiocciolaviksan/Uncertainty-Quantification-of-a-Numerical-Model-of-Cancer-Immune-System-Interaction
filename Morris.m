%% ------------------------------------------------------------------- %%
%                    Script that implements the 
%                   CALCULATION OF MORRIS INDICES.
% -------------------------------------------------------------------- %%

function [emme,mi,sigma] = Morris(S,Sr,I_on_Sr,domain,K)
    
   
    % number of variables to consider
    N = size(domain,2);                  
    
    % --------------------------------- % 
    %     Step delta_i 
    % --------------------------------- % 
    % Note that I may have different domains in x and y, and thus 
    % different discretization steps in x and y.
    d = zeros(1,N);
    for i = 1:N
            d(i) = ( domain(2,i)-domain(1,i) ) / 1E5;
            %d(i) = d(i).*1000;
    end
    
    % --------------------------------- % 
    %     Random Points P_k 
    % --------------------------------- % 
    % i.e., K random values of N elements within the given domain. 
    % Note that I may have two different domains in x and y.
    P = repmat(domain(1,:), [K 1]);  
    P = P + rand(K,N).*repmat(domain(2,:)-domain(1,:), [K,1]);
   
    % plot(P(:,1),P(:,2),'*')

    % --------------------------------- %
    %      perturbated points P_k 
    % --------------------------------- %
    Pert = zeros(K,1);
    P_pert = P;         
    P_di{1} = P; 

    for i = 1:N
        P_pert(:,i) = P_pert(:,i)+d(i);
        P_di{i+1} = P_pert;
    end 

    
    % EVAL_POINTS are the points at which the derivative must be evaluated. 
    % It is a matrix with points stored as columns.
    P = P';
   
    % --------------------------------- %
    %      EE_i calculation
    % --------------------------------- %
    for i = 2:(N+1)
        % ciclo sulle N variabili
        f_pk_1 = interpolate_on_sparse_grid(S,Sr,I_on_Sr,P_di{i-1}');
        f_pk = interpolate_on_sparse_grid(S,Sr,I_on_Sr,P_di{i}');
        
        if size(f_pk,1)==1
            % static case
            EE_i(i-1,:) = ( f_pk - f_pk_1 ) / d(i-1);
            emme(i-1,:) = 1/K*sum(abs(EE_i(i-1,:)))*(domain(1,i-1)+domain(2,i-1));     
            mi(i-1,:) = 1/K*sum(EE_i(i-1,:))*(domain(1,i-1)+domain(2,i-1));            

            % sigma(i-1,:) = 1/K*sum((abs(abs(EE_i(i-1,:))-mi(i-1,:)))^2);
            sigma(i-1,:) = 1/K*sum( (EE_i(i-1,:)-mi(i-1,:)).^2 );
         
        else
            % time dependent
            T = size(f_pk,1);   

            for t = 1:T
                EE_i(i-1,:) = ( f_pk(t,:) - f_pk_1(t,:) ) / d(i-1);
                emme(i-1,t) = 1/K*sum(abs(EE_i(i-1,:)))*(domain(1,i-1)+domain(2,i-1));              
                mi(i-1,t) = 1/K*sum(EE_i(i-1,:))*(domain(1,i-1)+domain(2,i-1));                     
                
                % sigma(i-1,t) = 1/K*sum((abs(EE_i(i-1,:)-mi(i-1,t))).^2);
                sigma(i-1,t) = 1/K*sum( (EE_i(i-1,:)-mi(i-1,t)).^2 );
            end 
            

        end 

    end
    
    % ------------------------------------------------- %
    %      GRAFICAL representation of the concept
    % ------------------------------------------------- %
    if N==2 

        z = zeros(K,1);
        figure('Name','Controllo')
        plot3(P_di{1,1}(:,1),P_di{1,1}(:,2),z,'*')                  % original points
        hold on
        plot3(P_di{1,2}(:,1),P_di{1,2}(:,2),z,'mo')                 % perturbed x
        hold on 
        plot3(P_di{1,3}(:,1),P_di{1,3}(:,2),z,'r*')                 % perturbed x e y
        hold on
        plot3(P_di{1,2}(:,1),P_di{1,2}(:,2),f_pk_1(5,:),'mo')       % f --  pert in x
        hold on
        plot3(P_di{1,3}(:,1),P_di{1,3}(:,2),f_pk(5,:),'r*')         % f |  pert in x e y
        hold on
        J = 10;       
        for i=1:J
            % pert in x
            xl = P_di{1,2}(i,1)*ones(K,1);
            yl = P_di{1,2}(i,2)*ones(K,1);
            zl = linspace(0,f_pk_1(5,i),K);
            plot3(xl,yl,zl,'-c')
            hold on
            % pert in x e y
            xl2 = P_di{1,3}(i,1)*ones(K,1);
            yl2 = P_di{1,3}(i,2)*ones(K,1);
            zl2 = linspace(0,f_pk(5,i),K);
            plot3(xl2,yl2,zl2,'-g')
            hold on
        end 
        grid on
        xlabel('y1')
        ylabel('y2')
        zlabel('QoI')
        legend('punti originali','pert. in x','pert. in x e y',...
            'f(pert in x)','f(pert in x e y)','proiezione','proiezione')
        title('Controllo')
    
    end

   
    
% end of function 
end 

