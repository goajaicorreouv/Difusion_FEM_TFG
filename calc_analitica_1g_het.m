function [phi_1, phi_2, phi_3] = calc_analitica_1g_het(x_eval)
% CALC_ANALITICA_1G_HET Calcula los 3 primeros modos espaciales de la solución 
% analítica teórica para un reactor de geometría de placa heterogéneo (1 Grupo).
% 
% Uso: 
%   [phi_1, phi_2, phi_3] = calc_analitica_1g_het(x_eval);
% Donde x_eval es un vector columna con las posiciones dadas en [0, L] cm.

    %% 1. Definición de Propiedades y Geometría
    Dc = 0.776;
    Sigma_ac = 0.0244;
    nu_Sigma_f = 0.0260;
    
    Dr = 1.446;
    Sigma_ar = 0.0077;
    
    a = 150.0;  % Semiespesor núcleo (cm)
    b = 25.0;   % Espesor reflectores (cm)
    
    %% 2. Parámetros dependientes del reflector
    kappa_r = sqrt(Sigma_ar / Dr);
    coth_kb = 1.0 / tanh(kappa_r * b);
    rhs = Dr * kappa_r * coth_kb; 
    
    %% 3. Ecuaciones de Criticidad
    f_even = @(Bc) Dc .* Bc .* tan(Bc * a) - rhs;
    f_odd  = @(Bc) -Dc .* Bc .* cot(Bc * a) - rhs;
    
    %% 4. Root Finding (Barriendo valores de Buckling)
    Bc_max = 5 * pi / a; 
    Bc_scan = linspace(1e-8, Bc_max, 20000);
    
    even_roots = encuentra_raices(f_even, Bc_scan);
    odd_roots  = encuentra_raices(f_odd, Bc_scan);
    
    %% 5. Unir y ordenar modos de mayor a menor k_eff
    all_modes = [];
    for i = 1:length(even_roots)
        k = nu_Sigma_f / (Dc * even_roots(i)^2 + Sigma_ac);
        all_modes = [all_modes; 1, even_roots(i), k]; 
    end
    for i = 1:length(odd_roots)
        k = nu_Sigma_f / (Dc * odd_roots(i)^2 + Sigma_ac);
        all_modes = [all_modes; 2, odd_roots(i), k];
    end
    
    % Ordenar matriz decrecientemente según k_eff (columna 3)
    all_modes = sortrows(all_modes, -3);
    
    if size(all_modes, 1) < 3
        error('MATLAB:RootsNotFound', 'No se detectaron suficientes modos en el escaneo');
    end
    
    %% 6. Evaluar la función en la malla de coordenadas x_eval
    x_sym = x_eval - (a + b);
    phis = zeros(length(x_eval), 3);
    
    for m = 1:3
        paridad = all_modes(m, 1);
        Bc = all_modes(m, 2);
        
        A_core = 1.0;
        phi_modal = zeros(size(x_sym));
        
        for k = 1:length(x_sym)
            x_val = x_sym(k);
            
            if x_val < -a % Reflector izquierdo
                if paridad == 1 
                    C_R = A_core * cos(Bc * a) / sinh(kappa_r * b);
                    C_L = C_R; 
                else 
                    C_R = A_core * sin(Bc * a) / sinh(kappa_r * b);
                    C_L = -C_R;
                end
                phi_modal(k) = C_L * sinh(kappa_r * (x_val + a + b));
                
            elseif x_val > a % Reflector derecho
                if paridad == 1
                    C_R = A_core * cos(Bc * a) / sinh(kappa_r * b);
                else
                    C_R = A_core * sin(Bc * a) / sinh(kappa_r * b);
                end
                phi_modal(k) = C_R * sinh(kappa_r * (a + b - x_val));
                
            else % Núcleo
                if paridad == 1
                    phi_modal(k) = A_core * cos(Bc * x_val);
                else
                    phi_modal(k) = A_core * sin(Bc * x_val);
                end
            end
        end
        
        phis(:, m) = phi_modal / max(abs(phi_modal));
    end
    
    phi_1 = phis(:, 1);
    phi_2 = phis(:, 2);
    phi_3 = phis(:, 3);
end

%% Sub-función auxiliar puramente abstracta para las raíces (CORREGIDA)
function roots = encuentra_raices(func, grid_vals)
    roots = [];
    f_vals = func(grid_vals);
    
    for i = 1:(length(grid_vals)-1)
        v1 = f_vals(i);
        v2 = f_vals(i+1);
        
        % EL TRUCO MAGISTRAL: 
        % Las raíces físicas reales SIEMPRE cruzan de negativo a positivo.
        % Las asíntotas (falsas raíces) saltan de positivo a negativo.
        if (v1 < 0) && (v2 > 0)
            try
                root_val = fzero(func, [grid_vals(i), grid_vals(i+1)]);
                
                % Guardar solo si el residuo es verdaderamente ~ 0
                if abs(func(root_val)) < 1e-4
                    roots = [roots, root_val];
                end
            catch
                % Ignorar si fzero falla
            end
        end
    end
end