function [phi_1, phi_2, phi_3, k_effs] = calc_analitica_2g_het(x_eval)
% CALC_ANALITICA_2G_HET Solución analítica completa para el reactor placa
% heterogéneo de 2 GRUPOS (Reflector - Núcleo - Reflector).
%
% Devuelve:
%   phi_1, phi_2, phi_3 : Matrices de tamaño [length(x_eval), 2].
%                         Columna 1 = Grupo 1 (Rápido), Columna 2 = Grupo 2 (Térmico)
%   k_effs              : Vector con los 3 autovalores fundamentales k_eff calculados.

    % =========================================================================
    % 1. PROPIEDADES DE MATERIALES Y GEOMETRÍA (2 Grupos)
    % =========================================================================
    a = 150.0; % Semiespesor del núcleo
    b = 25.0;  % Espesor del reflector

    % REFLECTOR (Material 1)
    D1_r = 2.0;   D2_r = 0.3;
    Sa1_r =  0.00; Sa2_r = 0.01;
    Ss12_r = 0.04;
    
    % NÚCLEO COMBUSTIBLE (Material 2)
    D1_c = 1.446;   D2_c = 0.40;
    Sa1_c = 0.010; Sa2_c = 0.085;
    Ss12_c = 0.020;
    nSf1_c = 0.000; nSf2_c = 0.1350;

    % =========================================================================
    % 2. CONSTANTES DEL REFLECTOR (Material no multiplicante)
    % =========================================================================
    k1_r = sqrt((Sa1_r + Ss12_r) / D1_r);
    k2_r = sqrt(Sa2_r / D2_r);
    % Constante de acoplamiento espectral en el reflector
    S3_r = Ss12_r / (D2_r * (k2_r^2 - k1_r^2));

    % =========================================================================
    % 3. BÚSQUEDA DEL ESPECTRO DE AUTOVALORES (BARRIDO EN K_EFF)
    % =========================================================================
    k_scan = linspace(0.1, 1.8, 1500);
    
    k_evens = search_roots(@(k) eq_det(k, 1), k_scan);
    k_odds  = search_roots(@(k) eq_det(k, 2), k_scan);

    % Consolidar y ordenar [k_eff, paridad]
    todos_modos = [[k_evens, ones(size(k_evens))]; [k_odds, 2 * ones(size(k_odds))]];
    todos_modos = sortrows(todos_modos, -1); % Orden descendente (mayor a menor reactividad)

    if size(todos_modos, 1) < 3
        error('MATLAB:Analitico2G:FaltanRaices', 'Aumenta la densidad del grid del scan de K_eff.');
    end
    k_effs = todos_modos(1:3, 1);

    % =========================================================================
    % 4. RECONSTRUCCIÓN VECTORIAL DE LOS MODOS
    % =========================================================================
    x_s = x_eval - (a + b); % Centrado en el núcleo [-a-b, a+b] asumiendo L = 2*(a+b) y origen en 0
    phi_all = cell(3, 1);

    for m = 1:3
        k = todos_modos(m, 1);
        paridad = todos_modos(m, 2);
        
        S1_c_val = Sa1_c + Ss12_c - nSf1_c / k;
        A_p = 0.5 * (S1_c_val / D1_c + Sa2_c / D2_c);
        B_p = 0.5 * sqrt((S1_c_val / D1_c - Sa2_c / D2_c)^2 + 4 * Ss12_c * nSf2_c / (k * D1_c * D2_c));
        
        mu = sqrt(-A_p + B_p);
        nu = sqrt(A_p + B_p);
        
        S1 = Ss12_c / (D2_c * mu^2 + Sa2_c);
        S2 = Ss12_c / (-D2_c * nu^2 + Sa2_c);

        % Recreamos la matriz singular
        if paridad == 1
            colA = [cos(mu*a); -D1_c*mu*sin(mu*a); S1*cos(mu*a); -D2_c*S1*mu*sin(mu*a)];
            colC = [cosh(nu*a); D1_c*nu*sinh(nu*a); S2*cosh(nu*a); D2_c*S2*nu*sinh(nu*a)];
        else
            colA = [sin(mu*a); D1_c*mu*cos(mu*a); S1*sin(mu*a); D2_c*S1*mu*cos(mu*a)];
            colC = [sinh(nu*a); D1_c*nu*cosh(nu*a); S2*sinh(nu*a); D2_c*S2*nu*cosh(nu*a)];
        end
        
        colE = [-sinh(k1_r*b); D1_r*k1_r*cosh(k1_r*b); -S3_r*sinh(k1_r*b); D2_r*S3_r*k1_r*cosh(k1_r*b)];
        colG = [0; 0; -sinh(k2_r*b); D2_r*k2_r*cosh(k2_r*b)];

        % Extraer el vector del kernel puramente analítico usando Descomposición de Valores Singulares (SVD)
        M = [colA, colC, colE, colG];
        [~, ~, V] = svd(M);
        consts = V(:, end); % El autovector correspondiente al "cero" numérico

        g1 = zeros(size(x_s));
        g2 = zeros(size(x_s));

        for j = 1:length(x_s)
            xv = x_s(j);
            if abs(xv) <= a % Zona del núcleo
                if paridad == 1
                    g1(j) = consts(1)*cos(mu*xv) + consts(2)*cosh(nu*xv);
                    g2(j) = consts(1)*S1*cos(mu*xv) + consts(2)*S2*cosh(nu*xv);
                else
                    g1(j) = consts(1)*sin(mu*xv) + consts(2)*sinh(nu*xv);
                    g2(j) = consts(1)*S1*sin(mu*xv) + consts(2)*S2*sinh(nu*xv);
                end
            else % Reflector exterior
                dr = a + b - abs(xv);
                sig = 1;
                if xv < -a && paridad == 2
                    sig = -1; % Parche antisimetría espacial
                end 
                
                g1(j) = sig * consts(3)*sinh(k1_r*dr);
                g2(j) = sig * (consts(3)*S3_r*sinh(k1_r*dr) + consts(4)*sinh(k2_r*dr));
            end
        end

        % Normalización dual respetando el ratio térmico/rápido intrínseco
        max_flujo = max(max(abs(g1)), max(abs(g2)));
        phi_all{m} = [g1 / max_flujo, g2 / max_flujo];
    end

    phi_1 = phi_all{1};
    phi_2 = phi_all{2};
    phi_3 = phi_all{3};


    % =========================================================================
    % FUNCIONES ANIDADAS (Tienen acceso a las variables principales)
    % =========================================================================
    
    function det_val = eq_det(k, paridad)
        % Coeficientes de absorción efectivos en núcleo (añadiendo fisión)
        S1_c_val = Sa1_c + Ss12_c - nSf1_c / k;
        
        A_p = 0.5 * (S1_c_val / D1_c + Sa2_c / D2_c);
        B_p = 0.5 * sqrt((S1_c_val / D1_c - Sa2_c / D2_c)^2 + 4 * Ss12_c * nSf2_c / (k * D1_c * D2_c));
        
        mu_sq = -A_p + B_p;
        nu_sq = A_p + B_p;
        
        % Si ocurre algún k que lance raíces complejas espurias, lo vetamos
        if mu_sq <= 0 || nu_sq <= 0
            det_val = 1e10;
            return;
        end
        
        mu_val = sqrt(mu_sq);
        nu_val = sqrt(nu_sq);
        
        % Coeficientes de acoplamiento modal en el núcleo
        S1_val = Ss12_c / (D2_c * mu_val^2 + Sa2_c);
        S2_val = Ss12_c / (-D2_c * nu_val^2 + Sa2_c);
        
        % Construcción de las 4 columnas de la matriz global de continuidad
        % Orden de filas: [Flux G1; Curr G1; Flux G2; Curr G2]
        if paridad == 1 % Modos Pares (Cosenos)
            colA_val = [cos(mu_val*a); -D1_c*mu_val*sin(mu_val*a); S1_val*cos(mu_val*a); -D2_c*S1_val*mu_val*sin(mu_val*a)];
            colC_val = [cosh(nu_val*a); D1_c*nu_val*sinh(nu_val*a); S2_val*cosh(nu_val*a); D2_c*S2_val*nu_val*sinh(nu_val*a)];
        else % Modos Impares (Senos)
            colA_val = [sin(mu_val*a); D1_c*mu_val*cos(mu_val*a); S1_val*sin(mu_val*a); D2_c*S1_val*mu_val*cos(mu_val*a)];
            colC_val = [sinh(nu_val*a); D1_c*nu_val*cosh(nu_val*a); S2_val*sinh(nu_val*a); D2_c*S2_val*nu_val*cosh(nu_val*a)];
        end
        
        % El decaimiento del reflector requiere continuidad en a
        colE_val = [-sinh(k1_r*b); D1_r*k1_r*cosh(k1_r*b); -S3_r*sinh(k1_r*b); D2_r*S3_r*k1_r*cosh(k1_r*b)];
        colG_val = [0; 0; -sinh(k2_r*b); D2_r*k2_r*cosh(k2_r*b)];
        
        % Calcular determinante analítico
        det_val = det([colA_val, colC_val, colE_val, colG_val]);
    end

    function roots = search_roots(func, grid_vals)
        roots = [];
        f_vals = zeros(length(grid_vals), 1);
        for i = 1:length(grid_vals)
            f_vals(i) = func(grid_vals(i));
        end
        
        for i = 1:(length(grid_vals)-1)
            v1 = f_vals(i);
            v2 = f_vals(i+1);
            if v1 * v2 <= 0 && isfinite(v1) && isfinite(v2)
                try
                    r = fzero(func, [grid_vals(i), grid_vals(i+1)]);
                    roots = [roots; r];
                catch
                    % Falla silenciosamente si fzero no puede converger en el intervalo
                end
            end
        end
    end

end