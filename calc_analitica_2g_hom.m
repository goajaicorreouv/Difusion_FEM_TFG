function [phi_1, phi_2, phi_3, k_effs] = calc_analitica_2g_hom(x_eval)
% CALC_ANALITICA_2G_HOM Solución analítica completa para el reactor placa
% HOMOGÉNEO de 2 GRUPOS (Longitud total L = 350 cm).
%
% Devuelve:
%   phi_1, phi_2, phi_3 : Matrices de tamaño [length(x_eval), 2].
%                         Columna 1 = Rápido, Columna 2 = Térmico
%   k_effs              : Vector con los 3 autovalores fundamentales k_eff.

    %% 1. PROPIEDADES DEL MATERIAL ÚNICO Y GEOMETRÍA
    L = 350.0;  % Contemplado el tamaño completo de 350 cm
    
    D1  = 1.446;   D2  = 0.40;
    Sa1 = 0.010;  Sa2 = 0.0850;
    Ss12 = 0.020;
    nSf1 = 0.000; nSf2 = 0.1350;


    %% 2. CÁLCULO ANALÍTICO DEL k_eff PARA CADA MODO
    n_modos = 3;
    k_effs = zeros(n_modos, 1);
    ratios = zeros(n_modos, 1); % Ratio de acoplamiento A2/A1
    
    for n = 1:n_modos
        Bn = n * pi / L;
        Bn2 = Bn^2;
        
        % Relación de dispersión acoplada de 2 grupos
        denom_fast  = D1 * Bn2 + Sa1 + Ss12;
        denom_therm = D2 * Bn2 + Sa2;
        
        k_effs(n) = (nSf1 + nSf2 * Ss12 / denom_therm) / denom_fast;
        ratios(n) = Ss12 / denom_therm;  
    end

    %% 3. RECONSTRUCCIÓN DE LOS PERFILES ESPACIALES
    x_eval = x_eval(:);
    phi_all = cell(n_modos, 1);
    
    for n = 1:n_modos
        % Grupo 1 (Rápido): sin(Bn * x)
        g1 = sin(n * pi * x_eval / L);
        
        % Grupo 2 (Térmico): ratios(n) * sin(Bn * x)
        g2 = ratios(n) * g1;
        
        % Normalización global conjunta (conserva la proporción física real G1/G2)
        max_flujo = max(max(abs(g1)), max(abs(g2)));
        
        % Guardamos como matriz de dos columnas
        phi_all{n} = [g1 / max_flujo, g2 / max_flujo];
    end

    %% 4. ASIGNACIÓN DE SALIDAS
    phi_1 = phi_all{1};
    phi_2 = phi_all{2};
    phi_3 = phi_all{3};
end