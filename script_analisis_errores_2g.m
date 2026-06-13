% DESCRIPCIÓN:
% Script centralizado para el cálculo y ploteo del Error Residual FEM frente a la Analítica
% PARA 2 GRUPOS DE ENERGÍA. 
% Preserva el ratio físico de los flujos mediante normalización global.

clearvars;
close all; 

% Forzar LaTeX globalmente en todas las gráficas
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');

% =========================================================================
% === ZONA DE CONFIGURACIÓN DEL USUARIO ===================================
% =========================================================================
% Tipo de simulación (debe ser de 2 grupos: '2g_homo' o '2g_heter')
tipo_caso = '2g_heter'; 

% Lista de archivos generados por h-refinamiento (Asegúrate de que sean de 2G)
archivos_convergencia = {'FEM_2G_HET_N14_p4.mat', ...
                         'FEM_2G_HET_N28_p4.mat', ...
                         'FEM_2G_HET_N56_p4.mat', ...
                         'FEM_2G_HET_N112_p4.mat'}; 

% El archivo base al que le queremos examinar el error espacial Absoluto
archivo_error_local = 'FEM_2G_HET_N28_p4.mat';
% =========================================================================

%% --- PREPARACIÓN DE FIGURAS ---
fig_rmse_g1 = figure('Name', 'RMSE - Grupo 1 (Rápido)', 'Position', [100 500 1200 400]);
fig_rmse_g2 = figure('Name', 'RMSE - Grupo 2 (Térmico)', 'Position', [100 50 1200 400]);

%% --- BUCLE PRINCIPAL DE CONVERGENCIA (RMSE) ---
for modo_i = 1:3
    ndofs_list = zeros(1, length(archivos_convergencia));
    rmse_g1_list = zeros(1, length(archivos_convergencia));
    rmse_g2_list = zeros(1, length(archivos_convergencia));
    
    for i = 1:length(archivos_convergencia)
        try
            data = load(archivos_convergencia{i}, 'malla', 'problema');
        catch
            if modo_i == 1, warning('No se encontró: %s', archivos_convergencia{i}); end
            continue;
        end
        n_pts = length(data.malla.x_nodos);
        
        % 1. Extracción FEM (Vector completo acoplado)
        if isprop(data.problema, 'phi_nodos') && size(data.problema.phi_nodos, 2) >= modo_i 
            phi_fem_full = data.problema.phi_nodos(:, modo_i);
        else 
            phi_fem_full = data.problema.phi_inc(:, modo_i);
        end 
        
        phi_fem_g1 = phi_fem_full(1:n_pts);
        phi_fem_g2 = phi_fem_full(n_pts+1 : 2*n_pts);
        
        % 2. Extracción Analítica (Ambos grupos)
        [phi_ex_g1, phi_ex_g2] = eval_analytic_on_nodes_2g(data.malla.x_nodos, data.malla.L, tipo_caso, modo_i);
        
        % 3. NORMALIZACIÓN GLOBAL ACORDE A LA FÍSICA
        max_fem = max(max(abs(phi_fem_g1)), max(abs(phi_fem_g2)));
        max_ex  = max(max(abs(phi_ex_g1)), max(abs(phi_ex_g2)));
        
        phi_fem_g1 = phi_fem_g1 / max_fem;
        phi_fem_g2 = phi_fem_g2 / max_fem;
        phi_ex_g1  = phi_ex_g1 / max_ex;
        phi_ex_g2  = phi_ex_g2 / max_ex;
        
        % 4. ALINEACIÓN DE SIGNOS (Usando G1 como ancla)
        idx_ref = max(1, floor(length(phi_fem_g1) / 5));
        if sign(phi_fem_g1(idx_ref)) ~= sign(phi_ex_g1(idx_ref)) && phi_ex_g1(idx_ref) ~= 0
            phi_fem_g1 = -phi_fem_g1;
            phi_fem_g2 = -phi_fem_g2; % G2 también debe invertirse
        end
        
        % 5. CÁLCULO DE ERRORES RMSE
        ndofs_list(i)   = n_pts; % Se puede usar length(phi_fem_full) si prefieres los dofs totales
        rmse_g1_list(i) = sqrt(mean((phi_fem_g1 - phi_ex_g1).^2));
        rmse_g2_list(i) = sqrt(mean((phi_fem_g2 - phi_ex_g2).^2));
    end
    
    % Filtrar datos válidos
    valid_idx = ndofs_list > 0;
    ndofs_list = ndofs_list(valid_idx);
    
    % PLOT GRUPO 1
    figure(fig_rmse_g1); subplot(1, 3, modo_i); hold on; grid on; box on; set(gca, 'XScale', 'log', 'YScale', 'log');
    if length(ndofs_list) >= 2
        plot(ndofs_list, rmse_g1_list(valid_idx), 'b-o', 'LineWidth', 2, 'MarkerFaceColor', 'b');
        p_poly = polyfit(log(ndofs_list), log(rmse_g1_list(valid_idx)), 1);
        text(min(ndofs_list)*1.2, max(rmse_g1_list(valid_idx))*0.5, sprintf('Pendiente $m \\approx %.2f$', p_poly(1)), 'Interpreter', 'latex', 'FontSize', 12, 'BackgroundColor', 'w', 'EdgeColor', 'k');
    end
    xlabel('$N_{nodos}$', 'Interpreter', 'latex', 'FontSize', 13);
    if modo_i == 1, ylabel('RMSE Rápido', 'Interpreter', 'latex', 'FontSize', 13); end
    title(sprintf('\\textbf{Modo %d - G1}', modo_i), 'Interpreter', 'latex');
    
    % PLOT GRUPO 2
    figure(fig_rmse_g2); subplot(1, 3, modo_i); hold on; grid on; box on; set(gca, 'XScale', 'log', 'YScale', 'log');
    if length(ndofs_list) >= 2
        plot(ndofs_list, rmse_g2_list(valid_idx), 'r-s', 'LineWidth', 2, 'MarkerFaceColor', 'r');
        p_poly2 = polyfit(log(ndofs_list), log(rmse_g2_list(valid_idx)), 1);
        text(min(ndofs_list)*1.2, max(rmse_g2_list(valid_idx))*0.5, sprintf('Pendiente $m \\approx %.2f$', p_poly2(1)), 'Interpreter', 'latex', 'FontSize', 12, 'BackgroundColor', 'w', 'EdgeColor', 'k');
    end
    xlabel('$N_{nodos}$', 'Interpreter', 'latex', 'FontSize', 13);
    if modo_i == 1, ylabel('RMSE Térmico', 'Interpreter', 'latex', 'FontSize', 13); end
    title(sprintf('\\textbf{Modo %d - G2}', modo_i), 'Interpreter', 'latex');
end


%% --- GRÁFICA 3 y 4: ERROR ABSOLUTO ESPACIAL ---
if exist(archivo_error_local, 'file')
    fig_err_g1 = figure('Name', 'Error Absoluto - Grupo 1', 'Position', [100 500 1200 400]);
    fig_err_g2 = figure('Name', 'Error Absoluto - Grupo 2', 'Position', [100 50 1200 400]);
    
    data_loc = load(archivo_error_local, 'malla', 'problema');
    n_pts = length(data_loc.malla.x_nodos);
    
    for modo_i = 1:3
        % Extracción FEM y separación
        if isprop(data_loc.problema, 'phi_nodos')
            phi_fem_full = data_loc.problema.phi_nodos(:, modo_i);
        else 
            phi_fem_full = data_loc.problema.phi_inc(:, modo_i);
        end 
        phi_fem_g1 = phi_fem_full(1:n_pts);
        phi_fem_g2 = phi_fem_full(n_pts+1 : 2*n_pts);
        
        % Extracción Analítica
        [phi_ex_g1, phi_ex_g2] = eval_analytic_on_nodes_2g(data_loc.malla.x_nodos, data_loc.malla.L, tipo_caso, modo_i);
        
        % Normalización global
        max_fem = max(max(abs(phi_fem_g1)), max(abs(phi_fem_g2)));
        max_ex  = max(max(abs(phi_ex_g1)), max(abs(phi_ex_g2)));
        phi_fem_g1 = phi_fem_g1 / max_fem; phi_fem_g2 = phi_fem_g2 / max_fem;
        phi_ex_g1  = phi_ex_g1 / max_ex;   phi_ex_g2  = phi_ex_g2 / max_ex;
        
        % Alineación de signos
        idx_ref = max(1, floor(length(phi_fem_g1) / 5));
        if sign(phi_fem_g1(idx_ref)) ~= sign(phi_ex_g1(idx_ref)) && phi_ex_g1(idx_ref) ~= 0
            phi_fem_g1 = -phi_fem_g1; phi_fem_g2 = -phi_fem_g2;
        end
        
        % Cálculo de error absoluto
        error_abs_g1 = abs(phi_ex_g1 - phi_fem_g1);
        error_abs_g2 = abs(phi_ex_g2 - phi_fem_g2);
        
        p_str = regexp(archivo_error_local, 'p(\d+)', 'tokens');
        if ~isempty(p_str), titulo_p = sprintf(' ($p=%s$)', p_str{1}{1}); else, titulo_p = ''; end
        
        % Plot Local G1
        figure(fig_err_g1); subplot(1, 3, modo_i); hold on; grid on; box on; set(gca, 'YScale', 'log');
        plot(data_loc.malla.x_nodos, error_abs_g1 + eps, 'b-', 'LineWidth', 1.5);
        xlim([0, data_loc.malla.L]);
        if contains(tipo_caso, 'heter')
            xline(25, 'k--', 'LineWidth', 0.8, 'HandleVisibility', 'off'); xline(data_loc.malla.L - 25, 'k--', 'LineWidth', 0.8, 'HandleVisibility', 'off');
        end
        xlabel('Posición $x$ (cm)', 'Interpreter', 'latex', 'FontSize', 12);
        if modo_i == 1, ylabel('$|\phi_{anal} - \phi_{fem}|$ (G1)', 'Interpreter', 'latex', 'FontSize', 12); end
        title(sprintf('\\textbf{Modo %d%s - G1}', modo_i, titulo_p), 'Interpreter', 'latex');
        
        % Plot Local G2
        figure(fig_err_g2); subplot(1, 3, modo_i); hold on; grid on; box on; set(gca, 'YScale', 'log');
        plot(data_loc.malla.x_nodos, error_abs_g2 + eps, 'r-', 'LineWidth', 1.5);
        xlim([0, data_loc.malla.L]);
        if contains(tipo_caso, 'heter')
            xline(25, 'k--', 'LineWidth', 0.8, 'HandleVisibility', 'off'); xline(data_loc.malla.L - 25, 'k--', 'LineWidth', 0.8, 'HandleVisibility', 'off');
        end
        xlabel('Posición $x$ (cm)', 'Interpreter', 'latex', 'FontSize', 12);
        if modo_i == 1, ylabel('$|\phi_{anal} - \phi_{fem}|$ (G2)', 'Interpreter', 'latex', 'FontSize', 12); end
        title(sprintf('\\textbf{Modo %d%s - G2}', modo_i, titulo_p), 'Interpreter', 'latex');
    end
else
    disp('Archivo Local No Encontrado. No se generan las gráficas de Error Local.');
end

% =========================================================================
% FUNCIÓN INTERNA: INYECTOR DINÁMICO ANALÍTICO PARA 2 GRUPOS
% =========================================================================
function [phi_ex_g1, phi_ex_g2] = eval_analytic_on_nodes_2g(x_nodos, L, tipo_caso, modo_i)
    x_nodos = x_nodos(:); % Forzar columna
    
    if strcmp(tipo_caso, '2g_homo')
        [phi_a1, phi_a2, phi_a3, ~] = calc_analitica_2g_hom(x_nodos);
        switch modo_i
            case 1, phi_ex_g1 = phi_a1(:,1); phi_ex_g2 = phi_a1(:,2);
            case 2, phi_ex_g1 = phi_a2(:,1); phi_ex_g2 = phi_a2(:,2);
            case 3, phi_ex_g1 = phi_a3(:,1); phi_ex_g2 = phi_a3(:,2);
            otherwise, error('Modo no soportado.');
        end
        
    elseif strcmp(tipo_caso, '2g_heter')
        [phi_a1, phi_a2, phi_a3, ~] = calc_analitica_2g_het(x_nodos);
        switch modo_i
            case 1, phi_ex_g1 = phi_a1(:,1); phi_ex_g2 = phi_a1(:,2);
            case 2, phi_ex_g1 = phi_a2(:,1); phi_ex_g2 = phi_a2(:,2);
            case 3, phi_ex_g1 = phi_a3(:,1); phi_ex_g2 = phi_a3(:,2);
            otherwise, error('Modo no soportado.');
        end
        
    else
        error('Este script de validación residual está adaptado exclusivamente para 2 grupos (2g_homo o 2g_heter).');
    end
end