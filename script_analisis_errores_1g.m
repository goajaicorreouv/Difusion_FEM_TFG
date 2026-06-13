% DESCRIPCIÓN:
% Script centralizado para el cálculo y ploteo del Error Residual FEM frente a la Analítica.
% Permite evaluar el RMSE (Root Mean Square Error) en log-log a lo largo de un p-refinamiento
% y graficar el error absoluto local para visualizar dónde sufre más la malla.

clearvars;
close all; 

% Forzar LaTeX globalmente en todas las gráficas
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');

% =========================================================================
% === ZONA DE CONFIGURACIÓN DEL USUARIO ===================================
% =========================================================================

% Tipo de simulación (elige: '1g_homo', '1g_heter', '2g_homo', '2g_heter')
tipo_caso = '1g_heter'; 

% Lista de archivos generados por p-refinamiento para medir la pendiente (Ndofs)
archivos_convergencia = {'FEM_1G_HET_N14_p5.mat', ...
                         'FEM_1G_HET_N28_p5.mat', ...
                         'FEM_1G_HET_N56_p5.mat','FEM_1G_HET_N112_p5.mat'}; 

% El archivo base al que le queremos examinar el error espacial Absoluto
archivo_error_local = 'FEM_1G_HET_N28_p5.mat';

% El script evaluará automáticamente los 3 primeros modos (Fundamental y 2 armónicos).

% =========================================================================
% =========================================================================

%% --- GRÁFICA 1: LOG-LOG RMSE vs NDOFS ---
figure('Name', 'Convergencia Espacial (RMSE)', 'Position', [100 200 1200 450]);

for modo_i = 1:3
    subplot(1, 3, modo_i);
    hold on; grid on; box on;
    set(gca, 'XScale', 'log', 'YScale', 'log');

    ndofs_list = zeros(1, length(archivos_convergencia));
    rmse_list = zeros(1, length(archivos_convergencia));

    for i = 1:length(archivos_convergencia)
        try
            data = load(archivos_convergencia{i}, 'malla', 'problema');
        catch
            if modo_i == 1
                warning('MATLAB:NoFile', 'No se ha encontrado el archivo %s.', archivos_convergencia{i});
            end
            continue;
        end
        n_pts = length(data.malla.x_nodos);
        
        if isprop(data.problema, 'phi_nodos') && size(data.problema.phi_nodos, 2) >= modo_i 
            phi_fem = data.problema.phi_nodos(:, modo_i);
        else 
            phi_fem = data.problema.phi_inc(1:n_pts, modo_i);
        end 
        
        if contains(tipo_caso, '2g')
            phi_fem = phi_fem(1:n_pts);
        end
        
        phi_exact = eval_analytic_on_nodes(data.malla.x_nodos, data.malla.L, tipo_caso, modo_i);
        
        phi_fem = phi_fem / max(abs(phi_fem));
        phi_exact = phi_exact / max(abs(phi_exact));
        
        idx_ref = max(1, floor(length(phi_fem) / 5));
        if sign(phi_fem(idx_ref)) ~= sign(phi_exact(idx_ref))
            phi_fem = -phi_fem;
        end
        
        ndofs_list(i) = length(phi_fem); 
        rmse_list(i) = sqrt(mean((phi_fem - phi_exact).^2));
    end

    valid_idx = ndofs_list > 0;
    ndofs_list = ndofs_list(valid_idx);
    rmse_list = rmse_list(valid_idx);

    if length(ndofs_list) >= 2
        plot(ndofs_list, rmse_list, 'b-o', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'b');
        p_poly = polyfit(log(ndofs_list), log(rmse_list), 1);
        pendiente = p_poly(1);
        
        % Posicionado dinámico para no estorbar mucho la traza
        text(min(ndofs_list)*1.2, max(rmse_list)*0.5, ...
             sprintf('Pendiente $m \\approx %.2f$', pendiente), ...
             'Interpreter', 'latex', 'FontSize', 12, 'BackgroundColor', 'w', 'EdgeColor', 'k');
         
        xlabel('$N_{dofs}$', 'Interpreter', 'latex', 'FontSize', 13);
        if modo_i == 1, ylabel('RMSE', 'Interpreter', 'latex', 'FontSize', 13); end
        title(sprintf('\\textbf{Modo %d}', modo_i), 'Interpreter', 'latex', 'FontSize', 14);
    else
        title(sprintf('\\textbf{Modo %d (No hay datos)}', modo_i), 'Interpreter', 'latex', 'FontSize', 12);
    end
    hold off;
end

%% --- GRÁFICA 2: ERROR ABSOLUTO ESPACIAL ---
figure('Name', 'Distribución Local del Error Absoluto', 'Position', [100 100 1200 450]);

if exist(archivo_error_local, 'file')
    data_loc = load(archivo_error_local, 'malla', 'problema');
    n_pts = length(data_loc.malla.x_nodos);

    for modo_i = 1:3
        subplot(1, 3, modo_i);
        hold on; grid on; box on;
        set(gca, 'YScale', 'log');

        if isprop(data_loc.problema, 'phi_nodos') && size(data_loc.problema.phi_nodos, 2) >= modo_i 
            phi_fem_loc = data_loc.problema.phi_nodos(:, modo_i);
        else 
            phi_fem_loc = data_loc.problema.phi_inc(1:n_pts, modo_i);
        end 
        if contains(tipo_caso, '2g')
            phi_fem_loc = phi_fem_loc(1:n_pts);
        end

        phi_exact_loc = eval_analytic_on_nodes(data_loc.malla.x_nodos, data_loc.malla.L, tipo_caso, modo_i);

        phi_fem_loc = phi_fem_loc / max(abs(phi_fem_loc));
        phi_exact_loc = phi_exact_loc / max(abs(phi_exact_loc));
        
        idx_ref = max(1, floor(length(phi_fem_loc) / 5));
        if sign(phi_fem_loc(idx_ref)) ~= sign(phi_exact_loc(idx_ref))
            phi_fem_loc = -phi_fem_loc;
        end

        error_abs = abs(phi_exact_loc - phi_fem_loc);

        plot(data_loc.malla.x_nodos, error_abs + eps, 'r-', 'LineWidth', 1.5);
        xlim([0, data_loc.malla.L]);

        if contains(tipo_caso, 'heter')
            xline(25, 'k--', 'LineWidth', 0.8, 'HandleVisibility', 'off');
            xline(data_loc.malla.L - 25, 'k--', 'LineWidth', 0.8, 'HandleVisibility', 'off');
        end

        xlabel('Posición $x$ (cm)', 'Interpreter', 'latex', 'FontSize', 12);
        if modo_i == 1
            ylabel('$|\phi_{analitico} - \phi_{fem}|$ (Log)', 'Interpreter', 'latex', 'FontSize', 12);
        end
        
        p_str = regexp(archivo_error_local, 'p(\d+)', 'tokens');
        if ~isempty(p_str), titulo_p = sprintf(' ($p=%s$)', p_str{1}{1}); else, titulo_p = ''; end
        
        title(sprintf('\\textbf{Error Absoluto - Modo %d%s}', modo_i, titulo_p), 'Interpreter', 'latex');
        hold off;
    end
else
    disp('Archivo Local No Encontrado. No se genera la Gráfica 2.');
end


% =========================================================================
% FUNCIÓN INTERNA: INYECTOR DINÁMICO ANALÍTICO
% =========================================================================
function phi_exact = eval_analytic_on_nodes(x_nodos, L, tipo_caso, modo_i)
    x_nodos = x_nodos(:); % Forzar columna (Malla1D devuelve fila)
    if contains(tipo_caso, 'homo')
        % Exact homogeneous sine wave for both groups
        phi_exact = sin(modo_i * pi * x_nodos / L);
        
    elseif strcmp(tipo_caso, '1g_heter')
        [phi_a1, phi_a2, phi_a3] = calc_analitica_1g_het(x_nodos);
        switch modo_i
            case 1, phi_exact = phi_a1;
            case 2, phi_exact = phi_a2;
            case 3, phi_exact = phi_a3;
            otherwise, error('Modo no soportado (analítica limitada a modos 1-3).');
        end
        
    elseif strcmp(tipo_caso, '2g_heter')
        [phi_a1, phi_a2, phi_a3, ~] = calc_analitica_2g_het(x_nodos);
        % Tomaremos el grupo r?pido (col 1) frente al FEM cortado.
        switch modo_i
            case 1, phi_exact = phi_a1(:,1);
            case 2, phi_exact = phi_a2(:,1);
            case 3, phi_exact = phi_a3(:,1);
            otherwise, error('Modo no soportado (analítica limitada a modos 1-3).');
        end
        
    else
        error('Variable tipo_caso mal configurada. Revisa las cabeceras.');
    end
    
    % Alinear a columna explícitamente por si acaso
    if size(phi_exact, 2) > 1
        phi_exact = phi_exact';
    end
end
