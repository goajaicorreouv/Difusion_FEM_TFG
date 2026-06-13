% DESCRIPCIÓN:
% Script principal para resolver difusión de 1 GRUPO en configuración HETEROGÉNEA.
% Utiliza FEM con polinomios de grado 4.
% Configuración "Reflector-Núcleo-Reflector":
% - Extremos (Celdas 1 y 14): Material 1 (Reflector).
% - Centro (Celdas 2-13): Material 2 (Combustible).
clearvars;
close all; 

% Forzar LaTeX globalmente en todas las gráficas
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');
%%
L = 350;

% Propiedades Físicas (Mat 1: Reflector, Mat 2: Combustible)
D = [1.446, 0.776];
sigma_a = [0.0077, 0.0244];
nu_sigma_f = [0, 0.0260]; % Mat 1 no tiene fisión

% Definición de la Malla y distribución de materiales
tamano_celdas=[3.125, 3.125 + zeros(1, 110), 3.125];
materiales = [1,1,1,1,1,1,1,1, 2 + zeros(1, 96),1,1,1, 1,1,1,1, 1]; % Mat 1 en extremos, Mat 2 en el centro
N = 112;
grado_l = 5;

output_file = 'FEM_1G_HET_N112_p5.mat';

% Inicialización de objetos
malla = Malla1D(L, N, grado_l, tamano_celdas);
materiales = Materiales1D1g(materiales, D, sigma_a, nu_sigma_f, malla);
elemento = ElementoFinito(grado_l);
problema = ProblemaDifusion1D1g(malla, materiales, elemento);

% Resolución del sistema
problema = problema.ensamblar_matrices_fem();
problema = problema.aplicar_cc();

format long
problema = problema.resolver_autovalor(3); % Calcular 3 primeros modos

% Resultados
problema.graficar([1,2,3]);
save(output_file);
disp(problema.keff);
load(output_file);

%% --- 1. Gráfica del Espectro (Modos 1, 2 y 3) ---
figure('Name', 'Espectro de Modos', 'Position', [100, 100, 800, 500]);
hold on; grid on; box on;

colores = {'b', 'g', 'm'};
nombres = {'Fundamental ($n=1$)', '1er Armónico ($n=2$)', '2do Armónico ($n=3$)'};
n_pts = length(malla.x_nodos);

for i = 1:min(3, size(problema.phi_inc, 2)) 
    if isprop(problema, 'phi_nodos') && size(problema.phi_nodos, 2) >= i 
        phi_plot = problema.phi_nodos(1:n_pts, i);
    else 
        phi_plot = problema.phi_inc(1:n_pts, i); 
    end
    
    phi_norm = phi_plot / max(abs(phi_plot));
    plot(malla.x_nodos, phi_norm, '-', 'Color', colores{i}, 'LineWidth', 2, 'DisplayName', nombres{i});
end

yline(0, 'k--', 'LineWidth', 1.2, 'HandleVisibility', 'off'); 
xlabel('Posición $x$ (cm)', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Flujo Normalizado $\phi_n(x)$', 'Interpreter', 'latex', 'FontSize', 14);
title('\textbf{Espectro Espacial}', 'Interpreter', 'latex', 'FontSize', 16);
legend('Location', 'best', 'Interpreter', 'latex', 'FontSize', 14);
hold off;

%% --- 2. Gráfica de Validación (FEM vs Analítico) ---
figure('Name', 'Validación FEM vs Analítico', 'Position', [150, 150, 800, 750]);

x_anal = linspace(0, malla.L, 500)';
[phi_a1, phi_a2, phi_a3] = calc_analitica_1g_het(x_anal);

analiticos = {phi_a1, phi_a2, phi_a3};
colores_anal = {'b', 'g', 'm'};

for i = 1:3 
    subplot(3, 1, i);
    hold on; grid on; box on;
    
    h_analitico_2 = plot(x_anal, analiticos{i}, '-', 'Color', colores_anal{i}, 'LineWidth', 2, 'DisplayName', sprintf('Analítico ($n=%d$)', i));
    
    if isprop(problema, 'phi_nodos') 
        phi_fem = problema.phi_nodos(:, i);
    else 
        phi_fem = problema.phi_inc(1 : length(malla.x_nodos), i);
    end 
    
    phi_fem_norm = phi_fem / max(abs(phi_fem));
    
    idx_eval = max(1, floor(length(phi_fem_norm) / 5));
    if sign(phi_fem_norm(idx_eval)) ~= sign(analiticos{i}(idx_eval)) 
        analiticos{i} = -analiticos{i};
        set(h_analitico_2, 'YData', analiticos{i});
    end
    
    plot(malla.x_nodos, phi_fem_norm, 'ro', 'MarkerSize', 5, 'MarkerFaceColor', 'r', 'DisplayName', sprintf('Nodos FEM ($p=%d$)', grado_l));
    
    ylabel('Flujo $\phi(x)$', 'Interpreter', 'latex', 'FontSize', 13);
    title(sprintf('\\textbf{Validación Modo %d}', i), 'Interpreter', 'latex', 'FontSize', 14);
    legend('Location', 'northeast', 'Interpreter', 'latex');
    hold off;
end 
xlabel('Posición $x$ (cm)', 'Interpreter', 'latex', 'FontSize', 14);

%% --- 3. Gráfica de Convergencia p - refinamiento (Los 3 Modos) ---
% Buscar dinámicamente qué grados tenemos ejecutados
p_archivos = dir('FEM_1G_HET_N14_p*.mat');
if ~isempty(p_archivos)
    % Ajustada la altura y anchura
    figure('Name', 'p-Refinamiento por Modos', 'Position', [100, 100, 800, 1000]);
    
    % Estilos sin negro
    estilos = {'b:o', 'g--x', 'r-.s', 'c-d', 'm-*'};
    
    for modo_i = 1:3 
        % 3 filas, 1 columna
        subplot(3, 1, modo_i);
        hold on; grid on; box on;
        
        % Curva analítica (se dibuja una vez y NO se modifica)
        plot(x_anal, analiticos{modo_i}, '-', 'Color', [0.7 0.7 0.7], 'LineWidth', 4, 'DisplayName', 'Analítico');
        
        for k = 1:length(p_archivos)
            f_data = load(p_archivos(k).name, 'malla', 'problema');
            p_str = regexp(p_archivos(k).name, 'p(\d+)', 'tokens');
            if ~isempty(p_str), p_actual = str2double(p_str{1}{1}); else, p_actual = k; end
            
            if isprop(f_data.problema, 'phi_nodos') && size(f_data.problema.phi_nodos, 2) >= modo_i 
                phi_fem = f_data.problema.phi_nodos(:, modo_i);
            else 
                phi_fem = f_data.problema.phi_inc(1 : length(f_data.malla.x_nodos), modo_i);
            end 
            
            % Normalización
            phi_fem = phi_fem / max(abs(phi_fem));
            
            % --- CORRECCIÓN DE FASE (IGUAL QUE EN EL CASO HOMOGÉNEO) ---
            idx_ref = max(2, floor(length(phi_fem) / 5)); 
            val_analitico_ref = interp1(x_anal, analiticos{modo_i}, f_data.malla.x_nodos(idx_ref));
            if sign(phi_fem(idx_ref)) ~= sign(val_analitico_ref) 
                phi_fem = -phi_fem; % Invertimos el FEM, no el analítico
            end
            % -----------------------------------------------------------
            
            marca_idx = max(1, floor(length(f_data.malla.x_nodos) / 20));
            plot(f_data.malla.x_nodos, phi_fem, estilos{mod(k-1,5)+1}, 'LineWidth', 1.5, ... 
                 'MarkerIndices', 1 : marca_idx : length(f_data.malla.x_nodos), ... 
                 'DisplayName', sprintf('FEM $p=%d$', p_actual));
        end
        
        xlim([0, f_data.malla.x_nodos(end)]);
        ylim([-1.1, 1.1]);
        
        % Líneas verticales en gris oscuro
        xline(25, '-', 'Color', [0.4 0.4 0.4], 'LineWidth', 0.8, 'HandleVisibility', 'off');
        xline(f_data.malla.x_nodos(end) - 25, '-', 'Color', [0.4 0.4 0.4], 'LineWidth', 0.8, 'HandleVisibility', 'off');
        
        xlabel('Posición $x$ (cm)', 'Interpreter', 'latex', 'FontSize', 12);
        ylabel(sprintf('Flujo $\\phi_%d(x)$', modo_i), 'Interpreter', 'latex', 'FontSize', 12);
        title(sprintf('\\textbf{Modo %d}', modo_i), 'Interpreter', 'latex');
        
        legend('Location', 'best', 'Interpreter', 'latex', 'FontSize', 11);
        
        hold off;
    end
end