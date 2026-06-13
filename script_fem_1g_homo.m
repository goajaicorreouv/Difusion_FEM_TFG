% DESCRIPCIÓN:
% Script principal para resolver difusión de 1 GRUPO mediante FEM (polinomios orden 3).
% Simula un reactor HOMOGÉNEO de 350 cm (todo el dominio se define como material 2).
% Calcula los 3 primeros k_eff, grafica modos y guarda resultados.
clearvars;
close all; 

% Forzar LaTeX globalmente en todas las gráficas
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');

L = 350;

% Propiedades físicas (definidas para 2 materiales posibles)
D = [1.446, 0.776];
sigma_a = [0.0077, 0.0244];
nu_sigma_f = [0, 0.0260];

% Configuración de Malla y Materiales
tamano_celdas=[25, 25 + zeros(1, 12), 25];
materiales = [2, 2 + zeros(1, 12), 2]; % Se asigna Material 2 a TODOS los elementos (Homogéneo)
N = 14;
grado_l = 3;

output_file = 'FEM_1G_HOM_N14_p4.mat';

% Inicialización de objetos (Clases de 1 Grupo)
malla = Malla1D(L, N, grado_l, tamano_celdas);
materiales = Materiales1D1g(materiales, D, sigma_a, nu_sigma_f, malla);
elemento = ElementoFinito(grado_l);
problema = ProblemaDifusion1D1g(malla, materiales, elemento);

% Resolución del sistema
problema = problema.ensamblar_matrices_fem();
problema = problema.aplicar_cc();

format long
problema = problema.resolver_autovalor(3);

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
title('\textbf{Espectro Espacial (1 Grupo Homogéneo)}', 'Interpreter', 'latex', 'FontSize', 16);
legend('Location', 'best', 'Interpreter', 'latex', 'FontSize', 14);
hold off;

%% --- 2. Gráfica de Validación (FEM vs Analítico) ---
figure('Name', 'Validación FEM vs Analítico', 'Position', [150, 150, 800, 750]);

x_anal = linspace(0, malla.L, 500)';
% Factor de traslación si el origin x=0 no es el borde (-L/2 a L/2). Aquí asumimos 0 a L.
phi_a1 = sin(1 * pi * x_anal / L);
phi_a2 = sin(2 * pi * x_anal / L);
phi_a3 = sin(3 * pi * x_anal / L);

analiticos = {phi_a1, phi_a2, phi_a3};
colores_anal = {'b', 'g', 'm'};

for i = 1:3 
    subplot(3, 1, i);
    hold on; grid on; box on;
    
    % En reactores homogéneos exactos, las raíces senoidales encajan perfecto con Dirichlet=0
    plot(x_anal, analiticos{i}, '-', 'Color', colores_anal{i}, 'LineWidth', 2, 'DisplayName', sprintf('Analítico Exacto ($n=%d$)', i));
    
    if isprop(problema, 'phi_nodos') 
        phi_fem = problema.phi_nodos(:, i);
    else 
        phi_fem = problema.phi_inc(1 : length(malla.x_nodos), i);
    end 
    
    phi_fem_norm = phi_fem / max(abs(phi_fem));
    
    idx_eval = max(1, floor(length(phi_fem_norm) / 5));
    if sign(phi_fem_norm(idx_eval)) ~= sign(analiticos{i}(idx_eval)) 
        phi_fem_norm = -phi_fem_norm;
    end
    
    plot(malla.x_nodos, phi_fem_norm, 'ro', 'MarkerSize', 5, 'MarkerFaceColor', 'r', 'DisplayName', sprintf('Nodos FEM ($p=%d$)', grado_l));
    
    ylabel('Flujo $\phi(x)$', 'Interpreter', 'latex', 'FontSize', 13);
    title(sprintf('\\textbf{Validación Modo %d}', i), 'Interpreter', 'latex', 'FontSize', 14);
    legend('Location', 'northeast', 'Interpreter', 'latex');
    hold off;
end 
xlabel('Posición $x$ (cm)', 'Interpreter', 'latex', 'FontSize', 14);

%% --- 3. Gráfica de Convergencia p - refinamiento (Los 3 Modos) ---
p_archivos = dir('FEM_1G_HOM_N14_p*.mat');
if ~isempty(p_archivos)
    % Ajustamos la ventana para que sea alta (vertical)
    figure('Name', 'p-Refinamiento por Modos', 'Position', [100, 50, 800, 1000]);
    
    % Paleta de colores sin negro: azul, verde, rojo, cian, magenta
    estilos = {'b:o', 'g--x', 'r-.s', 'c-d', 'm-*'};
    
    for modo_i = 1:3 
        % CAMBIO: Disposición en 3 filas y 1 columna
        subplot(3, 1, modo_i);
        hold on; grid on; box on;
        
        % Referencia analítica en gris claro
        plot(x_anal, analiticos{modo_i}, '-', 'Color', [0.7 0.7 0.7], 'LineWidth', 4, 'DisplayName', 'Analítico');
        
        for k = 1:length(p_archivos)
            f_data = load(p_archivos(k).name, 'malla', 'problema');
            p_str = regexp(p_archivos(k).name, 'p(\d+)', 'tokens');
            if ~isempty(p_str), p_actual = str2double(p_str{1}{1}); else, p_actual = k; end
            
            % Lógica de selección de flujo
            if isprop(f_data.problema, 'phi_nodos') && size(f_data.problema.phi_nodos, 2) >= modo_i 
                phi_fem = f_data.problema.phi_nodos(:, modo_i);
            else 
                phi_fem = f_data.problema.phi_inc(1 : length(f_data.malla.x_nodos), modo_i);
            end 
            
            % Normalización y corrección de fase
            phi_fem = phi_fem / max(abs(phi_fem));
            idx_ref = max(1, floor(length(phi_fem) / 5)); 
            val_analitico_ref = interp1(x_anal, analiticos{modo_i}, f_data.malla.x_nodos(idx_ref));
            if sign(phi_fem(idx_ref)) ~= sign(val_analitico_ref) 
                phi_fem = -phi_fem;
            end
            
            % Ploteo con el estilo correspondiente (sin negro)
            marca_idx = max(1, floor(length(f_data.malla.x_nodos) / 20));
            plot(f_data.malla.x_nodos, phi_fem, estilos{mod(k-1,5)+1}, 'LineWidth', 1.5, ... 
                 'MarkerIndices', 1 : marca_idx : length(f_data.malla.x_nodos), ... 
                 'DisplayName', sprintf('FEM $p=%d$', p_actual));
        end
        
        % Configuración de ejes y etiquetas
        xlim([0, f_data.malla.x_nodos(end)]);
        ylim([-1.1, 1.1]);
        
        % Si usas xline, cambiamos 'k' (negro) por gris oscuro [0.4 0.4 0.4]
        % xline(25, 'Color', [0.4 0.4 0.4], 'LineWidth', 0.8, 'HandleVisibility', 'off');
        
        xlabel('Posición $x$ (cm)', 'Interpreter', 'latex', 'FontSize', 12);
        ylabel(sprintf('Flujo $\\phi_%d(x)$', modo_i), 'Interpreter', 'latex', 'FontSize', 12);
        title(sprintf('\\textbf{Modo %d}', modo_i), 'Interpreter', 'latex');
        
        % CAMBIO: La leyenda ahora se muestra en las tres imágenes
        legend('Location', 'best', 'Interpreter', 'latex', 'FontSize', 10);
        
        hold off;
    end
end