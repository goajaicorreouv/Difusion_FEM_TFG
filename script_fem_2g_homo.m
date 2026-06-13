% DESCRIPCIÓN:
% Script principal para resolver el problema de difusión de 2 GRUPOS con FEM (orden 3).
% Configuración HOMOGÉNEA:
% - Todo el dominio (elementos 1 a N) se define como Material 1 (Núcleo).

clearvars;
close all; 

% Forzar LaTeX globalmente en todas las gráficas para un acabado unificado
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');

L = 350;
N = 14;
grado_l = 4;

% Definición de la geometría y distribución de materiales
tamano_celdas = [25, 25 + zeros(1, 12), 25];

%Todo el reactor es Material 1
materiales = ones(1, N); 

% Propiedades Físicas (dejamos las del caso heterogéno)
% (Filas: Grupo 1, 2 | Columnas: Material 1, 2)
D         = [1.446,    2.0    % G1 
            0.40,    0.3];   % G2 
sigma_a    = [0.010,  0.00  
    0.085,  0.01];
sigma_s12    = [0.020,  0.04]; % Scattering 1->2
nu_sigma_f = [0.000,   0.0   %  G1
 0.135,   0.0];              % G2 

output_file = 'FEM_2G_HOM_N14_p4.mat';

% Inicialización de objetos
malla = Malla1D(L, N, grado_l, tamano_celdas);
materiales = Materiales1D2g(materiales, D, sigma_a, nu_sigma_f, sigma_s12, malla);
elemento = ElementoFinito(grado_l);
problema = ProblemaDifusion1D2g(malla, materiales, elemento);

% Resolución
problema = problema.ensamblar_matrices_fem();
problema = problema.aplicar_cc();
problema = problema.resolver_autovalor(3);

% Resultados
problema.graficar([1,2,3]);
save(output_file);
disp(problema.keff);
load(output_file);

%% --- 1.A. Gráfica del Espectro - Grupo 1 Rápidos (Modos 1, 2 y 3) ---
figure('Name', 'Espectro de Modos - Rápidos', 'Position', [100, 100, 800, 500]);
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
    plot(malla.x_nodos, phi_plot, '-', 'Color', colores{i}, 'LineWidth', 2, 'DisplayName', nombres{i});
end

yline(0, 'k--', 'LineWidth', 1.2, 'HandleVisibility', 'off'); 
xlabel('Posición $x$ (cm)', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Flujo Rápido Normalizado $\phi_{1,n}(x)$', 'Interpreter', 'latex', 'FontSize', 14);
title('\textbf{Espectro Espacial: Grupo 1 (Neutrones Rápidos)}', 'Interpreter', 'latex', 'FontSize', 16);
legend('Location', 'best', 'Interpreter', 'latex', 'FontSize', 14);
hold off;

%% --- 1.B. Gráfica del Espectro - Grupo 2 Térmicos (Modos 1, 2 y 3) ---
figure('Name', 'Espectro de Modos - Térmicos', 'Position', [150, 150, 800, 500]);
hold on; grid on; box on;

for i = 1:min(3, size(problema.phi_inc, 2)) 
    if isprop(problema, 'phi_nodos') && size(problema.phi_nodos, 2) >= i 
        phi_plot = problema.phi_nodos(n_pts + 1 : 2 * n_pts, i);
    else 
        phi_plot = problema.phi_inc(n_pts + 1 : 2 * n_pts, i); 
    end
    
    phi_norm = phi_plot / max(abs(phi_plot));
    plot(malla.x_nodos, phi_plot, '-', 'Color', colores{i}, 'LineWidth', 2, 'DisplayName', nombres{i});
end

yline(0, 'k--', 'LineWidth', 1.2, 'HandleVisibility', 'off'); 
xlabel('Posición $x$ (cm)', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Flujo Térmico Normalizado $\phi_{2,n}(x)$', 'Interpreter', 'latex', 'FontSize', 14);
title('\textbf{Espectro Espacial: Grupo 2 (Neutrones Térmicos)}', 'Interpreter', 'latex', 'FontSize', 16);
legend('Location', 'best', 'Interpreter', 'latex', 'FontSize', 14);
hold off;

%% --- 2. Gráfica de Validación (FEM vs Analítico) - GRUPO 1 (Rápido - Homogéneo) ---
figure('Name', 'Validación FEM vs Analítico Homogéneo (Grupo 1)', 'Position', [150, 150, 800, 750]);

% 1. Calculamos la solución analítica exacta homogénea de los 2 grupos
x_anal = linspace(0, malla.L, 500)';
[phi_a1, phi_a2, phi_a3] = calc_analitica_2g_hom(x_anal);

% 2. Tomamos la banda de RÁPIDOS (columna 1) para plotear
phi_anal_1 = phi_a1(:, 1);
phi_anal_2 = phi_a2(:, 1);
phi_anal_3 = phi_a3(:, 1);
analiticos = {phi_anal_1, phi_anal_2, phi_anal_3};

% 3. Normalizar SIGNO y MAGNITUD de la solución analítica a 1
idx_quarter = floor(length(x_anal)/4);
for m_idx = 1:3
    if analiticos{m_idx}(idx_quarter) < 0
        analiticos{m_idx} = -analiticos{m_idx};
    end
    analiticos{m_idx} = analiticos{m_idx} / max(abs(analiticos{m_idx}));
end

colores_anal = {'b', 'g', 'm'};
M = length(malla.x_nodos); % Número de nodos espaciales en la malla

for i = 1:3 
    subplot(3, 1, i);
    hold on; grid on; box on;
    
    % Línea continua gruesa para la solución exacta del Grupo 1
    plot(x_anal, analiticos{i}, '-', 'Color', colores_anal{i}, 'LineWidth', 2, ...
        'DisplayName', sprintf('Analítico R. ($n=%d$)', i));
    
    % Extraemos la solución numérica del GRUPO 1 (filas 1 a M)
    if isprop(problema, 'phi_nodos') 
        phi_fem = problema.phi_nodos(1:M, i);
    else 
        phi_fem = problema.phi_inc(1:M, i);
    end 
    
    % Normalizar magnitud del FEM a 1
    phi_fem_norm = phi_fem / max(abs(phi_fem));
    
    % Compensar signo invertido del FEM solver respecto al analítico
    idx_eval = max(1, floor(length(phi_fem_norm) / 5));
    val_anal_ref = interp1(x_anal, analiticos{i}, malla.x_nodos(idx_eval));
    if sign(phi_fem_norm(idx_eval)) ~= sign(val_anal_ref) 
        phi_fem_norm = -phi_fem_norm;
    end
    
    % Puntos rojos para los nodos discretos FEM
    plot(malla.x_nodos, phi_fem_norm, 'ro', 'MarkerSize', 5, 'MarkerFaceColor', 'r', ...
        'DisplayName', sprintf('Nodos FEM Rápido ($p=%d$)', malla.grado_l));
    
    ylabel('Flujo Normalizado $\phi_1(x)$', 'Interpreter', 'latex', 'FontSize', 13);
    title(sprintf('\\textbf{Validación Modo %d - Grupo 1 (Rápido)}', i), 'Interpreter', 'latex', 'FontSize', 14);
    legend('Location', 'northeast', 'Interpreter', 'latex');
    hold off;
end 
xlabel('Posición $x$ (cm)', 'Interpreter', 'latex', 'FontSize', 14);

%% --- 2. Gráfica de Validación (FEM vs Analítico) - GRUPO 2 (Térmico - Homogéneo) ---
figure('Name', 'Validación FEM vs Analítico Homogéneo (Grupo 2)', 'Position', [150, 150, 800, 750]);

% 1. Calculamos la solución analítica exacta homogénea
x_anal = linspace(0, malla.L, 500)';
[phi_a1, phi_a2, phi_a3] = calc_analitica_2g_hom(x_anal);

% Extraemos ambos grupos para preservar el ratio físico
phi_anal_1_modos = {phi_a1(:, 1), phi_a2(:, 1), phi_a3(:, 1)}; % Rápidos
phi_anal_2_modos = {phi_a1(:, 2), phi_a2(:, 2), phi_a3(:, 2)}; % Térmicos

analiticos = cell(1,3);
idx_quarter = floor(length(x_anal)/4);

for m_idx = 1:3
    % NORMALIZACIÓN CONJUNTA: Buscamos el máximo absoluto de todo el modo (G1 + G2)
    max_global = max(max(abs(phi_anal_1_modos{m_idx})), max(abs(phi_anal_2_modos{m_idx})));
    
    % Dividimos el térmico por ese máximo global para mantener el ratio
    analiticos{m_idx} = phi_anal_2_modos{m_idx} / max_global;
    
    % Alineamos el signo tomando como referencia el grupo rápido (igual que FEM)
    if phi_anal_1_modos{m_idx}(idx_quarter) < 0
        analiticos{m_idx} = -analiticos{m_idx};
    end
end

colores_anal = {'b', 'g', 'm'};
M = length(malla.x_nodos); % Número de nodos espaciales en la malla

for i = 1:3 
    subplot(3, 1, i);
    hold on; grid on; box on;
    
    % Línea continua gruesa para la solución exacta térmica
    plot(x_anal, analiticos{i}, '-', 'Color', colores_anal{i}, 'LineWidth', 2, ...
        'DisplayName', sprintf('Analítico T. ($n=%d$)', i));
    
    % Extraemos la solución numérica del GRUPO 2 (filas M+1 a 2*M)
    if isprop(problema, 'phi_nodos') 
        phi_fem = problema.phi_nodos(M+1 : 2*M, i);
    else 
        phi_fem = problema.phi_inc(M+1 : 2*M, i);
    end 
    
    % ¡CORRECCIÓN!: NO normalizamos el térmico numérico. 
    % Ya viene con la escala correcta del vector acoplado del solver.
    phi_fem_plot = phi_fem; 
    
    % Compensar signo invertido del FEM solver respecto al analítico si fuese necesario
    idx_eval = max(1, floor(length(phi_fem_plot) / 5));
    val_anal_ref = interp1(x_anal, analiticos{i}, malla.x_nodos(idx_eval));
    if sign(phi_fem_plot(idx_eval)) ~= sign(val_anal_ref) && val_anal_ref ~= 0
        phi_fem_plot = -phi_fem_plot;
    end
    
    % Puntos rojos para los nodos discretos FEM
    plot(malla.x_nodos, phi_fem_plot, 'ro', 'MarkerSize', 5, 'MarkerFaceColor', 'r', ...
        'DisplayName', sprintf('Nodos FEM Térmico ($p=%d$)', malla.grado_l));
    
    ylabel('Flujo Térmico $\phi_2(x)$', 'Interpreter', 'latex', 'FontSize', 13);
    title(sprintf('\\textbf{Validación Modo %d - Grupo 2 (Térmico)}', i), 'Interpreter', 'latex', 'FontSize', 14);
    legend('Location', 'northeast', 'Interpreter', 'latex');
    hold off;
end 
xlabel('Posición $x$ (cm)', 'Interpreter', 'latex', 'FontSize', 14);

%% --- 3. Gráfica de Convergencia p - refinamiento (Los 3 Modos) ---
p_archivos = dir('FEM_2G_HOM_N14_p*.mat');
if ~isempty(p_archivos)
    figure('Name', 'p-Refinamiento por Modos', 'Position', [100, 200, 1200, 450]);
    estilos = {'b:o', 'g--x', 'r-.s', 'k-d', 'm-*'};
    
    for modo_i = 1:3 
        subplot(1, 3, modo_i);
        hold on; grid on; box on;
        
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
            
            % Evaluamos el p-refinamiento usando la mitad correspondiente a Rápidos
            phi_fem = phi_fem(1:n_pts);
            phi_fem = phi_fem / max(abs(phi_fem));
            
            idx_ref = max(1, floor(length(phi_fem) / 5)); 
            val_analitico_ref = interp1(x_anal, analiticos{modo_i}, f_data.malla.x_nodos(idx_ref));
            if sign(phi_fem(idx_ref)) ~= sign(val_analitico_ref) 
                phi_fem = -phi_fem;
            end
            
            marca_idx = max(1, floor(length(f_data.malla.x_nodos) / 20));
            plot(f_data.malla.x_nodos, phi_fem, estilos{mod(k-1,5)+1}, 'LineWidth', 1.5, ... 
                 'MarkerIndices', 1 : marca_idx : length(f_data.malla.x_nodos), ... 
                 'DisplayName', sprintf('FEM $p=%d$', p_actual));
        end
        
        xlim([0, f_data.malla.x_nodos(end)]);
        ylim([-1.1, 1.1]);
        
        xlabel('Posición $x$ (cm)', 'Interpreter', 'latex', 'FontSize', 12);
        ylabel(sprintf('Flujo $\\phi_%d(x)$', modo_i), 'Interpreter', 'latex', 'FontSize', 12);
        title(sprintf('\\textbf{Modo %d}', modo_i), 'Interpreter', 'latex');
        
        if modo_i == 2 
            legend('Location', 'best', 'Interpreter', 'latex', 'FontSize', 11);
        end
        hold off;
    end
end

%% --- 4. Gráfica de Comparación de Archivos (FEM de distintos grados vs Analítico) ---
% Definimos la lista de archivos a comparar y las etiquetas correspondientes
archivos_comp = { 'FEM_2G_HOM_N14_p2.mat','FEM_2G_HOM_N14_p3.mat', ...
                 'FEM_2G_HOM_N14_p4.mat', 'FEM_2G_HOM_N14_p5.mat', ...
                 'FEM_2G_HOM_N14_p6.mat'};
             
% Asumimos que el primer archivo corresponde a p=1 u otro valor base. 
% Extraeremos los grados para la leyenda de forma visual.
etiquetas_p = {'p=2', 'p=3', 'p=4', 'p=5', 'p=6'};
marcadores = {'o', 's', '^', 'd', 'v'}; % Marcadores distintos para diferenciar

% Crear la figura de comparación
figure('Name', 'Comparación Grados FEM', 'Position', [200, 100, 900, 900]);

% Puntos para la solución analítica continua
x_anal = linspace(0, L, 1000)';

for modo = 1:3 
    subplot(3, 1, modo);
    hold on; grid on; box on;
    
    % 1. Calcular y dibujar la solución analítica para el modo actual
    phi_anal_exacto = sin(modo * pi * x_anal / L);
    plot(x_anal, phi_anal_exacto, 'k-', 'LineWidth', 2.5, ...
         'DisplayName', sprintf('Analítico Exacto ($n=%d$)', modo));
    
    % 2. Iterar sobre cada archivo y superponer su solución
    for i = 1:length(archivos_comp)
        if isfile(archivos_comp{i})
            datos_comp = load(archivos_comp{i});
            x_nodos_comp = datos_comp.malla.x_nodos;
            n_pts_comp = length(x_nodos_comp);
            
            % Extraer el flujo (usamos el Grupo 1 - Rápido para la comparativa)
            if isprop(datos_comp.problema, 'phi_inc') && size(datos_comp.problema.phi_inc, 2) >= modo
                phi_comp = datos_comp.problema.phi_inc(1:n_pts_comp, modo);
            else
                continue; % Si el archivo no tiene este modo calculado, lo saltamos
            end
            
            % Normalizar el flujo
            phi_comp_norm = phi_comp / max(abs(phi_comp));
            
            % Alinear el signo con la analítica comprobando en un punto significativo
            % Para el modo n, evitamos los nodos en los ceros de la función
            idx_eval = max(2, floor(n_pts_comp / (2 * modo))); 
            analitico_eval = sin(modo * pi * x_nodos_comp(idx_eval) / L);
            
            if sign(phi_comp_norm(idx_eval)) ~= sign(analitico_eval) && abs(analitico_eval) > 1e-3
                phi_comp_norm = -phi_comp_norm;
            end
            
            % Trazar la solución del archivo
            plot(x_nodos_comp, phi_comp, '--', 'Marker', marcadores{i}, ...
                 'MarkerSize', 5, 'LineWidth', 1.2, ...
                 'DisplayName', sprintf('FEM %s', etiquetas_p{i}));
        else
            warning('El archivo %s no se encuentra en el directorio.', archivos_comp{i});
        end
    end
    
    % Configurar etiquetas y leyenda para cada subplot
    ylabel(sprintf('Flujo Normalizado $\\phi_{1,%d}(x)$', modo), 'Interpreter', 'latex', 'FontSize', 12);
    title(sprintf('\\textbf{Comparación Convergencia - Modo %d}', modo), 'Interpreter', 'latex', 'FontSize', 14);
    
    % Configurar la leyenda en dos columnas para que no ocupe demasiado espacio vertical
    legend('Location', 'best', 'Interpreter', 'latex', 'FontSize', 10, 'NumColumns', 2);
    hold off;
end 

% Etiqueta global para el eje X en el último subplot
xlabel('Posición $x$ (cm)', 'Interpreter', 'latex', 'FontSize', 14);