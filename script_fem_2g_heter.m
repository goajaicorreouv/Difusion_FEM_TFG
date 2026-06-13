% =========================================================================
% DESCRIPCIÓN: 
% Script principal para resolver el problema de difusión de 2 GRUPOS 
% con FEM (orden 3). Configuración HETEROGÉNEA (Reflector - Núcleo - Reflector)
% - Material 2 (Extremos) : Reflector.
% - Material 1 (Centro)   : Núcleo.
% =========================================================================

clearvars;
close all;

% Forzar LaTeX globalmente en todas las gráficas para un acabado unificado
set(groot, 'defaultLegendInterpreter', 'latex');
set(groot, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'defaultTextInterpreter', 'latex');

L = 350;
N = 112;
grado_l = 4;

% Definición de la geometría y distribución de materiales 
tamano_celdas = [3.125, 3.125 + zeros(1, 110), 3.125];
materiales = [2,2,2,2,2,2,2,2, 1 + zeros(1, 96),2, 2,2, 2,2, 2,2, 2]; % Mat 1 (Combustible), Mat 2 (Reflector) 


% (Filas: Grupo 1, 2 | Columnas: Material 1, 2)
D         = [1.446,    2.0    % G1 
            0.40,    0.3];   % G2 
sigma_a    = [0.010,  0.00  
    0.085,  0.01];
sigma_s12    = [0.020,  0.04]; % Scattering 1->2
nu_sigma_f = [0.000,   0.0   %  G1
 0.135,   0.0];              % G2 

output_file = 'FEM_2G_HET_N112_p4.mat';

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
problema.graficar([1, 2, 3]);
save(output_file);
disp(problema.keff);
load(output_file);

%% --- 1.A. Gráfica del Espectro - Grupo 1 Rápidos (Modos 1, 2 y 3) ---
% Cargar explícitamente el archivo FEM de 2 GRUPOS
load(output_file, 'malla', 'problema');
figure('Name', 'Espectro de Modos - Rápidos', 'Position', [100, 100, 800, 500]);
hold on; grid on; box on;

colores = {'b', 'g', 'm'};
nombres = {'Fundamental ($n=1$)', '1er Armónico ($n=2$)', '2do Armónico ($n=3$)'};
n_pts = length(malla.x_nodos);

for i = 1:min(3, size(problema.phi_inc, 2)) 
    % Extraer el flujo rápido (los primeros N nodos corresponden siempre al Grupo 1) 
    if isprop(problema, 'phi_nodos') && size(problema.phi_nodos, 2) >= i 
        phi_plot = problema.phi_nodos(1:n_pts, i);
    else 
        phi_plot = problema.phi_inc(1:n_pts, i); % Fallback 
    end
    
    % Normalizar 
    phi_norm = phi_plot / max(abs(phi_plot));
    plot(malla.x_nodos, phi_norm, '-', 'Color', colores{i}, 'LineWidth', 2, 'DisplayName', nombres{i});
end

yline(0, 'k--', 'LineWidth', 1.2, 'HandleVisibility', 'off'); % Línea base de 0 
xlabel('Posición $x$ (cm)', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Flujo Rápido Normalizado $\phi_{1,n}(x)$', 'Interpreter', 'latex', 'FontSize', 14);
title('\textbf{Espectro Espacial: Grupo 1 (Neutrones Rápidos)}', 'Interpreter', 'latex', 'FontSize', 16);
legend('Location', 'best', 'Interpreter', 'latex', 'FontSize', 14);
hold off;

%% --- 1.B. Gráfica del Espectro - Grupo 2 Térmicos (Modos 1, 2 y 3) ---
figure('Name', 'Espectro de Modos - Térmicos', 'Position', [150, 150, 800, 500]);
hold on; grid on; box on;

for i = 1:min(3, size(problema.phi_inc, 2)) 
    % Extraer el flujo térmico (los segundos N nodos corresponden al Grupo 2) 
    if isprop(problema, 'phi_nodos') && size(problema.phi_nodos, 2) >= i 
        phi_plot = problema.phi_nodos(n_pts + 1 : 2 * n_pts, i);
    else 
        phi_plot = problema.phi_inc(n_pts + 1 : 2 * n_pts, i); % Fallback 
    end
    
    % Normalizar 
    phi_norm = phi_plot / max(abs(phi_plot));
    
    % Aplicar un estilo de línea distinto (opcional)
    plot(malla.x_nodos, phi_plot, '-', 'Color', colores{i}, 'LineWidth', 2, 'DisplayName', nombres{i});
end

yline(0, 'k--', 'LineWidth', 1.2, 'HandleVisibility', 'off'); % Línea base de 0 
xlabel('Posición $x$ (cm)', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Flujo Térmico Normalizado $\phi_{2,n}(x)$', 'Interpreter', 'latex', 'FontSize', 14);
title('\textbf{Espectro Espacial: Grupo 2 (Neutrones Térmicos)}', 'Interpreter', 'latex', 'FontSize', 16);
legend('Location', 'best', 'Interpreter', 'latex', 'FontSize', 14);
hold off;

%% --- 2. Gráfica de Validación (FEM vs Analítico) - GRUPO 1 (Rápido) ---
figure('Name', 'Validación FEM vs Analítico (Grupo 1)', 'Position', [150, 150, 800, 750]);

% 1. Calculamos la solución analítica de los 2 grupos
x_anal = linspace(0, malla.L, 500)';
[phi_a1, phi_a2, phi_a3] = calc_analitica_2g_het(x_anal);

% 2. Tomamos la banda de RÁPIDOS (columna 1) para plotear
phi_anal_1 = phi_a1(:,1);
phi_anal_2 = phi_a2(:,1);
phi_anal_3 = phi_a3(:,1);
analiticos = {phi_anal_1, phi_anal_2, phi_anal_3};

% 3. Normalizar SIGNO y MAGNITUD de la solución analítica
% Usamos x = L/4 como referencia: ahí todos los modos tienen signo definido
idx_quarter = floor(length(x_anal)/4);
for m_idx = 1:3
    % Ajustar signo
    if analiticos{m_idx}(idx_quarter) < 0
        analiticos{m_idx} = -analiticos{m_idx};
    end
    % [FALLO CORREGIDO]: Normalizar a amplitud máxima absoluta = 1
    analiticos{m_idx} = analiticos{m_idx} / max(abs(analiticos{m_idx}));
end

colores_anal = {'b', 'g', 'm'};
M = length(malla.x_nodos); % Número de nodos espaciales (para extraer bien el G1)

for i = 1:3 
    subplot(3, 1, i);
    hold on; grid on; box on;
    
    % Línea continua gruesa para la solución exacta
    plot(x_anal, analiticos{i}, '-', 'Color', colores_anal{i}, 'LineWidth', 2, ...
        'DisplayName', sprintf('Analítico R. ($n=%d$)', i));
    
    % Extraemos matriz FEM nodal (Aseguramos coger solo el G1: de 1 a M)
    if isprop(problema, 'phi_nodos') 
        phi_fem = problema.phi_nodos(1:M, i);
    else 
        phi_fem = problema.phi_inc(1:M, i);
    end 
    
    % Normalizar magnitud FEM a 1
    phi_fem_norm = phi_fem / max(abs(phi_fem));
    
    % Compensar signo invertido del FEM solver (muy típico) 
    idx_eval = max(1, floor(length(phi_fem_norm) / 5));
    val_anal_ref = interp1(x_anal, analiticos{i}, malla.x_nodos(idx_eval));
    if sign(phi_fem_norm(idx_eval)) ~= sign(val_anal_ref) 
        phi_fem_norm = -phi_fem_norm;
    end
    
    % Puntos rojos para los nodos discretos FEM
    plot(malla.x_nodos, phi_fem_norm, 'ro', 'MarkerSize', 5, 'MarkerFaceColor', 'r', ...
        'DisplayName', 'Nodos FEM ($p=3$)');
    
    ylabel('Flujo Rápido $\phi_1(x)$', 'Interpreter', 'latex', 'FontSize', 13);
    title(sprintf('\\textbf{Validación Modo %d - Grupo 1 (Rápido)}', i), 'Interpreter', 'latex', 'FontSize', 14);
    legend('Location', 'northeast', 'Interpreter', 'latex');
    hold off;
end 
xlabel('Posición $x$ (cm)', 'Interpreter', 'latex', 'FontSize', 14);

%% --- 2. Gráfica de Validación (FEM vs Analítico) - GRUPO 2 (Térmico) ---
figure('Name', 'Validación FEM vs Analítico (Grupo 2)', 'Position', [150, 150, 800, 750]);

% 1. Calculamos la solución analítica de los 2 grupos
x_anal = linspace(0, malla.L, 500)';
[phi_a1, phi_a2, phi_a3] = calc_analitica_2g_het(x_anal);

% 2. Tomamos la banda de TÉRMICOS (columna 2) para plotear
phi_anal_1 = phi_a1(:, 2);
phi_anal_2 = phi_a2(:, 2);
phi_anal_3 = phi_a3(:, 2);
analiticos = {phi_anal_1, phi_anal_2, phi_anal_3};

% 3. Normalizar SIGNO y MAGNITUD de la solución analítica
% Usamos x = L/4 como referencia: ahí todos los modos tienen signo definido
idx_quarter = floor(length(x_anal)/4);
for m_idx = 1:3
    % Ajustar signo
    if analiticos{m_idx}(idx_quarter) < 0
        analiticos{m_idx} = -analiticos{m_idx};
    end
    % [FALLO CORREGIDO]: Normalizar a amplitud 1 para igualar a la gráfica FEM
    analiticos{m_idx} = analiticos{m_idx} / max(abs(analiticos{m_idx}));
end

colores_anal = {'b', 'g', 'm'};
M = length(malla.x_nodos); % Número de grados de libertad espaciales (por grupo)

for i = 1:3 
    subplot(3, 1, i);
    hold on; grid on; box on;
    
    % Línea continua gruesa para la solución exacta del Grupo 2
    plot(x_anal, analiticos{i}, '-', 'Color', colores_anal{i}, 'LineWidth', 2, ...
        'DisplayName', sprintf('Analítico T. ($n=%d$)', i));
    
    % Extraemos matriz FEM nodal para el GRUPO 2 (índices de M+1 hasta 2*M)
    if isprop(problema, 'phi_nodos') 
        phi_fem = problema.phi_nodos(M + 1 : 2*M, i);
    else 
        phi_fem = problema.phi_inc(M + 1 : 2*M, i);
    end 
    
    % Normalizar magnitud FEM
    phi_fem_norm = phi_fem / max(abs(phi_fem));
    
    % Compensar signo invertido del FEM solver respecto al analítico
    idx_eval = max(1, floor(length(phi_fem_norm) / 5));
    val_anal_ref = interp1(x_anal, analiticos{i}, malla.x_nodos(idx_eval));
    if sign(phi_fem_norm(idx_eval)) ~= sign(val_anal_ref) 
        phi_fem_norm = -phi_fem_norm;
    end
    
    % Puntos rojos para los nodos discretos FEM
    plot(malla.x_nodos, phi_fem_norm, 'ro', 'MarkerSize', 5, 'MarkerFaceColor', 'r', ...
        'DisplayName', 'Nodos FEM ($p=3$)');
    
    ylabel('Flujo Térmico $\phi_2(x)$', 'Interpreter', 'latex', 'FontSize', 13);
    title(sprintf('\\textbf{Validación Modo %d - Grupo 2 (Térmico)}', i), 'Interpreter', 'latex', 'FontSize', 14);
    legend('Location', 'northeast', 'Interpreter', 'latex');
    hold off;
end 
xlabel('Posición $x$ (cm)', 'Interpreter', 'latex', 'FontSize', 14);

%% --- 3. Gráfica de Convergencia p - GRUPO 1 (Rápido) ---
figure('Name', 'p-Refinamiento - Grupo 1 (Rápido)', 'Position', [150, 100, 600, 950]);

% 1. Calcular y aislar la solución analítica del GRUPO 1 (Columna 1)
x_anal = linspace(0, malla.L, 500)';
[phi_a1, phi_a2, phi_a3] = calc_analitica_2g_het(x_anal);
analiticos_G1 = {phi_a1(:, 1), phi_a2(:, 1), phi_a3(:, 1)};

% Normalizar signo y magnitud del analítico del Grupo 1
idx_quarter = floor(length(x_anal)/4);
for m_idx = 1:3
    if analiticos_G1{m_idx}(idx_quarter) < 0
        analiticos_G1{m_idx} = -analiticos_G1{m_idx};
    end
    analiticos_G1{m_idx} = analiticos_G1{m_idx} / max(abs(analiticos_G1{m_idx}));
end

% 2. Cargar los datos de los 4 casos numéricos FEM
data_p1 = load("FEM_2G_HET_N14_p2.mat", 'malla', 'problema');
data_p2 = load("FEM_2G_HET_N14_p3.mat", 'malla', 'problema');
data_p3 = load("FEM_2G_HET_N14_p4.mat", 'malla', 'problema');
data_p4 = load("FEM_2G_HET_N14_p5.mat", 'malla', 'problema');

mallas_p = {data_p1.malla, data_p2.malla, data_p3.malla, data_p4.malla};
probs_p  = {data_p1.problema, data_p2.problema, data_p3.problema, data_p4.problema};

% Añadido un cuarto estilo (magenta con diamantes) y nombre
estilos   = {'b:o', 'g--x', 'r-.s', 'm-d'};
nombres_p = {'FEM $p=1$', 'FEM $p=2$', 'FEM $p=3$', 'FEM $p=4$'};

for modo_i = 1:3 
    subplot(3, 1, modo_i);
    hold on; grid on; box on;
    
    % Graficar Analítico del Grupo 1
    plot(x_anal, analiticos_G1{modo_i}, 'k-', 'LineWidth', 2, 'DisplayName', 'Analítico R.');
    
    % Bucle ampliado a 4 iteraciones
    for p_idx = 1:4 
        m_local = mallas_p{p_idx};
        prob_local = probs_p{p_idx};
        M_local = length(m_local.x_nodos); % Número de nodos de esta malla
        
        % Extraer estrictamente el Grupo 1 (filas 1 a M_local)
        if isprop(prob_local, 'phi_nodos')
            phi_fem = prob_local.phi_nodos(1:M_local, modo_i);
        else 
            phi_fem = prob_local.phi_inc(1:M_local, modo_i);
        end 
        
        % Normalizar numérico
        phi_fem = phi_fem / max(abs(phi_fem));
        
        % Alinear signos respecto al analítico
        idx_ref = max(1, floor(length(phi_fem) / 5)); 
        val_analitico_ref = interp1(x_anal, analiticos_G1{modo_i}, m_local.x_nodos(idx_ref));
        if sign(phi_fem(idx_ref)) ~= sign(val_analitico_ref) 
            phi_fem = -phi_fem;
        end
        
        plot(m_local.x_nodos, phi_fem, estilos{p_idx}, 'LineWidth', 1.5, ... 
             'MarkerIndices', 1 : max(1, floor(length(m_local.x_nodos) / 20)) : length(m_local.x_nodos), ... 
             'DisplayName', nombres_p{p_idx});
    end
    
    xlim([0, m_local.x_nodos(end)]); ylim([-1.1, 1.1]);
    xline(25, 'k-', 'LineWidth', 0.8, 'HandleVisibility', 'off');
    xline(m_local.x_nodos(end) - 25, 'k-', 'LineWidth', 0.8, 'HandleVisibility', 'off');
    
    xlabel('Posición $x$ (cm)', 'Interpreter', 'latex', 'FontSize', 12);
    ylabel(sprintf('Flujo Rápido $\\phi_{1,%d}(x)$', modo_i), 'Interpreter', 'latex', 'FontSize', 12);
    title(sprintf('\\textbf{Modo %d - Grupo 1 (Rápido)}', modo_i), 'Interpreter', 'latex');
    
    if modo_i == 1, legend('Location', 'best', 'Interpreter', 'latex', 'FontSize', 11); end
    hold off;
end

%% --- 4. Gráfica de Convergencia p - GRUPO 2 (Térmico) ---
figure('Name', 'p-Refinamiento - Grupo 2 (Térmico)', 'Position', [750, 100, 600, 950]);

% 1. Calcular y aislar la solución analítica del GRUPO 2 (Columna 2)
x_anal = linspace(0, malla.L, 500)';
[phi_a1, phi_a2, phi_a3] = calc_analitica_2g_het(x_anal);
analiticos_G2 = {phi_a1(:, 2), phi_a2(:, 2), phi_a3(:, 2)};

% Normalizar signo y magnitud del analítico del Grupo 2
idx_quarter = floor(length(x_anal)/4);
for m_idx = 1:3
    if analiticos_G2{m_idx}(idx_quarter) < 0
        analiticos_G2{m_idx} = -analiticos_G2{m_idx};
    end
    analiticos_G2{m_idx} = analiticos_G2{m_idx} / max(abs(analiticos_G2{m_idx}));
end

% 2. Cargar los datos de los 4 casos numéricos FEM
data_p1 = load("FEM_2G_HET_N14_p2.mat", 'malla', 'problema');
data_p2 = load("FEM_2G_HET_N14_p3.mat", 'malla', 'problema');
data_p3 = load("FEM_2G_HET_N14_p4.mat", 'malla', 'problema');
data_p4 = load("FEM_2G_HET_N14_p5.mat", 'malla', 'problema');

mallas_p = {data_p1.malla, data_p2.malla, data_p3.malla, data_p4.malla};
probs_p  = {data_p1.problema, data_p2.problema, data_p3.problema, data_p4.problema};

% Añadido un cuarto estilo (magenta con diamantes) y nombre
estilos   = {'b:o', 'g--x', 'r-.s', 'm-d'};
nombres_p = {'FEM $p=1$', 'FEM $p=2$', 'FEM $p=3$', 'FEM $p=4$'};

for modo_i = 1:3 
    subplot(3, 1, modo_i);
    hold on; grid on; box on;
    
    % Graficar Analítico del Grupo 2
    plot(x_anal, analiticos_G2{modo_i}, 'k-', 'LineWidth', 2, 'DisplayName', 'Analítico T.');
    
    % Bucle ampliado a 4 iteraciones
    for p_idx = 1:4 
        m_local = mallas_p{p_idx};
        prob_local = probs_p{p_idx};
        M_local = length(m_local.x_nodos); % Nodos espaciales de esta malla
        
        % Extraer estrictamente el Grupo 2 (filas M_local+1 a 2*M_local)
        if isprop(prob_local, 'phi_nodos')
            phi_fem = prob_local.phi_nodos(M_local + 1 : 2*M_local, modo_i);
        else 
            phi_fem = prob_local.phi_inc(M_local + 1 : 2*M_local, modo_i);
        end 
        
        % Normalizar numérico
        phi_fem = phi_fem / max(abs(phi_fem));
        
        % Alinear signos respecto al analítico
        idx_ref = max(1, floor(length(phi_fem) / 5)); 
        val_analitico_ref = interp1(x_anal, analiticos_G2{modo_i}, m_local.x_nodos(idx_ref));
        if sign(phi_fem(idx_ref)) ~= sign(val_analitico_ref) 
            phi_fem = -phi_fem;
        end
        
        plot(m_local.x_nodos, phi_fem, estilos{p_idx}, 'LineWidth', 1.5, ... 
             'MarkerIndices', 1 : max(1, floor(length(m_local.x_nodos) / 20)) : length(m_local.x_nodos), ... 
             'DisplayName', nombres_p{p_idx});
    end
    
    xlim([0, m_local.x_nodos(end)]); ylim([-1.1, 1.1]);
    xline(25, 'k-', 'LineWidth', 0.8, 'HandleVisibility', 'off');
    xline(m_local.x_nodos(end) - 25, 'k-', 'LineWidth', 0.8, 'HandleVisibility', 'off');
    
    xlabel('Posición $x$ (cm)', 'Interpreter', 'latex', 'FontSize', 12);
    ylabel(sprintf('Flujo Térmico $\\phi_{2,%d}(x)$', modo_i), 'Interpreter', 'latex', 'FontSize', 12);
    title(sprintf('\\textbf{Modo %d - Grupo 2 (Térmico)}', modo_i), 'Interpreter', 'latex');
    
    if modo_i == 1, legend('Location', 'best', 'Interpreter', 'latex', 'FontSize', 11); end
    hold off;
end