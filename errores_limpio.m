%% ANÁLISIS DE ERRORES - SOLUCIONES ANALÍTICAS
tic;

%% ========================================================================
%  SOLUCIÓN 1G HOMOGÉNEO
%  ========================================================================
clear all
load('FEM_1G_HOM_N14_p3.mat')

% Solución analítica
m = length(problema.keff);
n_modos = 1:m;
k_eff_analitico_1g = nu_sigma_f(2) ./ (D(2) * (pi * n_modos / L).^2 + sigma_a(2));

phi_analitico_1g = @(x, n) sin(pi * n * x / L);

k_eff_numerico_1g = problema.keff(:);

ErA_1g = zeros(m, 1);
ErR_1g = zeros(m, 1);
x_nodos = malla.x_nodos;

for modo = 1:m
    phi_numerico_1g = problema.phi_inc(:, modo);
    phi_numerico_1g = phi_numerico_1g(:);

    phi_analitico_norm_1g = phi_analitico_1g(x_nodos, modo);
    phi_analitico_norm_1g = phi_analitico_norm_1g(:);

    if phi_analitico_norm_1g(2) < 0
        phi_analitico_norm_1g(:) = -phi_analitico_norm_1g(:);
    end
    phi_analitico_norm_1g(:) = phi_analitico_norm_1g(:) / max(phi_analitico_norm_1g(:));

    ErA_1g(modo) = abs(k_eff_analitico_1g(modo) - k_eff_numerico_1g(modo)) * 10^5;
    ErR_1g(modo) = rmse(phi_analitico_norm_1g, phi_numerico_1g);
end

figure(1);
hold on; grid on;
colores = {'b', 'g', 'm', 'c'};
n_puntos = length(problema.phi_inc(:, 1));
x_plot = linspace(0, L, n_puntos)';

for modo = 1:m
    phi_analitico_norm = phi_analitico_1g(x_plot, modo);
    phi_analitico_norm = phi_analitico_norm(:) / max(abs(phi_analitico_norm(:)));
    phi_numerico_norm = problema.phi_inc(:, modo) / max(abs(problema.phi_inc(:, modo)));

    plot(x_plot, phi_analitico_norm, 'Color', colores{modo}, 'LineStyle', '-', ...
         'LineWidth', 2, 'DisplayName', sprintf('Analítica Modo %d', modo));
    plot(x_plot, phi_numerico_norm, 'Color', colores{modo}, 'LineStyle', '--', ...
         'LineWidth', 1.5, 'Marker', 'o', 'MarkerSize', 4, 'DisplayName', sprintf('FEM Modo %d', modo));
end
xlabel('x (cm)');
ylabel('\phi(x)');
title('Solución 1G Homogéneo - Comparación Analítica vs FEM (Todos los modos)');
legend('Location', 'best');

fprintf('\n=== RESULTADOS 1G HOMOGÉNEO- POR MODO ===\n');
fprintf('%-10s | %-18s | %-18s | %-15s | %-15s\n', 'Modo', 'k_eff Analítico', 'k_eff Numérico', 'Error k_eff (×10^5)', 'RMSE Flujo');
fprintf('----------------------------------------------------------------------------------------\n');
for modo = 1:m
    fprintf('%-10d | %18.10f | %18.10f | %15.6f | %15.6e\n', modo, k_eff_analitico_1g(modo), k_eff_numerico_1g(modo), ErA_1g(modo), ErR_1g(modo));
end
fprintf('----------------------------------------------------------------------------------------\n');

fid = fopen('tabla_errores_1g_homogeneo.tex', 'w');
fprintf(fid, '\\begin{table}[h]\n');
fprintf(fid, '\\centering\n');
fprintf(fid, '\\caption{Resultados de errores para el caso 1G Homogéneo- Comparación analítica vs numérica}\n');
fprintf(fid, '\\label{tab:errores_1g}\n');
fprintf(fid, '\\begin{tabular}{c|cc|cc}\n');
fprintf(fid, '\\toprule\n');
fprintf(fid, '\\textbf{Modo} & \\textbf{$k_{eff}$ Analítico} & \\textbf{$k_{eff}$ Numérico} & \\textbf{Error $k_{eff}$ ($\\times 10^5$)} & \\textbf{RMSE Flujo} \\\\\n');
fprintf(fid, '\\midrule\n');
for modo = 1:m
    fprintf(fid, '%d & %.10f & %.10f & %.6f & %.6e \\\\\n', modo, k_eff_analitico_1g(modo), k_eff_numerico_1g(modo), ErA_1g(modo), ErR_1g(modo));
end
fprintf(fid, '\\bottomrule\n');
fprintf(fid, '\\end{tabular}\n');
fprintf(fid, '\\end{table}\n');
fclose(fid);
fprintf('\n Tabla LaTeX guardada en: tabla_errores_1g_homogeneo.tex\n');

m_hom = m;
ErA_1g_hom = ErA_1g;
ErR_1g_hom = ErR_1g;
k_eff_analitico_1g_hom = k_eff_analitico_1g;
k_eff_numerico_1g_hom = k_eff_numerico_1g;

%% ========================================================================
%  GENERACIÓN AVANZADA DE TABLA LATEX (Modos n=1, n=2, n=3) - 1G HOMOGÉNEO
%  ========================================================================

archivos = {'FEM_1G_HOM_N14_p1.mat','FEM_1G_HOM_N14_p2.mat','FEM_1G_HOM_N14_p3.mat'};
n_modos_tabla = 3;
n_archivos = length(archivos);

nombres_modos = {'Modo fundamental ($n=1$)', 'Primer armónico ($n=2$)', 'Segundo armónico ($n=3$)'};

fid = fopen('tabla_modos_1g_homogeneo.tex', 'w');
fprintf(fid, '\\begin{table}[h!]\n');
fprintf(fid, '\\centering\n');
fprintf(fid, '\\begin{tabular}{c c cccc}\n');
fprintf(fid, '\\toprule\n');

for n = 1:n_modos_tabla
    fprintf(fid, '\\multicolumn{6}{c}{\\textbf{%s}} \\\\ \\midrule\n', nombres_modos{n});
    fprintf(fid, '$N$ & $p$ & $N_{dof}(=N\\cdot p+1)$ & $k_{eff}$ & $\\epsilon_{k}$ (pcm) & $RMSE(\\phi)$ \\\\ \\midrule\n');

    for i = 1:n_archivos
        datos = load(archivos{i});
        N_elem = datos.malla.N;
        p_grado = datos.malla.grado_l;
        N_nodos = N_elem * p_grado + 1;

        k_num = datos.problema.keff(n);
        k_ana = datos.nu_sigma_f(2) / (datos.D(2) * (pi * n / datos.L)^2 + datos.sigma_a(2));

        phi_num = datos.problema.phi_inc(:, n);
        phi_num = phi_num(:) / max(abs(phi_num));

        x_nodos_tmp = datos.malla.x_nodos;
        phi_ana = sin(pi * n * x_nodos_tmp / datos.L);
        phi_ana = phi_ana(:);

        if phi_ana(2) * phi_num(2) < 0
            phi_ana = -phi_ana;
        end
        phi_ana = phi_ana / max(abs(phi_ana));

        err_k_pcm = abs(k_ana - k_num) * 1e5;
        rmse_val = rmse(phi_ana, phi_num);

        fprintf(fid, '%d & %d & %d & %.6f & %.2f & %.2e \\\\\n', ...
                N_elem, p_grado, N_nodos, k_num, err_k_pcm, rmse_val);
    end

    fprintf(fid, '\\midrule\n');
    fprintf(fid, '\\multicolumn{3}{c}{$k_{eff}$ analítico:} & %.6f & \\multicolumn{2}{c}{-} \\\\\n', k_ana);

    if n < n_modos_tabla
        fprintf(fid, '\\bottomrule\n\\vspace{0.1cm} \\\\\n');
    else
        fprintf(fid, '\\bottomrule\n');
    end
end

fprintf(fid, '\\end{tabular}\n');
fprintf(fid, '\\caption{Resultados de los 3 primeros modos comparados con el valor analítico.}\n');
fprintf(fid, '\\label{tab:modos_1g}\n');
fprintf(fid, '\\end{table}\n');
fclose(fid);
fprintf('\nTabla LaTeX con modos 1, 2 y 3 guardada en: tabla_modos_1g_homogeneo.tex\n');

%% ========================================================================
%  GENERACIÓN AVANZADA DE TABLA LATEX (Modos n=1, n=2, n=3) - 1G HETEROGÉNEO
%  ========================================================================

archivos = {'FEM_1G_HET_N14_p3.mat', 'FEM_1G_HET_N14_p4.mat', 'FEM_1G_HET_N14_p5.mat'};
n_modos_tabla = 3;
n_archivos = length(archivos);

nombres_modos = {'Modo fundamental ($n=1$)', 'Primer armónico ($n=2$)', 'Segundo armónico ($n=3$)'};
k_ana_het = [1.0621909229, 1.0521643703, 1.0358496988];

fid = fopen('tabla_modos_1g_heterogeneo.tex', 'w');
fprintf(fid, '\\begin{table}[h!]\n');
fprintf(fid, '\\centering\n');
fprintf(fid, '\\begin{tabular}{c c cccc}\n');
fprintf(fid, '\\toprule\n');

for n = 1:n_modos_tabla
    fprintf(fid, '\\multicolumn{6}{c}{\\textbf{%s}} \\\\ \\midrule\n', nombres_modos{n});
    fprintf(fid, '$N$ & $p$ & $N_{dof}(=N\\cdot p+1)$ & $k_{eff}$ & $\\epsilon_{k}$ (pcm) & $RMSE(\\phi)$ \\\\ \\midrule\n');

    for i = 1:n_archivos
        datos = load(archivos{i});
        N_elem = datos.malla.N;
        p_grado = datos.malla.grado_l;
        N_nodos = N_elem * p_grado + 1;

        k_num = datos.problema.keff(n);
        phi_num = datos.problema.phi_inc(:, n);
        phi_num = phi_num(:);

        [phi_a1, phi_a2, phi_a3] = calc_analitica_1g_het(datos.malla.x_nodos);
        phi_ana_cell = {phi_a1, phi_a2, phi_a3};
        phi_ana_nodos = phi_ana_cell{n};
        phi_ana_nodos = phi_ana_nodos(:);

        phi_num = phi_num / max(abs(phi_num));
        phi_ana_nodos = phi_ana_nodos / max(abs(phi_ana_nodos));

        idx_ref = max(2, floor(length(phi_num) / 5));
        if sign(phi_ana_nodos(idx_ref)) ~= sign(phi_num(idx_ref))
            phi_ana_nodos = -phi_ana_nodos;
        end

        err_k_pcm = abs(k_ana_het(n) - k_num) * 1e5;
        rmse_val = rmse(phi_ana_nodos, phi_num);

        fprintf(fid, '%d & %d & %d & %.6f & %.2f & %.2e \\\\\n', ...
                N_elem, p_grado, N_nodos, k_num, err_k_pcm, rmse_val);
    end

    fprintf(fid, '\\midrule\n');
    fprintf(fid, '\\multicolumn{3}{c}{$k_{eff}$ analítico:} & %.10f & \\multicolumn{2}{c}{-} \\\\\n', k_ana_het(n));

    if n < n_modos_tabla
        fprintf(fid, '\\bottomrule\n\\vspace{0.1cm} \\\\\n');
    else
        fprintf(fid, '\\bottomrule\n');
    end
end

fprintf(fid, '\\end{tabular}\n');
fprintf(fid, '\\caption{Resultados heterogéneos de los 3 primeros modos comparados con el valor analítico.}\n');
fprintf(fid, '\\label{tab:modos_1g_het}\n');
fprintf(fid, '\\end{table}\n');
fclose(fid);
fprintf('\nTabla LaTeX heterogénea guardada en: tabla_modos_1g_heterogeneo.tex\n');
%% =======================================================================
%% ========================================================================
%  GENERACIÓN AVANZADA DE TABLA LATEX (Modos n=1, n=2, n=3) - 2G HETEROGÉNEO
%  ========================================================================

archivos_2g_het = {'FEM_2G_HET_N14_p2.mat', 'FEM_2G_HET_N14_p3.mat', 'FEM_2G_HET_N14_p4.mat', 'FEM_2G_HET_N14_p5.mat'};
n_modos_tabla = 3;
n_archivos_2g = length(archivos_2g_het);

nombres_modos = {'Modo fundamental ($n=1$)', 'Primer armónico ($n=2$)', 'Segundo armónico ($n=3$)'};

datos_tmp = load(archivos_2g_het{1});
x_tmp = datos_tmp.malla.x_nodos;
[~, ~, ~, k_ana_2g_het] = calc_analitica_2g_het(x_tmp);

fprintf('\n--- k_eff analíticos 2G HET (calc_analitica_2g_het) ---\n');
for nn = 1:3
    fprintf('  Modo %d: k_eff = %.10f\n', nn, k_ana_2g_het(nn));
end

fid = fopen('tabla_modos_2g_heterogeneo.tex', 'w');
fprintf(fid, '\\begin{table}[h!]\n');
fprintf(fid, '\\centering\n');
fprintf(fid, '\\begin{tabular}{c c ccccc}\n');
fprintf(fid, '\\toprule\n');

for n = 1:n_modos_tabla
    fprintf(fid, '\\multicolumn{7}{c}{\\textbf{%s}} \\\\ \\midrule\n', nombres_modos{n});
    fprintf(fid, '$N$ & $p$ & $N_{dof}$ & $k_{eff}$ & $\\epsilon_{k}$ (pcm) & $RMSE(\\phi_1)$ & $RMSE(\\phi_2)$ \\\\ \\midrule\n');

    for i = 1:n_archivos_2g
datos = load(archivos_2g_het{i});
        N_elem = datos.malla.N;
        p_grado = datos.malla.grado_l;
        M_nod = N_elem * p_grado + 1;
        
        if length(datos.problema.keff) < n || size(datos.problema.phi_inc, 2) < n
            fprintf('  Aviso: %s solo tiene %d modos, saltando modo %d\n', archivos_2g_het{i}, length(datos.problema.keff), n);
            continue;
        end
        
        % 1. Extraer numérico
        k_num = datos.problema.keff(n);
        phi_num_g1 = datos.problema.phi_inc(1:M_nod, n);
        phi_num_g2 = datos.problema.phi_inc(M_nod+1:end, n);
        
        % Normalizar numérico
        max_num = max(max(abs(phi_num_g1)), max(abs(phi_num_g2)));
        phi_num_g1 = phi_num_g1 / max_num;
        phi_num_g2 = phi_num_g2 / max_num;

        % 2. Extraer analítico
        [phi_a1, phi_a2, phi_a3, ~] = calc_analitica_2g_het(datos.malla.x_nodos);
        phi_ana_cell = {phi_a1, phi_a2, phi_a3};
        phi_ana_modo = phi_ana_cell{n};
        
        % ¡EL PARCHE SALVAVIDAS!: 
        % Si la función devolvió un vector fila largo en lugar de matriz Nx2, lo partimos
        if size(phi_ana_modo, 2) ~= 2
            N_nodos = length(phi_ana_modo) / 2;
            phi_ana_g1 = phi_ana_modo(1:N_nodos);
            phi_ana_g2 = phi_ana_modo(N_nodos+1:end);
        else
            phi_ana_g1 = phi_ana_modo(:, 1);
            phi_ana_g2 = phi_ana_modo(:, 2);
        end
        
        % Forzar obligatoriamente a que sean vectores columna (Nx1)
        phi_ana_g1 = phi_ana_g1(:);
        phi_ana_g2 = phi_ana_g2(:);
        
        % Normalizar analítico
        max_ana = max(max(abs(phi_ana_g1)), max(abs(phi_ana_g2)));
        phi_ana_g1 = phi_ana_g1 / max_ana;
        phi_ana_g2 = phi_ana_g2 / max_ana;
        
        % 3. Corrección de signos (Ahora 100% segura)
        idx_ref = max(2, floor(M_nod / 5));
        if sign(phi_ana_g1(idx_ref)) ~= sign(phi_num_g1(idx_ref))
            phi_ana_g1 = -phi_ana_g1;
            phi_ana_g2 = -phi_ana_g2;
        end
        
        % 4. Calcular errores
        err_k_pcm = abs(k_ana_2g_het(n) - k_num) * 1e5;
        rmse_g1 = rmse(phi_ana_g1, phi_num_g1);
        rmse_g2 = rmse(phi_ana_g2, phi_num_g2);
        
        fprintf(fid, '%d & %d & %d & %.6f & %.2f & %.2e & %.2e \\\\\n', ...
                N_elem, p_grado, M_nod, k_num, err_k_pcm, rmse_g1, rmse_g2);    end

    fprintf(fid, '\\midrule\n');
    fprintf(fid, '\\multicolumn{3}{c}{$k_{eff}$ analítico:} & %.10f & \\multicolumn{3}{c}{-} \\\\\n', k_ana_2g_het(n));

    if n < n_modos_tabla
        fprintf(fid, '\\bottomrule\n\\vspace{0.1cm} \\\\\n');
    else
        fprintf(fid, '\\bottomrule\n');
    end
end

fprintf(fid, '\\end{tabular}\n');
fprintf(fid, '\\caption{Resultados 2G heterogéneos de los 3 primeros modos (Grupos 1 y 2).}\n');
fprintf(fid, '\\label{tab:modos_2g_het}\n');
fprintf(fid, '\\end{table}\n');
fclose(fid);
fprintf('\nTabla LaTeX 2G heterogénea guardada en: tabla_modos_2g_heterogeneo.tex\n');


%% ========================================================================
%  GENERACIÓN AVANZADA DE TABLA LATEX (Modos n=1, n=2, n=3) - 2G HOMOGÉNEO
%  ========================================================================

archivos_2g_hom = {'FEM_2G_HOM_N14_p2.mat', 'FEM_2G_HOM_N14_p3.mat', 'FEM_2G_HOM_N14_p4.mat', 'FEM_2G_HOM_N14_p5.mat'};
n_modos_tabla = 3;
n_archivos_2g_hom = length(archivos_2g_hom);

nombres_modos = {'Modo fundamental ($n=1$)', 'Primer armónico ($n=2$)', 'Segundo armónico ($n=3$)'};

% Obtener k_eff analíticos de la función
datos_tmp = load(archivos_2g_hom{4});
x_tmp = datos_tmp.malla.x_nodos;
[~, ~, ~, k_ana_2g_hom] = calc_analitica_2g_hom(x_tmp);

fprintf('\n--- k_eff analíticos 2G HOM (calc_analitica_2g_hom) ---\n');
for nn = 1:3
    fprintf('  Modo %d: k_eff = %.10f\n', nn, k_ana_2g_hom(nn));
end

% Abrir archivo para escribir
fid = fopen('tabla_modos_2g_homogeneo.tex', 'w');
fprintf(fid, '\\begin{table}[h!]\n');
fprintf(fid, '\\centering\n');
fprintf(fid, '\\begin{tabular}{c c ccccc}\n');
fprintf(fid, '\\toprule\n');

for n = 1:n_modos_tabla
    fprintf(fid, '\\multicolumn{7}{c}{\\textbf{%s}} \\\\ \\midrule\n', nombres_modos{n});
    fprintf(fid, '$N$ & $p$ & $N_{dof}$ & $k_{eff}$ & $\\epsilon_{k}$ (pcm) & $RMSE(\\phi_1)$ & $RMSE(\\phi_2)$ \\\\ \\midrule\n');

    for i = 1:n_archivos_2g_hom
        datos = load(archivos_2g_hom{i});
        N_elem = datos.malla.N;
        p_grado = datos.malla.grado_l;
        M_nod = N_elem * p_grado + 1;

        % Verificar que este archivo tiene suficientes modos calculados
        if length(datos.problema.keff) < n || size(datos.problema.phi_inc, 2) < n
            fprintf('  Aviso: %s solo tiene %d modos, saltando modo %d\n', archivos_2g_hom{i}, length(datos.problema.keff), n);
            continue;
        end

        k_num = datos.problema.keff(n);
        phi_num_g1 = datos.problema.phi_inc(1:M_nod, n);
        phi_num_g2 = datos.problema.phi_inc(M_nod+1:end, n);

        % Calcular solución analítica en los nodos de ESTA malla
        [phi_a1, phi_a2, phi_a3, ~] = calc_analitica_2g_hom(datos.malla.x_nodos);
        phi_ana_cell = {phi_a1, phi_a2, phi_a3};
        phi_ana_modo = phi_ana_cell{n}; % Matriz [Nnodos x 2]
        phi_ana_g1 = phi_ana_modo(:, 1);
        phi_ana_g2 = phi_ana_modo(:, 2);

        % Normalizar ambos (FEM y analítico) por su max absoluto global
        max_num = max(max(abs(phi_num_g1)), max(abs(phi_num_g2)));
        max_ana = max(max(abs(phi_ana_g1)), max(abs(phi_ana_g2)));
        phi_num_g1 = phi_num_g1 / max_num;
        phi_num_g2 = phi_num_g2 / max_num;
        phi_ana_g1 = phi_ana_g1 / max_ana;
        phi_ana_g2 = phi_ana_g2 / max_ana;

        % Alinear signos (usando G1 como referencia)
        idx_ref = max(2, min(floor(M_nod / 5), length(phi_ana_g1)));
        if sign(phi_ana_g1(idx_ref)) ~= sign(phi_num_g1(idx_ref))
            phi_ana_g1 = -phi_ana_g1;
            phi_ana_g2 = -phi_ana_g2;
        end

        % Calcular errores
        err_k_pcm = abs(k_ana_2g_hom(n) - k_num) * 1e5;
        rmse_g1 = rmse(phi_ana_g1, phi_num_g1);
        rmse_g2 = rmse(phi_ana_g2, phi_num_g2);

        fprintf(fid, '%d & %d & %d & %.6f & %.2f & %.2e & %.2e \\\\\n', ...
                N_elem, p_grado, M_nod, k_num, err_k_pcm, rmse_g1, rmse_g2);
    end

    fprintf(fid, '\\midrule\n');
    fprintf(fid, '\\multicolumn{3}{c}{$k_{eff}$ analítico:} & %.10f & \\multicolumn{3}{c}{-} \\\\\n', k_ana_2g_hom(n));

    if n < n_modos_tabla
        fprintf(fid, '\\bottomrule\n\\vspace{0.1cm} \\\\\n');
    else
        fprintf(fid, '\\bottomrule\n');
    end
end

fprintf(fid, '\\end{tabular}\n');
fprintf(fid, '\\caption{Resultados 2G homogéneos de los 3 primeros modos (Grupos 1 y 2).}\n');
fprintf(fid, '\\label{tab:modos_2g_hom}\n');
fprintf(fid, '\\end{table}\n');
fclose(fid);
fprintf('\nTabla LaTeX 2G homogénea guardada en: tabla_modos_2g_homogeneo.tex\n');



% =========================================================================
% FUNCIONES LOCALES (Deben ir estrictamente al final del script)
% =========================================================================

function d = calc_det_2g_het(k, a, b, L, pR, pC)
    [W, ~, ~, ~, ~, ~, ~] = build_matrix_2g_het(k, a, b, L, pR, pC);
    W_norm = W ./ max(abs(W), [], 1);
    d = real(det(W_norm));
end

function [W, V1, L1, V2, L2, V3, L3] = build_matrix_2g_het(k, a, b, L, pR, pC)
    a11_R = (pR.Sa1 + pR.Ss12 - pR.nSf1/k) / pR.D1;
    a12_R = -pR.nSf2 / (k * pR.D1);
    a21_R = -pR.Ss12 / pR.D2;
    a22_R = pR.Sa2 / pR.D2;
    M_R = [0 0 1 0; 0 0 0 1; a11_R a12_R 0 0; a21_R a22_R 0 0];
    a11_C = (pC.Sa1 + pC.Ss12 - pC.nSf1/k) / pC.D1;
    a12_C = -pC.nSf2 / (k * pC.D1);
    a21_C = -pC.Ss12 / pC.D2;
    a22_C = pC.Sa2 / pC.D2;
    M_C = [0 0 1 0; 0 0 0 1; a11_C a12_C 0 0; a21_C a22_C 0 0];
    [V1, Lam1] = eig(M_R); L1 = diag(Lam1);
    [V2, Lam2] = eig(M_C); L2 = diag(Lam2);
    [V3, Lam3] = eig(M_R); L3 = diag(Lam3);
    W = zeros(12, 12);
    for i=1:4
        W(1,i) = V1(1,i); W(2,i) = V1(2,i);
    end
    for i=1:4
        E1 = exp(L1(i)*a); E2 = 1.0;
        W(3,i) = V1(1,i)*E1;         W(3,i+4) = -V2(1,i)*E2;
        W(4,i) = V1(2,i)*E1;         W(4,i+4) = -V2(2,i)*E2;
        W(5,i) = pR.D1*V1(3,i)*E1;   W(5,i+4) = -pC.D1*V2(3,i)*E2;
        W(6,i) = pR.D2*V1(4,i)*E1;   W(6,i+4) = -pC.D2*V2(4,i)*E2;
    end
    for i=1:4
        E2 = exp(L2(i)*(b-a)); E3 = 1.0;
        W(7,i+4) = V2(1,i)*E2;        W(7,i+8) = -V3(1,i)*E3;
        W(8,i+4) = V2(2,i)*E2;        W(8,i+8) = -V3(2,i)*E3;
        W(9,i+4) = pC.D1*V2(3,i)*E2;  W(9,i+8) = -pR.D1*V3(3,i)*E3;
        W(10,i+4) = pC.D2*V2(4,i)*E2; W(10,i+8) = -pR.D2*V3(4,i)*E3;
    end
    for i=1:4
        E3 = exp(L3(i)*(L-b));
        W(11,i+8) = V3(1,i)*E3;
        W(12,i+8) = V3(2,i)*E3;
    end
end

function d = calc_det_2g_hom(k, L, D1, D2, sigma_a1, sigma_a2, nu_sigma_f1, nu_sigma_f2, sigma_s12)
    a11 = (sigma_a1 + sigma_s12 - nu_sigma_f1/k) / D1;
    a12 = -nu_sigma_f2 / (k * D1);
    a21 = -sigma_s12 / D2;
    a22 = sigma_a2 / D2;
    M = [0,0,1,0; 0,0,0,1; a11,a12,0,0; a21,a22,0,0];
    [V, Lambda] = eig(M);
    lambda_vals = diag(Lambda);
    W = zeros(4, 4);
    for i = 1:4
        v1 = V(1,i); v2 = V(2,i);
        W(1,i) = v1; W(2,i) = v2;
        W(3,i) = v1 * exp(lambda_vals(i)*L);
        W(4,i) = v2 * exp(lambda_vals(i)*L);
    end
    d = det(W);
end

tiempo = toc;
fprintf('\nTiempo total de ejecución: %.2f segundos\n', tiempo);
