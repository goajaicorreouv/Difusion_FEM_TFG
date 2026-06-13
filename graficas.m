%% 1. GRÁFICO DE CONVERGENCIA LOG-LOG (Ejemplo conceptual)
% Para esto, asume que has corrido tu FEM para distintas mallas y 
% has guardado los resultados en estos vectores:

N_dof_list = [10, 20, 40, 80, 160]; % Número de nodos/grados de libertad
L_reactor = 350; % Longitud total del reactor en cm
h_vals = L_reactor ./ (N_dof_list - 1); % Tamaño característico del elemento

% Supongamos que estos son los RMSE que obtuviste del Modo 1 para cada malla
rmse_vals = [1.2e-2, 3.1e-3, 7.8e-4, 1.9e-4, 4.8e-5]; 

% Calcular el orden de convergencia (pendiente de la recta log-log)
p_fit = polyfit(log(h_vals), log(rmse_vals), 1);
orden_p = p_fit(1);

figure('Name', 'Convergencia Espacial');
loglog(h_vals, rmse_vals, 'b-o', 'LineWidth', 2, 'MarkerFaceColor', 'b');
grid on;
xlabel('Tamaño del elemento h (cm)');
ylabel('Error RMSE del Flujo (\phi)');
title(sprintf('Análisis de Convergencia FEM (Orden p \\approx %.2f)', orden_p));

% Añadir una línea de referencia teórica (O(h^2))
hold on;
h_ref = linspace(min(h_vals), max(h_vals), 50);
error_ref = (rmse_vals(end)/h_vals(end)^2) * h_ref.^2; % Ajuste visual
loglog(h_ref, error_ref, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Referencia O(h^2)');
legend('Error FEM', 'Pendiente Teórica', 'Location', 'best');


%% 2. DISTRIBUCIÓN ESPACIAL DEL ERROR (Caso 1G Heterogéneo - Modo 1)
% Asume que ya ejecutaste tu bucle y tienes phi_analitico_nodos y phi_numerico
modo_plot = 1; % Analizamos el modo fundamental

% Usamos la normalización que ya tienes en tu código
phi_ana_norm = phi_analitico_het{modo_plot}; 
phi_num_norm = problema.phi_inc(:, modo_plot) / max(abs(problema.phi_inc(:, modo_plot)));

% Interpolar la analítica fina a los nodos de la malla numérica
phi_ana_interp = interp1(linspace(0, L, 1000)', phi_ana_norm, x_nodos_het, 'linear', 'extrap');

% Calcular el error absoluto local
error_espacial = abs(phi_num_norm - phi_ana_interp);

figure('Name', 'Error Espacial Heterogéneo');
plot(x_nodos_het, error_espacial, 'r-', 'LineWidth', 2);
hold on; grid on;

% Marcar las interfaces físicas
xline(a, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Interfaz Reflector-Núcleo');
xline(b, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Interfaz Núcleo-Reflector');

xlabel('x (cm)');
ylabel('| \phi_{FEM}(x) - \phi_{Ana}(x) |');
title(sprintf('Distribución del Error Absoluto Local (Modo %d)', modo_plot));
legend('Location', 'best');

%% 3. RATIO DE FLUJOS (Solo aplicable si tienes resultados 2G)
% Asumiendo que tienes phi_numerico_g1 (Rápido) y phi_numerico_g2 (Térmico)
% y que la malla tiene interfaces en x=a y x=b

% Simulación de variables (Reemplaza con tus variables reales de 2G)
% phi_numerico_g1 = problema.phi_inc(1:M_nodos, 1);
% phi_numerico_g2 = problema.phi_inc(M_nodos+1:end, 1);

% Ratio de flujo: Térmico / Rápido
% Añadimos eps para evitar división por cero en los bordes
ratio_espectro = phi_numerico_g2 ./ (phi_numerico_g1 + eps);

figure('Name', 'Espectro Neutrónico');
plot(x_nodos_het, ratio_espectro, 'g-', 'LineWidth', 2);
hold on; grid on;

% Marcar interfaces
xline(a, 'k--', 'LineWidth', 1.5);
xline(b, 'k--', 'LineWidth', 1.5);

xlabel('x (cm)');
ylabel('Ratio de Flujo (\phi_{termico} / \phi_{rapido})');
title('Evolución Espacial del Espectro Neutrónico');

%% 4. BARRIDO DEL DETERMINANTE PARA HALLAR MODOS (1G Heterogéneo)
% Rango de k_eff a explorar (de 0.1 a 1.5)
k_test = linspace(0.1, 1.5, 500);
det_vals = zeros(size(k_test));

for i = 1:length(k_test)
    k_t = k_test(i);
    
    % Evitar valores complejos si k_t es muy pequeño en el núcleo
    if (nu_sigma_f2/k_t - sigma_a2) < 0
        det_vals(i) = NaN;
        continue;
    end
    
    k1 = sqrt(sigma_a1 / D1);
    k2 = sqrt((nu_sigma_f2/k_t - sigma_a2) / D2);
    
    % Matriz de condiciones (tu mismo código, simplificado)
    M11 = sinh(k1*a); M12 = -sin(k2*a); M13 = -cos(k2*a); M14 = 0;
    M21 = D1*k1*cosh(k1*a); M22 = -D2*k2*cos(k2*a); M23 = D2*k2*sin(k2*a); M24 = 0;
    M31 = 0; M32 = sin(k2*b); M33 = cos(k2*b); M34 = -sinh(k1*(L-b));
    M41 = 0; M42 = D2*k2*cos(k2*b); M43 = -D2*k2*sin(k2*b); M44 = D1*k1*cosh(k1*(L-b));
    
    Matriz = [M11 M12 M13 M14; M21 M22 M23 M24; M31 M32 M33 M34; M41 M42 M43 M44];
    
    % Guardamos el determinante (escala logarítmica para visualizar mejor si es enorme)
    det_vals(i) = det(Matriz);
end

figure('Name', 'Búsqueda de Modos');
plot(k_test, det_vals, 'b-', 'LineWidth', 1.5);
hold on; grid on;
yline(0, 'r-', 'LineWidth', 1.5); % Línea del cero

% Marcar los k_eff numéricos que ya encontraste
for m_i = 1:m
    xline(k_eff_numerico_1g_het(m_i), 'g--', 'LineWidth', 1.5, ...
          'DisplayName', sprintf('Modo %d FEM', m_i));
end

xlabel('k_{eff} de prueba');
ylabel('Valor del Determinante');
title('Búsqueda Analítica de Autovalores (Ceros del Determinante)');
ylim([-max(abs(det_vals(~isnan(det_vals))))/10, max(abs(det_vals(~isnan(det_vals))))/10]); % Acotar el eje Y


%% 1. Gráfica: Barrido del Determinante (Búsqueda de Autovalores)
% Requiere: 
% - 'k_vec': vector de valores de k_eff probados (e.g. linspace(0.1, 1.5, 500))
% - 'det_vals': vector con el valor del determinante para cada 'k_vec'
% - 'problema.keff': autovalor(es) encontrado(s) por el método

% --- Configuración de Gráfica ---
figure('Name', 'Barrido del Determinante', 'Position', [100, 100, 700, 450]);
hold on; grid on; box on;

% Graficar la curva del determinante
% Usamos log o limitamos los ejes porque estas funciones suelen explotar exponencialmente
plot(k_vec, det_vals, 'k-', 'LineWidth', 1.5, 'DisplayName', '$\det(M(k_{eff}))$');

% Línea base de referencia en Y=0
yline(0, 'r--', 'LineWidth', 1.2, 'DisplayName', 'Cero ($y=0$)');

% Marcar las raíces halladas
if isfield(problema, 'keff')
    % Se asume que en keff el determinante es 0
    plot(problema.keff, zeros(size(problema.keff)), 'bo', 'MarkerSize', 8, ...
        'MarkerFaceColor', 'b', 'DisplayName', 'Autovalores $k_{eff}$');
end

% Acotación CLAVE del eje Y: evitamos que las colas exponenciales oculten las raíces.
% Calculamos un umbral basado en el rango intercuartílico u orden de magnitud
rango_y = prctile(abs(det_vals), 85); % Ajusta este percentil (ej. 80-95) si no cruza bien
ylim([-rango_y, rango_y]); 
xlim([min(k_vec), max(k_vec)]);

% Etiquetas con LaTeX
xlabel('Valores de prueba de $k_{eff}$', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Determinante de condiciones de contorno', 'Interpreter', 'latex', 'FontSize', 14);
title('\textbf{Aproximación de Autovalores: Barrido del Determinante}', 'Interpreter', 'latex', 'FontSize', 15);
legend('Location', 'best', 'Interpreter', 'latex', 'FontSize', 12);

hold off;


%% 2. Gráfica: Convergencia de Malla Log-Log
% Requiere:
% - 'h_vals': vector con los tamaños de malla (L / x_nodos) de las pruebas
% - 'error_vals': vector con la norma del error (e.g. norma L2 o RMSE)

% --- Configuración de Gráfica ---
figure('Name', 'Analisis de Convergencia', 'Position', [150, 150, 700, 450]);
hold on; grid on; box on; 
set(gca, 'XScale', 'log', 'YScale', 'log'); % Configuración Log-Log directa

% 1. Ajuste lineal Log-Log para encontrar pendiente (p) y constante (C)
% log(E) = p*log(h) + log(C)
fit_params = polyfit(log(h_vals), log(error_vals), 1);
pendiente_p = fit_params(1); 
C_fit = exp(fit_params(2));

% 2. Generar recta teórica de referencia (Asumiendo O(h^2) típico de FEM lineal)
p_teorico = 2; % Cambia a 3 si usas polinomios de base cuadráticos
% Para que se dibuje paralela y cercana sin estorbar los datos principales:
C_teorico = error_vals(1) / (h_vals(1)^p_teorico) * 0.8; 
error_teorico = C_teorico * h_vals.^p_teorico;

% Trazados
loglog(h_vals, error_vals, 's-', 'Color', [0, 0.447, 0.741], 'LineWidth', 2, ...
    'MarkerSize', 8, 'MarkerFaceColor', [0, 0.447, 0.741], ...
    'DisplayName', 'Error Calculado $||\phi_h - \phi_{an}||$');

loglog(h_vals, error_teorico, 'k--', 'LineWidth', 1.5, ...
    'DisplayName', sprintf('Referencia Teórica $\\mathcal{O}(h^{%d})$', p_teorico));

% Estética y Etiquetas
xlabel('Tama\~no característico del elemento $h$ (cm)', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Error RMS', 'Interpreter', 'latex', 'FontSize', 14);
% Insertamos orden calculado en el propio título
title(sprintf('\\textbf{Convergencia Espacial: Orden Calculado} $\\mathbf{p \\approx %.2f}$', pendiente_p), ...
      'Interpreter', 'latex', 'FontSize', 15);
legend('Location', 'best', 'Interpreter', 'latex', 'FontSize', 12);
xlim([min(h_vals)*0.7, max(h_vals)*1.4]); % Márgenes extra para que luzca bien en log

hold off;


%% 3. Gráfica: Error Espacial (FEM vs Analítica)
% Requiere:
% - 'malla.x_nodos': Vector espacial (posiciones)
% - 'phi_FEM': Flujo escalar resuelto
% - 'phi_analitico': Flujo puramente analítico en los mismos nodos
% - 'malla.L': Longitud total
% - 'pos_interfaces': (Opcional) Array con las coordenadas x donde cambia el material 
%                     e.g. pos_interfaces = [30, 60];

error_local = abs(phi_FEM - phi_analitico);

% --- Configuración de Gráfica ---
figure('Name', 'Error Espacial', 'Position', [200, 200, 750, 400]);
hold on; grid on; box on;

% Gráfica principal del error 
plot(malla.x_nodos, error_local, '-', 'Color', [0.850, 0.325, 0.098], ...
     'LineWidth', 1.5, 'DisplayName', '$|\phi_{FEM}(x) - \phi_{an}(x)|$');

% Representar interfaces de materiales (Líneas verticales)
if exist('pos_interfaces', 'var') && ~isempty(pos_interfaces)
    y_limits = ylim; % Coger límites para cubrir todo el alto
    for idx_i = 1:length(pos_interfaces)
        % Truco: Solo añadir a la leyenda la primera línea para que no se duplique
        if idx_i == 1
            plot([pos_interfaces(idx_i), pos_interfaces(idx_i)], y_limits, 'k-.', ...
                 'LineWidth', 1.2, 'DisplayName', 'Interfase de Material');
        else
            plot([pos_interfaces(idx_i), pos_interfaces(idx_i)], y_limits, 'k-.', ...
                 'LineWidth', 1.2, 'HandleVisibility', 'off');
        end
    end
end

% Ajuste de visor
xlim([0, malla.L]);
ylim([0, max(error_local)*1.1]); % Dejar en 0 como mínimo ya que es valor absoluto

% Textos
xlabel('Coordenada espacial $x$ (cm)', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Error Local Absoluto', 'Interpreter', 'latex', 'FontSize', 14);
title('\textbf{Distribución Espacial del Error de Aproximación}', 'Interpreter', 'latex', 'FontSize', 15);
legend('Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 12);

hold off;


%% 4. Gráfica: Espectro Neutrónico (Solo 2 Grupos)
% Requiere:
% - 'malla.x_nodos': Vector espacial
% - 'phi_g1_FEM': Flujo rápido
% - 'phi_g2_FEM': Flujo térmico
% - 'pos_interfaces': Vector con fronteras del núcleo (e.g. fronteras del combustible)

% Calculamos ratio evitando divisiones por ceros exactos con (eps)
ratio_espectro = phi_g2_FEM ./ (phi_g1_FEM + eps);

% --- Configuración de Gráfica ---
figure('Name', 'Espectro Neutronico', 'Position', [250, 250, 750, 400]);
hold on; grid on; box on;

y_lims = [0, max(ratio_espectro) * 1.1];

% Sombrear la región del Núcleo Combustible (asumimos que la pos_interfaces las delimitan)
% Esto requiere un array "pos_interfaces = [inicio_nucleo, fin_nucleo]"
if exist('pos_interfaces', 'var') && length(pos_interfaces) >= 2
    % patch crea una región rectangular sombreada en el fondo
    patch([pos_interfaces(1), pos_interfaces(2), pos_interfaces(2), pos_interfaces(1)], ...
          [y_lims(1), y_lims(1), y_lims(2), y_lims(2)], ...
          [0.9, 0.9, 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.6, ...
          'DisplayName', 'Zonas: Reflector / Núcleo / Reflector'); 
end

% Curva de termalización
plot(malla.x_nodos, ratio_espectro, '-', 'Color', [0.494, 0.184, 0.556], 'LineWidth', 2, ...
     'DisplayName', 'Ratio $\phi_{1,term}/\phi_{2,rap}$');

% Líneas de interfaz
if exist('pos_interfaces', 'var') && ~isempty(pos_interfaces)
    for idx_i = 1:length(pos_interfaces)
        plot([pos_interfaces(idx_i), pos_interfaces(idx_i)], y_lims, 'k--', ...
             'LineWidth', 1, 'HandleVisibility', 'off');
    end
end

% Ajuste de visor
xlim([0, malla.L]);
ylim(y_lims);

xlabel('Coordenada espacial $x$ (cm)', 'Interpreter', 'latex', 'FontSize', 14);
ylabel('Cociente $\phi_2(x) / \phi_1(x)$', 'Interpreter', 'latex', 'FontSize', 14);
title('\textbf{Espectro Neutrónico Espacial: Termalización en el Reflector}', 'Interpreter', 'latex', 'FontSize', 15);
legend('Location', 'northeast', 'Interpreter', 'latex', 'FontSize', 12);

hold off;
