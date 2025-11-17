%% Randles + (opcional) Warburg (semi-infinito) con concentraciones en mg/dL (glucosa)
%  ZW = sigma / sqrt(j*w)
%  sigma = (R*T) / (n^2*F^2*A*C*sqrt(2*D))
%  Rct   = (R*T) / (k0*n^2*F^2*A*C)   [opcional con k0]
clear; clc; close all;

%% ---------------- FISICOQUÍMICA ----------------
T = 298.15;                 % K
n = 1;                      % #
F = 96485.33212;            % C/mol
R = 8.314462618;            % J/mol/K

A = 0.1;                    % cm^2  (área electroactiva)
D = 1e-5;                   % cm^2/s

%% ---------------- CONCENTRACIONES (tu formato) ----------------
concs_mg_dL = [50 100 200 400];   % típico clínico (mg/dL)

% Conversión a mM usando tu factor:
% 1 mM glucosa ≈ 18.06 mg/dL  (PM ≈ 180.16 g/mol)
mgdL_to_mM = @(x) x ./ 18.06;
concs_mM = mgdL_to_mM(concs_mg_dL);

% mM -> mol/cm^3  (1 mM = 1e-3 mol/L = 1e-6 mol/cm^3)
mM_to_mol_cm3 = @(c_mM) c_mM * 1e-6;

%% ---------------- ELEMENTOS DE CIRCUITO ----------------
Rs  = 50;                    % ohm
Cdl = 2e-6;                  % F

% ¿Rct(C) con k0 (ec. 57) o fijo?
USE_k0_FOR_Rct = true;       % <-- cambia a false para Rct fijo
k0  = 1e-3;                  % cm/s
Rct_fixed = 1000;            % ohm (usado si USE_k0_FOR_Rct=false)

% *** NUEVO: usar o no la difusión de Warburg ***
USE_WARBURG = false;          % <-- pon false para Randles sin Warburg

%% ---------------- FRECUENCIAS ----------------
fmin = 1e-1; fmax = 1e5; Npts = 80;
f  = logspace(log10(fmin), log10(fmax), Npts).';
w  = 2*pi*f; jw = 1j*w;

%% ---------------- SIMULACIÓN ----------------
Z_all = zeros(Npts, numel(concs_mg_dL));
leg   = strings(1, numel(concs_mg_dL));

for k = 1:numel(concs_mg_dL)
    C_mg_dL = concs_mg_dL(k);
    C_mM    = concs_mM(k);
    C_cm3   = mM_to_mol_cm3(C_mM);      % mol/cm^3

    % Rct
    if USE_k0_FOR_Rct && C_cm3 > 0
        Rct = (R*T) / (k0 * n^2 * F^2 * A * C_cm3);        % ohm (depende de C)
    else
        Rct = Rct_fixed;                                    % ohm (fijo o C=0)
    end

    % Warburg (opcional)
    if USE_WARBURG && C_cm3 > 0
        sigma = (R*T) / (n^2 * F^2 * A * C_cm3 * sqrt(2*D)); % ohm*s^(-1/2)
        Zw = sigma ./ sqrt(jw);
    else
        Zw = zeros(size(w));                                % sin difusión
    end

    % Casos límite si C=0: sin especie redox => sólo Cdl (y Rct→infinito idealmente)
    if C_cm3 <= 0
        Z = Rs + 1./(jw*Cdl);
        Z_all(:,k) = Z;
        leg(k) = sprintf('C = %g mg/dL (sin redox)', C_mg_dL);
        continue
    end

    % Rama faradaica: ZF = Rct + Zw
    ZF = Rct + Zw;

    % Paralelo con Cdl
    Ypar = jw.*Cdl + 1./ZF;

    % Impedancia total
    Z = Rs + 1./Ypar;

    Z_all(:,k) = Z;
    if USE_WARBURG
        leg(k) = sprintf('C = %g mg/dL (%.3g mM) + W', C_mg_dL, C_mM);
    else
        leg(k) = sprintf('C = %g mg/dL (%.3g mM) sin W', C_mg_dL, C_mM);
    end
end

%% ---------------- GRÁFICAS ----------------

% --- Nyquist ---
figure('Color','w'); hold on; grid on; axis equal;
for k = 1:numel(concs_mg_dL)
    plot(real(Z_all(:,k)), -imag(Z_all(:,k)), 'LineWidth',3);
end
xlabel('Z'' (\Omega)','FontSize',14); ylabel('-Z'''' (\Omega)','FontSize',14);
title('Randles (sin Warburg): efecto de la glucosa (mg/dL)','FontSize',25);
legend(leg, 'Location','best');

% --- Bode |Z| (log-log forzado) ---
figure('Color','w');

hold on; grid on;
for k = 1:numel(concs_mg_dL)
    loglog(f, abs(Z_all(:,k)), 'LineWidth',1.6);
end
% Fuerza explícitamente ambas escalas en log
set(gca, 'XScale','log', 'YScale','log');
xlim([min(f(f>0)) max(f)]);           % evita problemas si hubiera f=0
xlabel('f (Hz)','FontSize',16); ylabel('|Z| (\Omega)','FontSize',16);
title('Bode – Magnitud','FontSize',20);
legend(leg, 'Location','best');

% --- Bode Fase (semilog-x forzado) ---
figure();
hold on; grid on;
for k = 1:numel(concs_mg_dL)
    % Desenvuelve y convierte a grados
    phi_deg = unwrap(angle(Z_all(:,k)))*180/pi;
    semilogx(f, phi_deg, 'LineWidth',1.6);
end
set(gca, 'XScale','log');             % asegura eje X en log
xlim([min(f(f>0)) max(f)]);
xlabel('f (Hz)','FontSize',16); ylabel('\angle Z (deg)','FontSize',16);
title('Bode – Fase','FontSize',20);
legend(leg, 'Location','best');


%% ---------------- LECTURA RÁPIDA ----------------
fprintf('Ru=%.3g ohm, Cdl=%.3g F, A=%.3g cm^2, D=%.3g cm^2/s\n', Rs, Cdl, A, D);
if USE_k0_FOR_Rct
    fprintf('Rct(C) calculado con k0=%.3g cm/s\n', k0);
else
    fprintf('Rct fijo=%.3g ohm\n', Rct_fixed);
end
for k = 1:numel(concs_mg_dL)
    C_mg_dL = concs_mg_dL(k);
    C_mM    = concs_mg_dL(k)/18.06;
    C_cm3   = mM_to_mol_cm3(C_mM);
    if C_cm3<=0
        fprintf('C=%g mg/dL -> sin rama faradaica (sólo Cdl)\n', C_mg_dL);
    else
        if USE_WARBURG
            sigma = (R*T)/(n^2*F^2*A*C_cm3*sqrt(2*D));
        else
            sigma = 0;
        end
        if USE_k0_FOR_Rct
            Rct = (R*T)/(k0*n^2*F^2*A*C_cm3);
        else
            Rct = Rct_fixed;
        end
        fp = 1/(2*pi*Rct*Cdl);
        fprintf('C=%g mg/dL (%.3g mM): sigma=%.3g ohm*s^-1/2, Rct=%.3g ohm, fp~%.3g Hz\n', ...
                 C_mg_dL, C_mM, sigma, Rct, fp);
    end
end
 