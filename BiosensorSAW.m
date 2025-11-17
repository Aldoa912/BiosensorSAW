%% SAW delay line (Impulse Response model) + Ventana sensora Randles
% Dos IDTs idénticos (Tx/Rx) y una ventana sensora eléctrica (Randles)
%
% IL se calcula con Y_total = Ga + j(ωCt + Ba) + eta_sens*Ysens(f,c),
% donde Ysens es la admitancia del circuito de Randles evaluada en f_RF.
%

clear; clc; close all;

%% ================= PARTE SAW (base) ====================================
% -------- Parameters (EDIT HERE) ---------------------------------------
wavelength = 1e-3;                 % [m]
v_saw   = 3290;                    % [m/s]
f0      = v_saw/wavelength;        % [Hz]
NBW     = 5.5e5;                   % [Hz]
k       = sqrt(0.0075);            % [-]
Rg      = 50;                      % [Ohm]
Rin     = 50;                      % [Ohm]

% --------- Cs from dielectric constants: Eqs. (36)-(37) ----------------
eps0  = 8.854187817e-12;           
% <<< PON AQUÍ LOS VALORES DE TU TABLA >>>
eps11_rel = 53.6;
eps33_rel = 43.4;
eps13_rel = 0.0;

eps_r = sqrt( (eps11_rel/eps33_rel) - (eps13_rel^2)/(eps33_rel^2) );   % (37)
Cs_base  = 20;                     % [pF/cm] típico en literatura SAW
denom    = 0.5 * ( 1 + eps_r * (eps33_rel) );
Cs_pF_per_cm_calc = Cs_base * denom;     % [pF/cm] (36)
unit_conv = (1e-12) / (1e-2);            % pF/cm -> F/m
Cs = Cs_pF_per_cm_calc * unit_conv;      % [F/m]
fprintf('Cs (calc, eqs.36-37) = %.6f pF/cm  ->  %.3e F/m\n', ...
        Cs_pF_per_cm_calc, Cs);

lambda = v_saw / f0;                      % [m]
Np     = round((f0*2)/(NBW));             % pares de dedos (heurística)
Wa     = Wa_optimized(f0,Np,Cs,k,Rin);    % (10) apertura [m]
Ct     = Np * Wa * Cs;                    % (7) capacitancia estática [F]

% Sweep de frecuencia alrededor de f0
fmin = 0.5*f0;  fmax = 1.5*f0;  Npts = 8001;
f = linspace(fmin,fmax,Npts).';

% Núcleo SAW
Hf    = H_resp(f,f0,Np,Cs,k);                           
Hf0   = H_resp(f0,f0,Np,Cs,k);
FRdB  = 20*log10( abs((Hf.^2) ./ (Hf0.^2)) );           

Gaf   = Ga(f,f0,Np,Cs,k,Wa);                            
Baf   = Ba(f,f0,Np,Cs,k,Wa);                            

Y_saw = Gaf + 1j*(2*pi*f*Ct + Baf);                     
Zf    = 1 ./ Y_saw;                                     

% IL base (sin sensor): usando Re/Im(Y_saw)
ILdB_base = -10*log10( (2*Gaf*Rg) ./ ( (1 + Rg*real(Y_saw)).^2 + (Rg*imag(Y_saw)).^2 ) );
IL_base   = (2*Gaf*Rg) ./ ( (1 + Rg*real(Y_saw)).^2 + (Rg*imag(Y_saw)).^2 );

% Scalars
lambda_um = lambda*1e6; Wa_um = Wa*1e6; Ct_pF = Ct*1e12;
[ minIL_base,idxb ] = min(ILdB_base); f_at_minIL_base = f(idxb);

fprintf('\n<<< SAW (sin sensor) >>>\n');
fprintf('Wavelength                = %.3f µm\n', lambda_um);
fprintf('f0 (synchronous)          = %.6e Hz\n', f0);
fprintf('Np (finger pairs)         = %d\n', Np);
fprintf('Aperture (Wa)             = %.1f µm\n', Wa_um);
fprintf('Ct (IDT)                  = %.3f pF\n', Ct_pF);
fprintf('Min IL (base)             = %.3f dB @ %.6e Hz\n', minIL_base, f_at_minIL_base);

%% ================= PARTE SENSOR (ventana Randles) ======================
% Randles clásico con Cdl ideal y (opcional) Warburg semi-infinito:
%  Zw = sigma / sqrt(jw)
%  sigma = (R*T) / (n^2*F^2*A*C*sqrt(2*D))
%  Rct   = (R*T) / (k0*n^2*F^2*A*C)

% --- Fisicoquímica ----------------
T = 298.15;                 % K
n = 1;                      % #
F = 96485.33212;            % C/mol
R = 8.314462618;            % J/mol/K

A_cm2 = 0.10;               % cm^2  (área electroactiva)
D     = 1e-5;               % cm^2/s (coef. difusión)

% ------------------ CONCENTRACIONES (mg/dL) ----------------------------
concs_mg_dL = [210];                  % puedes pasar un vector: [90 120 150 180 210]

% Conversión a mM:
% 1 mM glucosa ≈ 18.06 mg/dL  (PM ≈ 180.16 g/mol)
mgdL_to_mM = @(x) x./18.06;
concs_mM = mgdL_to_mM(concs_mg_dL);

% mM -> mol/cm^3  (1 mM = 1e-3 mol/L = 1e-6 mol/cm^3)
mM_to_mol_cm3 = @(c_mM) c_mM * 1e-6;

% ------------------ Rct físico (Butler–Volmer simplificado) ------------
% Constantes físicas para Rct = RT/(k0 n^2 F^2 A C)
k0_cms = 1e-3;             % cm/s (constante de velocidad heterogénea) <-- AJUSTA
C_floor = 1e-12;           % mol/cm^3 (piso numérico muy pequeño)

Rct_phys = @(C_mol_cm3) (R*T) ./ (k0_cms * (n^2) * F^2 * A_cm2 .* max(C_mol_cm3, C_floor));
Rct_of_c = @(c_mM) Rct_phys( mM_to_mol_cm3(c_mM) );

% --- Flags de modelo electroquímico
USE_WARBURG    = false;   % true: incluye difusión de Warburg dependiente de C

% --- Fracción de apertura cubierta por el sensor (0..1)
eta_sens = 0.5;           % ej., 50% de la apertura Wa está "cargada"

% --- Cdl (único modelo, ya no se usa CPE)
Cdl    = 2e-6;            % [F] (ajústalo según tu sistema)

% --- Resistencia serie opcional para no shuntear el IDT completo
INCLUDE_Rs_SERIE = false;
Rs_elec = 50;            % Ω

% --- Reservas y fase base para comparación -----------------------------
ILdB_sense     = zeros(numel(f), numel(concs_mM));
IL_sense_all   = zeros(numel(f), numel(concs_mM));
S21_phase_deg  = zeros(numel(f), numel(concs_mM));
f_min_IL       = zeros(numel(concs_mM),1);
IL_min         = zeros(numel(concs_mM),1);

Y_total_base = Y_saw;
S21_base = (Gaf) ./ ( (1 + Rg*real(Y_total_base)) + 1j*(Rg*imag(Y_total_base)) );
phi_base_deg = rad2deg(unwrap(angle(S21_base)));

omega = 2*pi*f;
jw    = 1j*omega;

% --- IL con sensor para cada concentración (Randles + Warburg opcional) -
for ic = 1:numel(concs_mM)
    c_mM   = concs_mM(ic);
    C_cm3  = mM_to_mol_cm3(c_mM);          % mol/cm^3

    % Rct físico desde Butler–Volmer
    Rct_k  = Rct_of_c(c_mM);               % Ω (depende de C)

    % Warburg dependiente de concentración (Zw = sigma / sqrt(jw))
    if USE_WARBURG && (C_cm3 > 0)
        sigma_k = (R*T) / ( (n^2) * F^2 * A_cm2 * C_cm3 * sqrt(2*D) ); % Ω·s^(-1/2)
        Zw_k    = sigma_k ./ sqrt(jw);                                  % Ω
        Yw_k    = 1 ./ Zw_k;                                            % S
    else
        Yw_k    = zeros(size(jw));                                      % sin difusión
    end

    % Doble capa (Cdl ideal)
    Ydl_k = jw .* Cdl;          % S

    % Núcleo Randles: Y = jωCdl + 1/Rct + Yw
    Ys_k = Ydl_k + (1./Rct_k) + Yw_k;   % S

    % Rs eléctrico en serie (opcional)
    if INCLUDE_Rs_SERIE
        Ys_k = 1 ./ ( Rs_elec + 1./Ys_k );
    end

    % Admitancia total en paralelo (SAW + ventana sensora)
    Y_total = Gaf + 1j*(2*pi*f*Ct + Baf) + eta_sens*(Ys_k);

    % IL (lineal y dB) con sensor
    IL_sense_k = (2*Rg*Gaf) ./ ( (1 + Rg*real(Y_total)).^2 + (Rg*imag(Y_total)).^2 );
    ILdB_k     = -10*log10(IL_sense_k);

    ILdB_sense(:,ic)   = ILdB_k;
    IL_sense_all(:,ic) = IL_sense_k;

    % Fase de S21 (en grados), con desenrollado
    S21_k = (Gaf) ./ ( (1 + Rg*real(Y_total)) + 1j*(Rg*imag(Y_total)) );
    S21_phase_deg(:,ic) = rad2deg(unwrap(angle(S21_k)));

    % Mínimo de IL y su frecuencia
    [IL_min(ic), idx] = min(ILdB_k);
    f_min_IL(ic)      = f(idx);
end

%% ==================== GRÁFICA IL NORMALIZADA + ZOOM ======================

% Figura ya maximizada (importante para que las annotations no se muevan)
fig = figure('Name','IL Normalizada – Base vs Sensor (glucosa)', ...
             'Color','w', ...
             'WindowState','maximized');

cols = lines(numel(concs_mM));

% --- Normalización respecto al valor máximo de IL_base ---
IL_base_norm = IL_base ./ max(IL_base);
IL_sense_norm_all = IL_sense_all ./ max(IL_sense_all,[],1);  % cada curva

% --- Eje principal --------------------------------------------------------
axMain = axes(fig); hold(axMain,'on'); grid(axMain,'on'); box(axMain,'on');

fMHz = f*1e-6; % frecuencia en MHz

plot(axMain, fMHz, IL_base_norm, 'k-', 'LineWidth',2);

for ic = 1:numel(concs_mM)
    plot(axMain, fMHz, IL_sense_norm_all(:,ic), '-', ...
        'LineWidth',2, 'Color', cols(ic,:));
end

xlabel(axMain,'f [MHz]','FontSize',25);
ylabel(axMain,'IL normalizado (adimensional)','FontSize',25);
title(axMain, sprintf('Insertion Loss Normalizado (η=%.2f, Randles-Cdl, %s, Rs_{elec}=%d Ω), Δf=%6.1f kHz', ...
    eta_sens, tern(USE_WARBURG,'con Warburg','sin Warburg'), Rs_elec, (f_min_IL(end)-f_at_minIL_base)/1e3 ), ...
    'FontSize',25);

leg = [{'Base (sin sensor)'} arrayfun(@(x) sprintf('%.0f mg/dL (%.1f mM)', ...
    concs_mg_dL(x), concs_mM(x)), 1:numel(concs_mM), 'uni', false)];
legend(axMain, leg, 'Location','best');

xlim(axMain,[min(fMHz) max(fMHz)]);    % opcional
ylim(axMain,[0 1.2]);                  % opcional

%% ==================== INSET DE ZOOM + RECUADRO ===========================
% Rango de zoom (en MHz)
xzoom = [3.15 3.40];

% Máscara en ese rango
mask = (fMHz >= min(xzoom)) & (fMHz <= max(xzoom));
if ~any(mask)
    warning('xzoom está fuera del rango de f.');
    return
end

% Datos en el rango para definir y-zoom
y_all_zoom = [IL_base_norm(mask); IL_sense_norm_all(mask,:)];
vals = y_all_zoom(:);
vals = vals(~isnan(vals) & ~isinf(vals));

ymin = min(vals);
ymax = max(vals) + 0.02;
if ymax <= ymin
    ymax = ymin + 0.05;
end
yzoom = [ymin ymax];

% --- Rectángulo en el eje principal --------------------------------------
rectPos = [min(xzoom), ymin, diff(xzoom), ymax-ymin];
rectangle(axMain,'Position',rectPos,'EdgeColor','k','LineWidth',1.2);

% --- Eje inset (a la derecha) --------------------------------------------
insetAx = axes('Parent',fig, ...
               'Position',[0.65 0.55 0.22 0.30]);  % [left bottom width height]
hold(insetAx,'on'); grid(insetAx,'on'); box(insetAx,'on');

plot(insetAx, fMHz, IL_base_norm, 'k-', 'LineWidth',1.4);
for ic = 1:numel(concs_mM)
    plot(insetAx, fMHz, IL_sense_norm_all(:,ic), '-', ...
        'LineWidth',1.4, 'Color', cols(ic,:));
end

xlim(insetAx,xzoom);
ylim(insetAx,yzoom);
set(insetAx,'FontSize',9);
title(insetAx,'Zoom');

%% --------- Líneas conectoras rectángulo <-> inset ------------------------
set(fig,'Units','normalized');
set(axMain,'Units','normalized');
set(insetAx,'Units','normalized');

ds2nfu = @(axh,xp,yp) deal( ...
    axh.Position(1) + (xp-axh.XLim(1))/diff(axh.XLim)*axh.Position(3), ...
    axh.Position(2) + (yp-axh.YLim(1))/diff(axh.YLim)*axh.Position(4) );

xRectRight = rectPos(1) + rectPos(3);
yRectTop   = rectPos(2) + rectPos(4);
yRectBot   = rectPos(2);

[xRT, yRT] = ds2nfu(axMain, xRectRight, yRectTop);  % sup. derecha
[xRB, yRB] = ds2nfu(axMain, xRectRight, yRectBot);  % inf. derecha

xInsetLeft = insetAx.Position(1);
yInsetTop  = insetAx.Position(2) + insetAx.Position(4);
yInsetBot  = insetAx.Position(2);

annotation(fig,'line',[xRT xInsetLeft],[yRT yInsetTop],'Color','k');
annotation(fig,'line',[xRB xInsetLeft],[yRB yInsetBot],'Color','k');

%% ==================== GRÁFICAS DE FASE (S21) ============================
% 3) Fase absoluta de S21: base vs. sensor
figure('Name','Fase(S21) – Base vs Sensor','Color','w');
plot(f*1e-6, phi_base_deg, 'k-', 'LineWidth',1.5); hold on; grid on;
cols = lines(numel(concs_mM));
for ic = 1:numel(concs_mM)
    plot(f*1e-6, S21_phase_deg(:,ic), '-', 'LineWidth',1.4, 'Color', cols(ic,:));
end
xlabel('f [MHz]');
ylabel('Fase(S_{21}) [grados]');
title(sprintf('Fase de S_{21} (η=%.2f, Randles-Cdl, %s, Rs_{elec}=%d Ω)', ...
    eta_sens, tern(USE_WARBURG,'con Warburg','sin Warburg'), Rs_elec));
leg = [{'Base (sin sensor)'} arrayfun(@(x) sprintf('%.0f mg/dL (%.1f mM)', ...
    concs_mg_dL(x), concs_mM(x)), 1:numel(concs_mM), 'uni', false)];
legend(leg, 'Location','best');

% 4) Fase diferencial: sensor – base
figure('Name','ΔFase(S21) – Sensor vs Base','Color','w');
for ic = 1:numel(concs_mM)
    dphi = S21_phase_deg(:,ic) - phi_base_deg;
    dphi = mod(dphi+180, 360) - 180;
    plot(f*1e-6, dphi, '-', 'LineWidth',1.4, 'Color', cols(ic,:)); hold on;
end
grid on;
xlabel('f [MHz]');
ylabel('ΔFase(S_{21}) [grados]');
title('Diferencia de fase (sensor – base)');
legend(arrayfun(@(x) sprintf('%.0f mg/dL (%.1f mM)', ...
    concs_mg_dL(x), concs_mM(x)), 1:numel(concs_mM), 'uni', false), 'Location','best');

%% ==================== Chequeo de escalas @ f0 ==========================
[~,i0] = min(abs(f - f0));
Y_saw_f0   = Y_saw(i0);

c_min_mM  = concs_mM(1);
c_max_mM  = concs_mM(end);
Cmin_cm3  = mM_to_mol_cm3(c_min_mM);
Cmax_cm3  = mM_to_mol_cm3(c_max_mM);

Ydl_all = jw .* Cdl;

% --- Ys para conc. mínima ---
Rct_min = Rct_of_c(c_min_mM);
if USE_WARBURG && (Cmin_cm3 > 0)
    sigma_min = (R*T)/( (n^2)*F^2*A_cm2*Cmin_cm3*sqrt(2*D) );
    Zw_min    = sigma_min ./ sqrt(jw);
    Yw_min    = 1 ./ Zw_min;
else
    Yw_min    = zeros(size(jw));
end
Ys_min_vec = Ydl_all + (1./Rct_min) + Yw_min;
if INCLUDE_Rs_SERIE
    Ys_min_vec = 1 ./ ( Rs_elec + 1./Ys_min_vec );
end

% --- Ys para conc. máxima ---
Rct_max = Rct_of_c(c_max_mM);
if USE_WARBURG && (Cmax_cm3 > 0)
    sigma_max = (R*T)/( (n^2)*F^2*A_cm2*Cmax_cm3*sqrt(2*D) );
    Zw_max    = sigma_max ./ sqrt(jw);
    Yw_max    = 1 ./ Zw_max;
else
    Yw_max    = zeros(size(jw));
end
Ys_max_vec = Ydl_all + (1./Rct_max) + Yw_max;
if INCLUDE_Rs_SERIE
    Ys_max_vec = 1 ./ ( Rs_elec + 1./Ys_max_vec );
end

Ys_min_f0  = Ys_min_vec(i0);
Ys_max_f0  = Ys_max_vec(i0);

fprintf('\n--- Chequeo de escalas a f0 ---\n');
fprintf('|Y_saw(f0)|   = %.3e S  (|Ga|=%.3e, |jωCt+Ba|=%.3e)\n', ...
    abs(Y_saw_f0), Gaf(i0), abs(1j*(2*pi*f(i0)*Ct)+Baf(i0)) );
fprintf('eta*|Y_sens(f0)| (min conc) = %.3e S\n', eta_sens*abs(Ys_min_f0));
fprintf('eta*|Y_sens(f0)| (max conc) = %.3e S\n', eta_sens*abs(Ys_max_f0));

%% ==================== CONSOLA: RESUMEN =================================
fprintf('\n<<< Ventana sensora (Rct físico: BV) >>>\n');
fprintf('eta_sens (Ws/Wa)          = %.2f (Ws = %.1f µm)\n', eta_sens, eta_sens*Wa_um);
fprintf('Modelo EC: Randles (Cdl ideal) | Difusión: %s\n', tern(USE_WARBURG,'Warburg semi-inf','No'));
fprintf('Cdl                       = %.3g F\n', Cdl);
fprintf('Params Rct: k0=%.2e cm/s, A=%.3f cm^2, n=%d, T=%g K\n', k0_cms, A_cm2, n, T);
for ic = 1:numel(concs_mM)
    C_mM   = concs_mM(ic);
    C_cm3  = mM_to_mol_cm3(C_mM);
    Rct_here = Rct_of_c(C_mM);

    if USE_WARBURG && (C_cm3 > 0)
        sigma_here = (R*T)/( (n^2)*F^2*A_cm2*C_cm3*sqrt(2*D) );
    else
        sigma_here = 0;
    end

    fprintf('  %3.0f mg/dL (=%4.1f mM): C=%.2e mol/cm^3  ->  Rct = %8.2f Ω, sigma=%.3g Ω·s^{-1/2} | Min IL=%.2f dB @ %.6f Hz | Δf=%6.1f kHz\n', ...
        concs_mg_dL(ic), concs_mM(ic), mM_to_mol_cm3(concs_mM(ic)), ...
        Rct_here, sigma_here, IL_min(ic), f_min_IL(ic), (f_min_IL(ic)-f_at_minIL_base)/1e3 );
end
fprintf('-----------------------------------------------------------------\n');

%% ====================== Helper functions ==============================
function val = x_of(f,f0,Np)
    X = pi*Np*((f - f0)./f0);
    X(abs(X) < 1e-18) = 1e-12;
    val = X;
end

function H = H_resp(f,f0,Np,Cs,k)
    X    = x_of(f,f0,Np);
    sinc_like = sin(X)./X;
    H = 2*k*sqrt(Cs*f0)*Np .* sinc_like .* exp(-1j*pi*f*Np./(2*f0)); % (1)
end

function G = Ga(f,f0,Np,Cs,k,Wa)
    X    = x_of(f,f0,Np);
    G = 8*(k^2)*Cs*f0*(Np^2)*Wa .* abs(sin(X)./X).^2;               % (4)
end

function B = Ba(f,f0,Np,Cs,k,Wa)
    X  = x_of(f,f0,Np);
    g0 = Ga(f0,f0,Np,Cs,k,Wa);
    B  = g0 .* ( sin(2*X) - 2*X ) ./ ( 2*(X.^2) );                  % (5)
end

function Wa = Wa_optimized(f0,Np,Cs,k,Rin)
    num = 4*Np*(k^2);
    den = (4*Np*(k^2))^2 + (pi^2);
    Wa  = (1/(2*Rin*Np*f0*Cs)) * (num/den);                         % (10)
end

function s = tern(cond,a,b)
    if cond
        s = a;
    else
        s = b;
    end
end

