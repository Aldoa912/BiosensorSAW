%% Línea de retardo SAW (modelo de respuesta al impulso) – Versión MATLAB
% Implementa H(f), G_a(f), B_a(f), Y(f), Z(f), IL(f) y optimización de Wa.

clear; clc;

%% -------- Parámetros (EDITAR AQUÍ) ---------------------------------------
% Ejemplo genérico (puedes sustituir por tu material).
wavelenght = 1e-3;                    % longitud de onda [m]
v_saw   = 3290;                       % velocidad acústica [m/s]
f0      = v_saw / wavelenght;         % frecuencia sincrónica [Hz]
NBW     = 5.5e5;                      % ancho de banda nulo [Hz]
k       = sqrt(0.0075);               % coeficiente de acoplamiento piezoeléctrico [-]
Rg      = 50;                         % resistencia de fuente/carga [Ohm]
Rin     = 50;                         % resistencia de entrada deseada para acoplamiento [Ohm]

%% --------- Cs a partir de constantes dieléctricas ------------------------
% Los valores de tabla son PERMITIVIDADES RELATIVAS (adimensionales)
eps0  = 8.854187817e-12;              % permitividad del vacío [F/m]

% <<< PON AQUÍ LOS VALORES DE TU TABLA >>>
eps11_rel = 53.6;                     % (T11)  -> relativa
eps33_rel = 43.4;                     % (T33)  -> relativa
eps13_rel = 0.0;                      % si no está reportado, usa 0.0

% Cálculo de epsilon_r efectivo
eps_r = sqrt( (eps11_rel/eps33_rel) - (eps13_rel^2)/(eps33_rel^2) );

% Cálculo de Cs con control de unidades
% Nota: “20” suele estar en pF/cm en literatura SAW; ajusta si tu fuente usa otra base.
Cs_base  = 20;                        % valor base en [pF/cm]
denom    = 0.5 * ( 1 + eps_r * (eps33_rel) );   % (ε33/ε0) = ε33_rel

Cs_pF_per_cm_calc = Cs_base * denom;           % resultado en pF/cm
% Conversión a SI [F/m]:
unit_conv = (1e-12) / (1e-2);                   % (pF -> F) / (cm -> m)
Cs = Cs_pF_per_cm_calc * unit_conv;             % [F/m]

% Mostrar Cs calculado
fprintf('Cs (calculado) = %.6f pF/cm  ->  %.3e F/m\n', ...
        Cs_pF_per_cm_calc, Cs);

%% -------- Dimensionamiento del IDT ---------------------------------------
lambda = v_saw / f0;                            % longitud de onda [m]
Np     = round((f0*2)/(NBW));                   % número de pares de dedos
% Np = 2;                                       % (opcional: forzar un valor)
Wa     = Wa_optimized(f0, Np, Cs, k, Rin);      % apertura [m]

% Capacitancia estática total (un solo IDT)
Ct     = Np * Wa * Cs;                          % [F]

% Barrido de frecuencia alrededor de f0 (~0.5 f0 .. 1.5 f0)
fmin = 0.5*f0;  fmax = 1.5*f0;  Npts = 8001;
f = linspace(fmin, fmax, Npts).';

%% -------- Cálculos principales -------------------------------------------
% Respuesta en frecuencia (dos IDTs idénticos): H_T(f) ~ [H(f)]^2 -> dB normalizados
Hf    = H_resp(f,  f0, Np, Cs, k);
Hf0   = H_resp(f0, f0, Np, Cs, k);
FRdB  = 20*log10( abs((Hf.^2) ./ (Hf0.^2)) );

% Conductancia de radiación G_a(f) y susceptancia acústica B_a(f)
Gaf   = Ga(f, f0, Np, Cs, k, Wa);
Baf   = Ba(f, f0, Np, Cs, k, Wa);
Bnf   = Baf / Gaf;                               

% Admitancia e impedancia de entrada del IDT
Yf    = Gaf + 1j*(2*pi*f*Ct + Baf);
Zf    = 1 ./ Yf;
Zre   = real(Zf);  
Zim   = imag(Zf);

% Pérdida de inserción (mínimo esperado cerca de f0)
ILdB  = -10*log10( (2*Gaf*Rg) ./ ((1 + Gaf*Rg).^2 + (Rg*(2*pi*f.*Ct + Baf)).^2) );

% Magnitudes clave para comprobación
lambda_um = lambda*1e6;
Wa_um     = Wa*1e6;
Ct_pF     = Ct*1e12;
[minIL, minIdx] = min(ILdB);
f_at_minIL = f(minIdx);

%% -------- Resumen en consola ---------------------------------------------
fprintf('\n<<< Resumen SAW Delay Line >>>\n');
fprintf('Longitud de onda          = %.6e m (%.3f µm)\n', lambda, lambda_um);
fprintf('Frecuencia sincrónica     = %.6e Hz\n', f0);
fprintf('Velocidad acústica        = %.3f m/s\n', v_saw);
fprintf('Pares de dedos (Np)       = %d\n', Np);
fprintf('Apertura (Wa)             = %.6e m (%.1f µm)\n', Wa, Wa_um);
fprintf('Capacitancia total (IDT)  = %.6e F (%.3f pF)\n', Ct, Ct_pF);
fprintf('Mínima IL                 = %.3f dB @ %.6e Hz\n', minIL, f_at_minIL);

%% -------- Gráficas rápidas (opcional) ------------------------------------
figure; plot(f*1e-6, FRdB, 'LineWidth',1.7);
xlabel('Frecuencia [MHz]','FontSize',25); ylabel('|\it{H_T}(f)| normalizada [dB]','FontSize',25); grid on;
title('Respuesta en frecuencia (dos IDTs idénticos)','FontSize',25);

figure; plot(f*1e-6, Gaf/max(Gaf), 'LineWidth',1.2); grid on;
xlabel('Frecuencia [MHz]'); ylabel('G_a(f) / G_a(f_0)'); 
title('Conductancia de radiación normalizada');

figure; plot(f*1e-6, Baf/max(Gaf), 'LineWidth',1.2); grid on;
xlabel('Frecuencia [MHz]'); ylabel('B_a(f) [S]'); 
title('Susceptancia acústica');

figure; plot(f*1e-6, ILdB, 'LineWidth',1.2); grid on;
xlabel('Frecuencia [MHz]'); ylabel('Insertion Loss [dB]'); 
title('Pérdida de inserción vs frecuencia');

figure; plot(f*1e-6, Zre, f*1e-6, Zim, 'LineWidth',1.2); grid on;
xlabel('Frecuencia [MHz]'); ylabel('Z(f) [\Omega]'); legend('Re\{Z\}','Im\{Z\}');
title('Impedancia de entrada (IDT único)');

%% ====================== Funciones auxiliares =============================

function val = x_of(f, f0, Np)
% Variable auxiliar X = pi*Np*((f-f0)/f0) con manejo seguro cerca de f=f0
    X = pi*Np*((f - f0)./f0);
    % Evitar división por cero cerca de X=0
    X(abs(X) < 1e-18) = 1e-12;
    val = X;
end

function H = H_resp(f, f0, Np, Cs, k)
% Respuesta en frecuencia de un IDT (compleja); para dos IDTs se usa H^2
% H(f) = 2*k*sqrt(Cs*f0)*Np * sinc_like * exp(-j*pi*f*Np/(2*f0))
    X    = x_of(f, f0, Np);
    sinc_like = sin(X)./X;
    H = 2*k*sqrt(Cs*f0)*Np .* sinc_like .* exp(-1j*pi*f*Np./(2*f0));
end

function G = Ga(f, f0, Np, Cs, k, Wa)
% Conductancia de radiación
    X = x_of(f, f0, Np);
    G = 8*(k^2)*Cs*f0*(Np^2)*Wa .* abs(sin(X)./X).^2;
end

function B = Ba(f, f0, Np, Cs, k, Wa)
% Susceptancia acústica
    X  = x_of(f, f0, Np);
    g0 = Ga(f0, f0, Np, Cs, k, Wa);
    B  = g0 .* ( sin(2*X) - 2*X ) ./ ( 2*(X.^2) );
end

function Wa = Wa_optimized(f0, Np, Cs, k, Rin)
% Apertura que aproxima el acoplamiento deseado a Rin
% Wa = [ 1/(2 Rin Np f0 Cs) ] * [ 4 Np k^2 / ( (4 Np k^2)^2 + pi^2 ) ]
    num = 4*Np*(k^2);
    den = (4*Np*(k^2))^2 + (pi^2);
    Wa  = (1/(2*Rin*Np*f0*Cs)) * (num/den);
end
