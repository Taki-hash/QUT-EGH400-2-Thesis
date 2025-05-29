clear;close all;clc

% Component Values
L = 50e-6; % [uF]
C = 50e-6; % [uH]
R = 5; % [Ohms]
Vpv = 48; % [V] 
D = 0.5;

% State-space matrices
A_matrix = [0, -1/L;1/C, -1/(R*C)];
B_matrix = [1/L;0];
C_matrix = [0 1];
D_matrix = [0];

% Simulation time
t_end = 30e-2; % need to be large (to keep most of analysis in steady state)
dt = 1e-6;
Fs = 1 / dt;
num_harmonics = 5;
f_pwm = 80e3; % increasing can help reduce THD?
t_vect = 0:dt:t_end;

% Time Vector
%t_vect = linspace(0,10e-3,10000); % [s]

% Sinusoidal Perturbation 
epsilon_V = 1; % [V]
f_tilde = 1e3; % [Hz]
v_tilde = epsilon_V * sin(2*pi*f_tilde*t_vect);

% Control Input Initialisation (used in all injection simualtions) 
u = ones(size(t_vect)) * D * Vpv; % constant duty cycle

% Transfer Function - Control Input
[num_u, den_u] = ss2tf(A_matrix, B_matrix, C_matrix, 0);
H_u = tf(num_u, den_u);

% Frequency Response
[mag_u, phase_u, w_u] = bode(H_u);
mag_u = squeeze(mag_u);
phase_u = squeeze(phase_u);

% Desired maximum output amplitude (e.g., 10 mV)
Vout_max = 10e-3;

%% Input Injection (Dual Input Model) 
% Define B and E matrices
E_matrix = [1/L; 0];

% Create dual-input system
sys_IN = ss(A_matrix, [B_matrix, E_matrix], C_matrix, [0 0]);

% Simulate
[y_IN, tout_IN] = lsim(sys_IN, [u(:), v_tilde(:)], t_vect);

% Extract small-signal TF
[num_IN, den_IN] = ss2tf(A_matrix, E_matrix, C_matrix, 0);
H_vtilde_IN = tf(num_IN, den_IN); % perturbation input TF

% Frequency Analysis - Small Signal
[mag_vtilde_IN, phase_vtilde_IN, ~] = bode(H_vtilde_IN, w_u);
mag_vtilde_IN = squeeze(mag_vtilde_IN);
phase_vtilde_IN = squeeze(phase_vtilde_IN);

% Cutoff frequency (3dB from DC)
cutoff_mag_db_IN = 20*log10(mag_vtilde_IN(1)) - 3;
idx_cutoff_IN = find(20*log10(mag_vtilde_IN) <= cutoff_mag_db_IN, 1, 'first');
cutoff_freq_IN = w_u(idx_cutoff_IN);
f_cutoff_IN = cutoff_freq_IN / (2*pi);

% Plots
figure;

subplot(3,1,1);
semilogx(w_u, 20*log10(mag_u), 'k', 'LineWidth', 1.5); hold on;
semilogx(w_u, 20*log10(mag_vtilde_IN), 'b--', 'LineWidth', 1.5);
xline(cutoff_freq_IN, 'k--', 'LineWidth', 1.2);
text(cutoff_freq_IN, cutoff_mag_db_IN - 5, sprintf('f_c = %.1f Hz', f_cutoff_IN),'Color', 'red', 'HorizontalAlignment', 'right');
xlabel('Frequency (rad/s)'), ylabel('Magnitude (dB)');
title('Overlay: Bode Magnitude of H_u(s) vs. H_vtilde(s)');
legend('H_u(s): Control Input', 'H_vtilde(s): Input Injection');
grid on;

subplot(3,1,2);
semilogx(w_u, phase_u, 'k', 'LineWidth', 1.5); hold on;
semilogx(w_u, phase_vtilde_IN, 'b--', 'LineWidth', 1.5);
xlabel('Frequency (rad/s)'), ylabel('Phase (deg)');
title('Overlay: Phase of H_u(s) vs. H_vtilde(s)');
legend('Phase of H_u(s)', 'Phase of H_vtilde(s)');
grid on;

subplot(3,1,3);
plot(tout_IN, y_IN, 'LineWidth', 1.2); grid on;
xlabel('Time (s)'), ylabel('Output Voltage');
title('Output Response to Input Injection: u(t) = D, ṽ(t) = εsin(ωt)');

%% Input Injection - Amplitude & Frequency Limits

% Bandwidth
fprintf('3dB Bandwidth, Input Injection: %.2f Hz\n', f_cutoff_IN);  % in Hz

% Attenuation
max_attenuation_db_IN = 20*log10(max(mag_vtilde_IN));
fprintf('Maximum Attenuation: (%.2f dB)\n', max_attenuation_db_IN);  

% Frequency idx 
f_idx = 1:idx_cutoff_IN;
f_IN = w_u(f_idx)/(2*pi);

% Amplitude's within bandwidth
epsilon_IN = Vout_max ./ mag_vtilde_IN(f_idx);

min_epsilon_IN = min(epsilon_IN);
max_epsilon_IN = max(epsilon_IN);

% Print results
fprintf('Input Injections:\n');
fprintf('Min amplitude Input: %.6f mV\n', min_epsilon_IN*100);
fprintf('Max amplitude Input: %.6f mV\n', max_epsilon_IN*100);
fprintf('-------------------------------------------\n');

% Plot
figure;
semilogx(f_IN, epsilon_IN*100, 'b', 'LineWidth', 1.5); hold on;
yline(min_epsilon_IN*100, 'k--', 'Min', 'LabelHorizontalAlignment','left');
yline(max_epsilon_IN*100, 'r--', 'Max', 'LabelHorizontalAlignment','left');
xlabel('Frequency [Hz]');
ylabel('Amplitude [mV]');
title('Input Injection, Amplitude Sweep Range');
grid on;

%% Output Injection (Dual Input Model) 
% Define E matrix for output injection
E_matrix = [0; 1/(C*R)];  % Perturbation input path (capacitor current)

% Create dual-input system
sys_OUT = ss(A_matrix, [B_matrix, E_matrix], C_matrix, [0 0]);

% Simulate
[y_OUT, tout_OUT] = lsim(sys_OUT, [u(:), v_tilde(:)], t_vect);

% Extract small-signal TF
[num_OUT, den_OUT] = ss2tf(A_matrix, E_matrix, C_matrix, 0);
H_vtilde_OUT = tf(num_OUT, den_OUT);  % perturbation input TF

% Frequency Analysis - Small Signal
[mag_vtilde_OUT, phase_vtilde_OUT, ~] = bode(H_vtilde_OUT, w_u);  % match freq vector
mag_vtilde_OUT = squeeze(mag_vtilde_OUT);
phase_vtilde_OUT = squeeze(phase_vtilde_OUT);

% Cutoff frequency (3dB from DC or peak)
cutoff_mag_db_OUT = 20 * log10(mag_vtilde_OUT);
[peak_mag_db_OUT, idx_peak_OUT] = max(cutoff_mag_db_OUT);
cutoff_level_db_OUT = peak_mag_db_OUT - 3;
idx_lower_OUT = find(cutoff_mag_db_OUT(1:idx_peak_OUT) <= cutoff_level_db_OUT, 1, 'last');
idx_upper_OUT = find(cutoff_mag_db_OUT(idx_peak_OUT:end) <= cutoff_level_db_OUT, 1, 'first') + idx_peak_OUT - 1;
f_lower_OUT = w_u(idx_lower_OUT) / (2*pi);
f_upper_OUT = w_u(idx_upper_OUT) / (2*pi);
bandwidth_OUT = f_upper_OUT - f_lower_OUT;

% Plots
figure;

subplot(3,1,1);
semilogx(w_u, 20*log10(mag_u), 'k', 'LineWidth', 1.5); hold on;
semilogx(w_u, cutoff_mag_db_OUT, 'b--', 'LineWidth', 1.5);
xline(w_u(idx_lower_OUT), 'k--', 'LineWidth', 1.2);
xline(w_u(idx_upper_OUT), 'k--', 'LineWidth', 1.2);
text(w_u(idx_lower_OUT), cutoff_level_db_OUT - 5,sprintf('f_L = %.1f Hz', f_lower_OUT),'Color', 'red', 'HorizontalAlignment', 'right');
text(w_u(idx_upper_OUT), cutoff_level_db_OUT - 5,sprintf('f_H = %.1f Hz', f_upper_OUT),'Color', 'red', 'HorizontalAlignment', 'left');
xlabel('Frequency (rad/s)'), ylabel('Magnitude (dB)');
title('Overlay: Bode Magnitude of H_u(s) vs. H_vtilde(s)');
legend('H_u(s): Control Input', 'H_vtilde(s): Input Injection');
grid on;

subplot(3,1,2);
semilogx(w_u, phase_u, 'k', 'LineWidth', 1.5); hold on;
semilogx(w_u, phase_vtilde_OUT, 'b--', 'LineWidth', 1.5);
xlabel('Frequency (rad/s)'), ylabel('Phase (deg)');
title('Overlay: Phase of H_u(s) vs. H_vtilde(s)');
legend('Phase of H_u(s)', 'Phase of H_vtilde(s)');
grid on;

subplot(3,1,3);
plot(tout_OUT, y_OUT, 'LineWidth', 1.2); grid on;
xlabel('Time (s)'), ylabel('Output Voltage');
title('Output Response to Output Injection: u(t) = D, ṽ(t) = εsin(ωt)');

%% Output Injection - Amplitude & Frequency Limits

% Bandwidth
fprintf('3dB Bandwidth, OUT Injection Upper: %.2f Hz\n', f_lower_OUT);  % in Hz
fprintf('3dB Bandwidth, OUT Injection Lower: %.2f Hz\n', f_upper_OUT);  % in Hz

% Attenuation
max_attenuation_db_OUT = 20*log10(max(mag_vtilde_OUT));
fprintf('Maximum Attenuation: (%.2f dB)\n', max_attenuation_db_OUT);  

% Frequency idx 
f_idx = idx_lower_OUT:idx_upper_OUT;
f_OUT = w_u(f_idx)/(2*pi);

% Amplitude's within bandwidth
epsilon_OUT = Vout_max ./ mag_vtilde_OUT(f_idx);

min_epsilon_OUT = min(epsilon_OUT);
max_epsilon_OUT = max(epsilon_OUT);

% Print results
fprintf('Input Injections:\n');
fprintf('Min amplitude Input: %.6f mV\n', min_epsilon_OUT*100);
fprintf('Max amplitude Input: %.6f mV\n', max_epsilon_OUT*100);
fprintf('-------------------------------------------\n');

% Plot
figure;
semilogx(f_OUT, epsilon_OUT*100, 'b', 'LineWidth', 1.5); hold on;
yline(min_epsilon_OUT*100, 'k--', 'Min', 'LabelHorizontalAlignment','left');
yline(max_epsilon_OUT*100, 'r--', 'Max', 'LabelHorizontalAlignment','left');
xlabel('Frequency [Hz]');
ylabel('Amplitude [mV]');
title('Input Injection, Amplitude Sweep Range');
grid on;

%% PWM Injection (Dual Input Model)
% Define E matrix for PWM comparator modulation
E_matrix = [1/L; 0];  % Perturbation input path (same as B)

% Create dual-input system
sys_PWM = ss(A_matrix, [B_matrix, E_matrix], C_matrix, [0 0]);

% Simulate
[y_PWM, tout_PWM] = lsim(sys_PWM, [u(:), v_tilde(:)], t_vect);

% Extract small-signal TF
[num_PWM, den_PWM] = ss2tf(A_matrix, E_matrix, C_matrix, 0);
H_vtilde_PWM = tf(num_PWM, den_PWM);

% Frequency Analysis - Small Signal
[mag_vtilde_PWM, phase_vtilde_PWM, ~] = bode(H_vtilde_PWM, w_u);
mag_vtilde_PWM = squeeze(mag_vtilde_PWM);
phase_vtilde_PWM = squeeze(phase_vtilde_PWM);

% Cutoff frequency (3dB from DC gain)
cutoff_mag_db_PWM = 20*log10(mag_vtilde_PWM(1)) - 3;
idx_cutoff_PWM = find(20*log10(mag_vtilde_PWM) <= cutoff_mag_db_PWM, 1, 'first');
cutoff_freq_PWM = w_u(idx_cutoff_PWM);
f_cutoff_PWM = cutoff_freq_PWM / (2*pi);  % Hz

% Plots
figure;

subplot(3,1,1);
semilogx(w_u, 20*log10(mag_u), 'k', 'LineWidth', 1.5); hold on;
semilogx(w_u, 20*log10(mag_vtilde_PWM), 'b--', 'LineWidth', 1.5);
xline(cutoff_freq_PWM, 'k--', 'LineWidth', 1.2);
text(cutoff_freq_PWM, cutoff_mag_db_PWM - 5, sprintf('f_c = %.1f Hz', f_cutoff_PWM),'Color', 'red', 'HorizontalAlignment', 'right');
xlabel('Frequency (rad/s)'), ylabel('Magnitude (dB)');
title('Overlay: Bode Magnitude of H_u(s) vs. H_vtilde(s)');
legend('H_u(s): Control Input', 'H_vtilde(s): Input Injection');
grid on;

subplot(3,1,2);
semilogx(w_u, phase_u, 'k', 'LineWidth', 1.5); hold on;
semilogx(w_u, phase_vtilde_PWM, 'b--', 'LineWidth', 1.5);
xlabel('Frequency (rad/s)'), ylabel('Phase (deg)');
title('Overlay: Phase of H_u(s) vs. H_vtilde(s)');
legend('Phase of H_u(s)', 'Phase of H_vtilde(s)');
grid on;

subplot(3,1,3);
plot(tout_PWM, y_PWM, 'LineWidth', 1.2); grid on;
xlabel('Time (s)'), ylabel('Output Voltage');
title('Output Response to PWM Injection: u(t) = D, ṽ(t) = \epsilon sin(\omega t)');

%% PWM Injection - Amplitude & Frequency Limits

% Bandwidth
fprintf('3dB Bandwidth, PWM Injection: %.2f Hz\n', cutoff_freq_PWM);  % in Hz

% Attenuation
max_attenuation_db_PWM = 20*log10(max(mag_vtilde_PWM));
fprintf('Maximum Attenuation: (%.2f dB)\n', max_attenuation_db_PWM);  

% Frequency idx 
f_idx = 1:idx_cutoff_PWM;
f_PWM = w_u(f_idx)/(2*pi);

% Amplitude's within bandwidth
epsilon_PWM = Vout_max ./ mag_vtilde_PWM(f_idx);

min_epsilon_PWM = min(epsilon_PWM);
max_epsilon_PWM = max(epsilon_PWM);

% Print results
fprintf('Input Injections:\n');
fprintf('Min amplitude Input: %.6f mV\n', min_epsilon_PWM*100);
fprintf('Max amplitude Input: %.6f mV\n', max_epsilon_PWM*100);
fprintf('-------------------------------------------\n');

% Plot
figure;
semilogx(f_PWM, epsilon_PWM*100, 'b', 'LineWidth', 1.5); hold on;
yline(min_epsilon_PWM*100, 'k--', 'Min', 'LabelHorizontalAlignment','left');
yline(max_epsilon_PWM*100, 'r--', 'Max', 'LabelHorizontalAlignment','left');
xlabel('Frequency [Hz]');
ylabel('Amplitude [mV]');
title('Input Injection, Amplitude Sweep Range');
grid on;
