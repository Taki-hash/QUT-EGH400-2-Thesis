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

% Simulation
t_end = 0.5; % need to be large (to keep most of analysis in steady state)
dt = 1e-6;
Fs = 1 / dt;
num_harmonics = 5;
f_pwm = 40e3; % increasing can help reduce THD?
t_vect = 0:dt:t_end;

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
Vout_max = 15e-3;

%% Input Injection (Dual Input Model) 
% Define B and E matrices
E_matrix = [1/L; 0];

% Create dual-input system
sys_IN = ss(A_matrix, [B_matrix, E_matrix], C_matrix, [0 0]);

% Simulate
N_10pct = round(0.04 * length(t_vect));
t_vect_10 = t_vect(1:N_10pct);u_10 = u(1:N_10pct);v_tilde_10 = v_tilde(1:N_10pct);
[y_IN, tout_IN] = lsim(sys_IN, [u_10(:), v_tilde_10(:)], t_vect_10);

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

figure;sgtitle('Linear Simulation - Input Injection','Interpreter','latex')

subplot(3,1,1);
semilogx(w_u, 20*log10(mag_u), 'k', 'LineWidth', 1.5); hold on;
semilogx(w_u, 20*log10(mag_vtilde_IN), 'b--', 'LineWidth', 1.5);
xline(cutoff_freq_IN, 'k--', 'LineWidth', 1.2);
text(cutoff_freq_IN, cutoff_mag_db_IN - 5,sprintf('f_c = %.1f Hz', f_cutoff_IN),'Color', 'red', 'HorizontalAlignment', 'right');
xlabel('Frequency (rad/s)', 'Interpreter', 'latex');
ylabel('Magnitude (dB)', 'Interpreter', 'latex');
title('Overlay: Bode Magnitude of $H_u(s)$ vs. $H_{\tilde{v}}(s)$', 'Interpreter', 'latex');
legend({'$H_u(s)$: Control Input', '$H_{\tilde{v}}(s)$: Input Injection'}, 'Interpreter', 'latex');
grid on;

subplot(3,1,2);
semilogx(w_u, phase_u, 'k', 'LineWidth', 1.5); hold on;
semilogx(w_u, phase_vtilde_IN, 'b--', 'LineWidth', 1.5);
xlabel('Frequency (rad/s)', 'Interpreter', 'latex');
ylabel('Phase (deg)', 'Interpreter', 'latex');
title('Overlay: Phase of $H_u(s)$ vs. $H_{\tilde{v}}(s)$', 'Interpreter', 'latex');
legend({'Phase of $H_u(s)$', 'Phase of $H_{\tilde{v}}(s)$'}, 'Interpreter', 'latex');
grid on;

subplot(3,1,3);
plot(tout_IN, y_IN, 'LineWidth', 1.2); grid on;
xlabel('Time (s)', 'Interpreter', 'latex');
ylabel('Output Voltage', 'Interpreter', 'latex');
title(sprintf(['Time Domain Response to Input Injection:\n','$u(t) = D$, $\\tilde{v}(t) = %.2f\\sin(2\\pi \\times %.0f\\ t)$'], epsilon_V, f_tilde),'Interpreter', 'latex');

%% Input Injection - Amplitude & Frequency Limits

% Bandwidth
fprintf('3dB Bandwidth, Input Injection: %.2f Hz\n', f_cutoff_IN);  % in Hz

% Attenuation
max_attenuation_db_IN = 20*log10(max(mag_vtilde_IN));
fprintf('Maximum Gain: (%.2f dB)\n', max_attenuation_db_IN);  

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

figure;
semilogx(f_IN, epsilon_IN*100, 'b', 'LineWidth', 1.5); hold on;
yline(min_epsilon_IN*100, 'k--', 'Min', 'LabelHorizontalAlignment','left');
yline(max_epsilon_IN*100, 'r--', 'Max', 'LabelHorizontalAlignment','left');
xlabel('Frequency [Hz]', 'Interpreter', 'latex');
ylabel('Injection Amplitude $|\tilde{v}|$ [mV]', 'Interpreter', 'latex');
title('Input Injection: Amplitude Sweep Range', 'Interpreter', 'latex');
set(gca, 'FontSize', 10); % for better readability
grid on;


%% Input Injection Extension - Calculate \sum THD 

tot_THD_IN = zeros(length(f_IN));

for i = 1:length(f_IN)

    v_tilde = epsilon_IN(i) * sin(2*pi*f_IN(i)*t_vect);
    [~, x_in, ~, fn, ~] = simulateNonlinearBuck(Vpv, R, L, C, D, f_pwm, t_end, dt, v_tilde, 1, false, false);
    thd_sf_in = computeTHD(x_in(:,2), Fs, f_IN(i), num_harmonics);
    thd_sw_in = computeTHD(x_in(:,2), Fs, fn, num_harmonics);
    tot_THD_IN(i) = thd_sf_in + thd_sw_in;

end

figure;
plot3(f_IN, epsilon_IN, tot_THD_IN, 'b.-', 'LineWidth', 1.5, 'MarkerSize', 12);
xlabel('Frequency [Hz]', 'Interpreter', 'latex');
ylabel('Injection Amplitude $|\tilde{v}|$ [V]', 'Interpreter', 'latex');
zlabel('$\Sigma~\mathrm{THD}$ [\%]', 'Interpreter', 'latex');
title('Total Harmonic Distortion vs. Input Injection', 'Interpreter', 'latex');
set(gca, 'FontSize', 10);
zlim([0 5]);
grid on;
view(20, 25);

