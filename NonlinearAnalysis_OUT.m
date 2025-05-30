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
Vout_max = 10e-3;

%% Output Injection (Dual Input Model) 
% Define E matrix for output injection
E_matrix = [0; 1/(C*R)];  % Perturbation input path (capacitor current)

% Create dual-input system
sys_OUT = ss(A_matrix, [B_matrix, E_matrix], C_matrix, [0 0]);

% Simulate
N_10pct = round(0.04 * length(t_vect));
t_vect_10 = t_vect(1:N_10pct);u_10 = u(1:N_10pct);v_tilde_10 = v_tilde(1:N_10pct);
[y_OUT, tout_OUT] = lsim(sys_OUT, [u_10(:), v_tilde_10(:)], t_vect_10);

% Extract small-signal TF
[num_OUT, den_OUT] = ss2tf(A_matrix, E_matrix, C_matrix, 0);
H_vtilde_OUT = tf(num_OUT, den_OUT);  % perturbation input TF

% Frequency Analysis - Small Signal
[mag_vtilde_OUT, phase_vtilde_OUT, ~] = bode(H_vtilde_OUT, w_u);  % match freq vector
mag_vtilde_OUT = squeeze(mag_vtilde_OUT)*0.01; % scale by injection amplitde (only for output injection
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
figure;sgtitle('Linear Simulation - Output Injection','Interpreter','latex')

subplot(3,1,1);
semilogx(w_u, 20*log10(mag_u), 'k', 'LineWidth', 1.5); hold on;
semilogx(w_u, cutoff_mag_db_OUT, 'b--', 'LineWidth', 1.5);
xline(w_u(idx_lower_OUT), 'k--', 'LineWidth', 1.2);
xline(w_u(idx_upper_OUT), 'k--', 'LineWidth', 1.2);
text(w_u(idx_lower_OUT)-500, cutoff_level_db_OUT - 30,sprintf('f_L = %.1f Hz', f_lower_OUT),'Color', 'red', 'HorizontalAlignment', 'right');
text(w_u(idx_upper_OUT)+500, cutoff_level_db_OUT - 30,sprintf('f_H = %.1f Hz', f_upper_OUT),'Color', 'red', 'HorizontalAlignment', 'left');
xlabel('Frequency (rad/s)', 'Interpreter', 'latex');
ylabel('Magnitude (dB)', 'Interpreter', 'latex');
title('Overlay: Bode Magnitude of $H_u(s)$ vs. $H_{\tilde{v}}(s)$', 'Interpreter', 'latex');
legend({'$H_u(s)$: Control Input', '$H_{\tilde{v}}(s)$: Output Injection'}, 'Interpreter', 'latex');
grid on;

subplot(3,1,2);
semilogx(w_u, phase_u, 'k', 'LineWidth', 1.5); hold on;
semilogx(w_u, phase_vtilde_OUT, 'b--', 'LineWidth', 1.5);
xlabel('Frequency (rad/s)', 'Interpreter', 'latex');
ylabel('Phase (deg)', 'Interpreter', 'latex');
title('Overlay: Phase of $H_u(s)$ vs. $H_{\tilde{v}}(s)$', 'Interpreter', 'latex');
legend({'Phase of $H_u(s)$', 'Phase of $H_{\tilde{v}}(s)$'}, 'Interpreter', 'latex');
grid on;

subplot(3,1,3);
plot(tout_OUT, y_OUT, 'LineWidth', 1.2); grid on;
xlabel('Time (s)', 'Interpreter', 'latex');
ylabel('Output Voltage', 'Interpreter', 'latex');
title(sprintf(['Time Domain Response to Output Injection:\n','$u(t) = D$, $\\tilde{v}(t) = %.2f\\sin(2\\pi \\times %.0f\\ t)$'],epsilon_V, f_tilde),'Interpreter', 'latex');


%% Output Injection - Amplitude & Frequency Limits

% Bandwidth
fprintf('3dB Bandwidth, OUT Injection Upper: %.2f Hz\n', f_lower_OUT);  % in Hz
fprintf('3dB Bandwidth, OUT Injection Lower: %.2f Hz\n', f_upper_OUT);  % in Hz

% Frequency idx 
f_idx = idx_lower_OUT:idx_upper_OUT;
f_OUT = w_u(f_idx)/(2*pi);

% Attenuation
max_attenuation_db_OUT = 20*log10(min(mag_vtilde_OUT(f_idx))); 
fprintf('Maximum Gain: (%.2f dB)\n', max_attenuation_db_OUT);  

% Amplitude's within bandwidth
epsilon_OUT = Vout_max ./ mag_vtilde_OUT(f_idx);

min_epsilon_OUT = min(epsilon_OUT);
max_epsilon_OUT = max(epsilon_OUT);

% Print results
fprintf('Input Injections:\n');
fprintf('Min amplitude OUT: %.6f mV\n', min_epsilon_OUT*100);
fprintf('Max amplitude OUT: %.6f mV\n', max_epsilon_OUT*100);
fprintf('-------------------------------------------\n');

figure;
semilogx(f_OUT, epsilon_OUT*100, 'b', 'LineWidth', 1.5); hold on;
yline(min_epsilon_OUT*100, 'k--', 'Min', 'LabelHorizontalAlignment','left');
yline(max_epsilon_OUT*100, 'r--', 'Max', 'LabelHorizontalAlignment','left');
xlabel('Frequency [Hz]', 'Interpreter', 'latex');
ylabel('Injection Amplitude $|\tilde{v}|$ [mV]', 'Interpreter', 'latex');
title('Output Injection: Amplitude Sweep Range', 'Interpreter', 'latex');
set(gca, 'FontSize', 10);
grid on;


%% Output Injection Extension - Calculate \sum THD 

tot_THD_OUT = zeros(length(f_OUT));

for i = 1:length(f_OUT)

    v_tilde = epsilon_OUT(i) * sin(2*pi*f_OUT(i)*t_vect);
    [~, x_out, ~, fn, ~] = simulateNonlinearBuck(Vpv, R, L, C, D, f_pwm, t_end, dt, v_tilde, 2, false, false);
    thd_sf_out = computeTHD(x_out(:,2), Fs, f_OUT(i), num_harmonics);
    thd_sw_out = computeTHD(x_out(:,2), Fs, fn, num_harmonics);
    tot_THD_OUT(i) = thd_sf_out + thd_sw_out;

end

figure;
plot3(f_OUT, epsilon_OUT, tot_THD_OUT, 'b.-', 'LineWidth', 1.5, 'MarkerSize', 12);
xlabel('Frequency [Hz]', 'Interpreter', 'latex');
ylabel('Injection Amplitude $|\tilde{v}|$ [V]', 'Interpreter', 'latex');
zlabel('$\Sigma~\mathrm{THD}$ [\%]', 'Interpreter', 'latex');
title('Total Harmonic Distortion vs. Output Injection', 'Interpreter', 'latex');
set(gca, 'FontSize', 10);
grid on;
view(35, 30);
