function [t_vect, x_log, u_log, fn, f_ripple_sim] = simulateNonlinearBuck(Vpv, R, L, C, D, f_pwm, t_end, dt, v_tilde, injection, plot_flag, print_flag)

% Time vector
t_vect = 0:dt:t_end;
n_steps = length(t_vect);
Fs = 1 / dt;

% Preallocate logs
x_log = zeros(n_steps, 2);  % [i_L, v_C]
u_log = zeros(n_steps, 1);  % switching signal
carrier_log = zeros(n_steps, 1);
duty_log = zeros(n_steps, 1);

% Initial state
x = [0; 0];

% PWM carrier
carrier_period = 1 / f_pwm;

for k = 1:n_steps
    t = t_vect(k);
    i_L = x(1);
    v_C = x(2);

    % Sawtooth carrier (0 to 1)
    carrier = mod(t, carrier_period) / carrier_period;
    carrier_log(k) = carrier;

    % Duty command: perturbed only for PWM injection
    if injection == 3
        d_control = D + v_tilde(k);
    else
        d_control = D;
    end

    d_control = max(min(d_control, 1), 0);
    duty_log(k) = d_control;

    % PWM switching logic
    switch_on = d_control >= carrier;
    u = Vpv * switch_on;
    u_log(k) = u;

    % Dynamics (based on injection location)
    if injection == 1      % Input injection
        diL_dt = (u - v_C) / L + (1/L) * v_tilde(k);
        dvC_dt = (i_L - v_C / R) / C;

    elseif injection == 2  % Output injection
        diL_dt = (u - v_C) / L;
        dvC_dt = (i_L - v_C / R) / C + (1/(R*C)) * v_tilde(k);

    elseif injection == 3  % PWM injection (already applied in duty)
        diL_dt = (u - v_C) / L;
        dvC_dt = (i_L - v_C / R) / C;

    elseif injection == 0  % No injection
        diL_dt = (u - v_C) / L;
        dvC_dt = (i_L - v_C / R) / C;

    else
        error('Invalid injection mode. Use 0 (none), 1 (IN), 2 (OUT), or 3 (PWM).');
    end

    % Euler integration
    x(1) = x(1) + diL_dt * dt;
    x(2) = x(2) + dvC_dt * dt;

    % Log state
    x_log(k,:) = x;
end

% === Natural frequency of LC filter ===
w_n = 1 / sqrt(L*C);      % rad/s
fn = w_n / (2*pi);        % Hz

% === Output ripple frequency via FFT ===
v_C = x_log(:,2);
V_fft = abs(fft(v_C));
f_axis = (0:length(V_fft)-1) * Fs / length(V_fft);
[~, idx] = max(V_fft(2:end));  % skip DC
f_ripple_sim = f_axis(idx + 1);

% === Optional Plot ===
if plot_flag
    figure;
    subplot(3,1,1)
    plot(t_vect, x_log(:,1), 'b')
    ylabel('i_L [A]')
    title('Inductor Current')
    grid on

    subplot(3,1,2)
    plot(t_vect, x_log(:,2), 'r')
    ylabel('v_C [V]')
    title('Capacitor Voltage (Output)')
    grid on

    subplot(3,1,3)
    plot(t_vect, u_log, 'k')
    ylabel('u(t) [V]')
    xlabel('Time [s]')
    title('Switching Signal (PWM Output)')
    grid on
end

% === Optional Print ===
if print_flag
    fprintf('Natural frequency (LC): %.2f Hz\n', fn);
    fprintf('Simulated ripple frequency: %.2f Hz\n', f_ripple_sim);
end

end
