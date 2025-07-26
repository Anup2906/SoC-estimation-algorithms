clc; clear; close all;

%% Load Discharge Data
load("03-09-17_17.59 3349_Dis1C_1.mat")
time = meas.Time;
Vt = meas.Voltage;
I = -meas.Current;  % Invert if needed
nominal_capacity = 2.9;  % Ah

%% Compute SOC
remaining_capacity = meas.Ah;
remaining_capacity = -remaining_capacity(end) + remaining_capacity;
actual_soc = remaining_capacity / nominal_capacity;

%% Load OCV-SOC Curve from HPPC Data
load("03-11-17_08.47 25degC_5Pulse_HPPC_Pan18650PF.mat");
rest_idx = abs(meas.Current) < 0.01;
soc_all = (meas.Ah(rest_idx) - min(meas.Ah)) / nominal_capacity;
ocv_all = meas.Voltage(rest_idx);

[soc_unique, ia] = unique(soc_all, 'stable');
ocv_unique = ocv_all(ia);

p = polyfit(soc_unique, ocv_unique, 5);
OCV = @(soc) polyval(p, soc);

fprintf('OCV(SOC) = ');
deg = length(p) - 1;
for i = 1:length(p)
    coeff = p(i);
    pow = deg - (i - 1);
    if coeff >= 0 && i > 1, fprintf(' + ');
    elseif coeff < 0, fprintf(' - '); coeff = abs(coeff); end
    if pow == 0
        fprintf('%.10f', coeff);
    elseif pow == 1
        fprintf('%.10f*SOC', coeff);
    else
        fprintf('%.10f*SOC^%d', coeff, pow);
    end
end
fprintf('\n');

%% ECM Voltage Residual
VOC = OCV(actual_soc);
VRRC = VOC - Vt;

%% Time Step
%time_step = mean(diff(time));
time_step = 10;

%% Initial ECM Parameters (Good Estimates)
Rs = 0.1;     % within limit
Rp = 0.1;
Cp = 1500;

a1 = 1 - time_step / (Cp * Rp);
b0 = Rs;
b1 = -Rs + (time_step / Cp) + (time_step * Rs) / (Cp * Rp);

Rs(1) = b0; Rp(1) = Rp; Cp(1) = Cp;
B0(1) = b0; B1(1) = b1; A1(1) = a1;

%% Normalize Inputs for Stability
I_max = max(abs(I));
VRRC_max = max(abs(VRRC));
I = I / I_max;
VRRC = VRRC / VRRC_max;

%% RLS-VDF Initialization
lambda(1) = 0.95;
lambdamin = 0.90;
thetai = [a1; b0; b1];
Pi = 1e4 * eye(3);
epsilon = 1e-6;
Alpha = 0.98;
sigma(1) = 1e-3;
L(1) = 0;
y = VRRC;

%% RLS-VDF Loop
for k = 2:length(time)
    phie = [VRRC(k-1); I(k); I(k-1)];
    Pf = Pi - (Pi * phie * phie' * Pi) / (1 + phie' * Pi * phie);
    Pi = Pf;
    thetaf = thetai + Pf * phie * (y(k) - phie' * thetai);
    thetai = thetaf;

    A1(k) = thetaf(1);
    B0(k) = thetaf(2);
    B1(k) = thetaf(3);

    phietheta = phie' * thetaf;
    e(k) = y(k) - phietheta;

    sigma(k) = Alpha * sigma(k-1)^2 + (1 - Alpha) * e(k)^2;
    L(k) = (e(k)^2) / (sigma(k)^2);
    lambda(k) = lambdamin + (1 - lambdamin) * 0.5^L(k);

    if abs(phie) < epsilon
        Pi = Pi;
    else
        Pi = Pi + ((1 - lambda(k)) / lambda(k)) * ((phie * phie') / (phie' * inv(Pi) * phie));
    end

    Rs(k) = max(1e-4, min(B0(k), 0.02));
    Rp(k) = max(1e-4, min((B0(k) + A1(k) * B1(k)) / (1 - A1(k)), 0.02));
    Cp(k) = max(100, min(time_step / (B1(k) + A1(k) * B0(k)), 2000));
end

%% Un-normalize Voltage for Tracking
VOC = VOC;  % already in original scale
I = I * I_max;       % undo normalization
VRRC = VRRC * VRRC_max;

%% Time Vector
k_time = (1:length(time)) * time_step;

%% Final Results
fprintf('\n--- Final Estimated Parameters ---\n');
fprintf('Final Estimated Rs = %.6f Ohm\n', Rs(end));
fprintf('Final Estimated Rp = %.6f Ohm\n', Rp(end));
fprintf('Final Estimated Cp = %.6f Farad\n', Cp(end));

%% Save parameters
save('estimated_params.mat', 'Rs', 'Rp', 'Cp');

%% Tracking Voltage Using Final RC Parameters with Real-Time Adaptive Bias Correction

N = 50;
Rp_const = mean(Rp(end-N:end));
Cp_const = mean(Cp(end-N:end));
Rs_const = mean(Rs(end-N:end));

a1_const = 1 - (time_step / (Rp_const * Cp_const));
b1_const = time_step / Cp_const;

V_RC = zeros(length(time), 1);
V_track_RC = zeros(length(time), 1);
b = zeros(length(time), 1);         % Adaptive bias term
gamma = 0.5;                        % Learning rate (tune if needed)

% Initial values
V_RC(1) = 0;
V_track_RC(1) = VOC(1) - I(1)*Rs_const - V_RC(1) + b(1);

for k = 2:length(time)
    V_RC(k) = a1_const * V_RC(k-1) + b1_const * I(k-1);
    V_track_RC(k) = VOC(k) - I(k)*Rs_const - V_RC(k) + b(k-1);
    
    % Update adaptive bias
    b(k) = b(k-1) + gamma * (Vt(k-1) - V_track_RC(k-1));
end

%% Plotting
figure;
plot(time, Vt, 'b', 'LineWidth', 1.5); hold on;
plot(time, V_track_RC, 'r--', 'LineWidth', 1.5);
legend('Measured V_t', 'Tracked V_{track}');
xlabel('Time (s)'); ylabel('Voltage (V)');
title('Measured vs Tracked Terminal Voltage (Constant RC Model)');
grid on;

figure;
subplot(3,1,1); plot(time, Rs); ylabel('R_s');
subplot(3,1,2); plot(time, Rp); ylabel('R_p');
subplot(3,1,3); plot(time, Cp); ylabel('C_p');
xlabel('Time'); sgtitle('Estimated Parameters');


%% Tracking Error
tracking_error = Vt - V_track_RC;
valid_idx = ~isnan(V_track_RC);

figure;
plot(time(valid_idx), abs(tracking_error(valid_idx)), 'm', 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('|Voltage Error| (V)');
title('Absolute Voltage Tracking Error');
grid on;

%% Error Stats
fprintf('\n--- Tracking Error Stats ---\n');
fprintf('Mean Absolute Error = %.6f V\n', mean(abs(tracking_error(valid_idx))));
fprintf('Max Absolute Error = %.6f V\n', max(abs(tracking_error(valid_idx))));
% Calculate percentage error point-wise, then take the mean
percentage_error = 100 * abs(tracking_error(valid_idx)) ./ abs(Vt(valid_idx));

figure;
plot(time(valid_idx), percentage_error, 'g', 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('Percentage Error (%)');
title('Voltage Tracking Percentage Error');
grid on;

% Error Stats
fprintf('\n--- Percentage Tracking Error Stats ---\n');
fprintf('Mean Percentage Error = %.6f %%\n', mean(percentage_error));
fprintf('Max Percentage Error = %.6f %%\n', max(percentage_error));


figure;
plot(time(valid_idx), percentage_error, 'g', 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('Percentage Error (%)');
title('Voltage Tracking Percentage Error');
grid on;

fprintf('\n--- Percentage Tracking Error Stats ---\n');
fprintf('Mean Percentage Error = %.6f %%\n', mean(abs(percentage_error)));
fprintf('Max Percentage Error = %.6f %%\n', max(abs(percentage_error)));

%% EKF for SoC Estimation
fprintf('\n--- Starting EKF for SoC Estimation ---\n');

% Battery model: Vt = OCV(SOC) - I*Rs - V_RC, where V_RC is modeled as an RC circuit

% Initialize EKF parameters
Q = 1e-6;   % Process noise covariance
R = 5e-4;   % Measurement noise covariance
x = [actual_soc(1); 0];  % Initial state: [SoC; V_RC]
P = eye(2);              % Initial covariance

Rs_ekf = Rs_const;
Rp_ekf = Rp_const;
Cp_ekf = Cp_const;
Voc_ekf = @(soc) polyval(p, min(max(soc, 0), 1));  % Ensure SOC in [0,1]

a1_ekf = 1 - (time_step / (Rp_ekf * Cp_ekf));
b1_ekf = time_step / Cp_ekf;

soc_est = zeros(length(time), 1);
soc_est(1) = x(1);

for k = 2:length(time)
    % Control input
    Ik = I(k-1);

    % Predict Step
    soc_pred = x(1) - (time_step / (3600 * nominal_capacity)) * Ik;
    vrc_pred = a1_ekf * x(2) + b1_ekf * Ik;
    x_pred = [soc_pred; vrc_pred];

    % Jacobian of state transition (F)
    F = [1, 0;
         0, a1_ekf];

    % Predict covariance
    P = F * P * F' + Q * eye(2);

    % Measurement prediction
    V_pred = Voc_ekf(x_pred(1)) - Ik * Rs_ekf - x_pred(2);

    % Jacobian of measurement model (H)
    dVoc_dSoC = polyval(polyder(p), x_pred(1));
    H = [dVoc_dSoC, -1];

    % Kalman Gain
    S = H * P * H' + R;
    K = P * H' / S;

    % Update step
    yk = Vt(k) - V_pred;
    x = x_pred + K * yk;

    % Update covariance
    P = (eye(2) - K * H) * P;

    % Clip SoC to [0,1]
    x(1) = max(0, min(1, x(1)));

    soc_est(k) = x(1);
end

%% Plot SoC Estimation
figure;
plot(time, actual_soc, 'b', 'LineWidth', 1.5); hold on;
plot(time, soc_est, 'r--', 'LineWidth', 1.5);
legend('True SoC', 'Estimated SoC (EKF)');
xlabel('Time (s)'); ylabel('SoC');
title('Extended Kalman Filter Based SoC Estimation');
grid on;

%% SoC Estimation Error Metrics
soc_error = actual_soc - soc_est;
RMSE_soc = sqrt(mean(soc_error.^2));
MAE_soc = mean(abs(soc_error));
MPAE_soc = mean(abs(soc_error))*100;
fprintf('\n--- EKF SoC Estimation Error Metrics ---\n');
fprintf('RMSE = %.6f\n', RMSE_soc);
fprintf('MAE  = %.6f\n', MAE_soc);
fprintf('MPAE  = %.6f\n', MPAE_soc);
