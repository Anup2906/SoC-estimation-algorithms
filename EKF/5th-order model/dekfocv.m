clear; clc; close all;

%% Load Discharge Data
load("soc_estimatedresult.mat");
time = 0:10.01:3800;
current = results.current;
voltage = results.voltage;
nominal_capacity = 2.9;
soc_true = results.soc_est;
soc_true(1)= 0.9641;
dt = 10;
N = length(time);

%% Estimate OCV–SOC Polynomial from HPPC Data
load("03-11-17_08.47 25degC_5Pulse_HPPC_Pan18650PF.mat");

rest_idx = abs(meas.Current) < 0.01;
soc_all = (meas.Ah(rest_idx) - min(meas.Ah)) / nominal_capacity;
ocv_all = meas.Voltage(rest_idx);

[soc_unique, ia] = unique(soc_all, 'stable');
ocv_unique = ocv_all(ia);

p = polyfit(soc_unique, ocv_unique, 5);
dp = polyder(p);

%% Initialize Model
Q = 2.9 * 3600;
x = [soc_true(1); 0];
theta = [0.02; 0.02; 2000];

Px = diag([1e-5, 1e-5]);
Ptheta = diag([0.002, 0.005, 100]);
Qx = diag([1e-5, 1e-5]);
Qtheta = diag([1e-5, 1e-5, 1]);
R = 1e-3;

x_store = zeros(2, N);
theta_store = zeros(3, N);
V_store = zeros(N, 1);

%% DEKF Loop
for k = 2:N
    I = current(k-1);
    V = voltage(k);

    % Prediction
    SoC_pred = x(1) - dt * I / Q;
    V_RC_pred = x(2) + dt * (-x(2)/(theta(2)*theta(3)) + I/theta(3));
    x_pred = [SoC_pred; V_RC_pred];

    A = [1, 0; 0, 1 - dt/(theta(2)*theta(3))];

    W = zeros(2,3);
    W(2,2) = dt * x(2) / (theta(2)^2 * theta(3));
    W(2,3) = dt * x(2) / (theta(2) * theta(3)^2) - dt * I / (theta(3)^2);

    Px = A * Px * A' + Qx;
    Ptheta = Ptheta + Qtheta;

    % Measurement Update
    OCV = polyval(p, x_pred(1));
    dOCV_dSoC = polyval(dp, x_pred(1));
    V_pred = OCV - theta(1)*I - x_pred(2);
    y_res = V - V_pred;

    Hx = [dOCV_dSoC, -1];
    Htheta = [-I, -x_pred(2), -1e-3*dt * I / (theta(3)^2)];

    Sx = Hx * Px * Hx' + R;
    Kx = Px * Hx' / Sx;
    x = x_pred + Kx * y_res;
    Px = (eye(2) - Kx * Hx) * Px;

    Stheta = Htheta * Ptheta * Htheta' + R;
    Ktheta = Ptheta * Htheta' / Stheta;
    theta = theta + Ktheta * y_res;
    theta(1) = max(theta(1), 1e-4);
    theta(2) = max(theta(2), 1e-3);
    theta(3) = max(theta(3), 100);
    Ptheta = (eye(3) - Ktheta * Htheta) * Ptheta;

    x_store(:,k) = x;
    theta_store(:,k) = theta;
    V_store(k) = V_pred;
end

%% Error Metrics
soc_est = x_store(1,:)';
soc_err = soc_true(2:end) - soc_est(2:end);
voltage_err = voltage(2:end) - V_store(2:end);

RMSE_soc = sqrt(mean(soc_err.^2));
MAE_soc = mean(abs(soc_err));
MAPE_soc = MAE_soc * 100;

N_clip = 342;
soc_err_clipped = soc_err(1:N_clip);
RMSE_soc_clipped = sqrt(mean(soc_err_clipped.^2));
MAE_soc_clipped = mean(abs(soc_err_clipped));
MAPE_soc_clipped = 100 * MAE_soc_clipped;

RMSE_V = sqrt(mean(voltage_err.^2));
MAPE_V = mean(abs(voltage_err ./ voltage(2:end))) * 100;

fprintf('--- SoC Estimation Errors ---\n');
fprintf('RMSE : %.4f\nMAE  : %.4f\nMAPE : %.2f %%\n\n', RMSE_soc, MAE_soc, MAPE_soc);
fprintf('--- SoC Estimation Errors (First %d Samples) ---\n', N_clip);
fprintf('RMSE : %.4f\nMAE  : %.4f\nMAPE : %.2f %%\n\n', RMSE_soc_clipped, MAE_soc_clipped, MAPE_soc_clipped);
fprintf('--- Voltage Tracking Errors ---\n');
fprintf('RMSE : %.4f V\nMAPE : %.2f %%\n\n', RMSE_V, MAPE_V);
fprintf('--- Final Estimated Parameters ---\n');
fprintf('R_0 = %.5f Ohm\nR_1 = %.5f Ohm\nC_1 = %.2f F\n', theta(1), theta(2), theta(3));

%% Plots
% SoC Estimation
figure;
plot(time, soc_true, 'k', time, soc_est, 'r--');
xlabel('Time (s)'); ylabel('SoC');
legend('True','Estimated'); title('SoC Estimation'); grid on;

% First 342 Samples
figure;
plot(time(1:N_clip), soc_true(1:N_clip), 'k', time(1:N_clip), soc_est(1:N_clip), 'r--');
xlabel('Time (s)'); ylabel('SoC');
legend('True','Estimated'); title(sprintf('SoC Estimation (First %d Samples)', N_clip)); grid on;

% Parameter Evolution
figure;
subplot(3,1,1); plot(time, theta_store(1,:)); ylabel('R_0 (\Omega)'); title('Estimated R_0'); grid on;
subplot(3,1,2); plot(time, theta_store(2,:)); ylabel('R_1 (\Omega)'); title('Estimated R_1'); grid on;
subplot(3,1,3); plot(time, theta_store(3,:)); ylabel('C_1 (F)'); xlabel('Time (s)'); title('Estimated C_1'); grid on;

% SoC Error
figure;
plot(time(2:end), soc_err);
xlabel('Time (s)'); ylabel('SoC Error');
title(sprintf('SoC Error | RMSE: %.4f | MAPE: %.2f %%', RMSE_soc, MAPE_soc));
grid on;

% Clipped SoC Error
figure;
plot(time(2:N_clip+1), soc_err_clipped, 'm');
xlabel('Time (s)'); ylabel('SoC Error');
title(sprintf('SoC Error (First %d Samples) | RMSE: %.4f | MAPE: %.2f %%', N_clip, RMSE_soc_clipped, MAPE_soc_clipped));
grid on;

% Voltage Tracking
figure;
plot(time, voltage, 'k', time, V_store, 'b--');
xlabel('Time (s)'); ylabel('Voltage (V)');
legend('Measured','Predicted');
title(sprintf('Voltage Tracking | RMSE: %.4f V | MAPE: %.2f %%', RMSE_V, MAPE_V));
grid on;

%% --- Compare OCV–SOC: Actual vs Polynomial Fit (using HPPC rest data) ---

% 1. Evaluate the polynomial OCV values at the actual SoC points from HPPC
ocv_fit = polyval(p, soc_unique);  % Fitted OCV values

% 2. Compute residuals
ocv_fit_error = ocv_unique - ocv_fit;

% 3. Error Metrics
RMSE_ocv = sqrt(mean(ocv_fit_error.^2));
MAE_ocv = mean(abs(ocv_fit_error));
MAPE_ocv = mean(abs(ocv_fit_error ./ ocv_unique)) * 100;

% 4. Print results
fprintf('--- OCV–SOC Polynomial Fit Errors (Using HPPC rest data) ---\n');
fprintf('RMSE : %.4f V\nMAE  : %.4f V\nMAPE : %.2f %%\n\n', RMSE_ocv, MAE_ocv, MAPE_ocv);

% 5. Plot actual vs fitted OCV–SOC curve
figure;
plot(soc_unique, ocv_unique, 'ko', 'DisplayName', 'Actual OCV (rest data)');
hold on;
plot(soc_unique, ocv_fit, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Polynomial Fit');
xlabel('SoC'); ylabel('OCV (V)');
legend(); grid on;
title(sprintf('OCV–SOC Fit | RMSE: %.4f V | MAPE: %.2f %%', RMSE_ocv, MAPE_ocv));

% 6. Plot residuals
figure;
plot(soc_unique, ocv_fit_error, 'b');
xlabel('SoC'); ylabel('OCV Error (V)');
title('OCV Fit Residuals (Actual - Fitted)');
grid on;

%results.soc_est= soc_est;
%results.current= zeros(380,1);
%results.voltage = polyval(p,soc_est);
%save('soc_estimatedresult.mat',"results");