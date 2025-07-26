    clear; clc; close all;
    
    %% Load measurement data
    load("03-11-17_08.47 25degC_5Pulse_HPPC_Pan18650PF.mat");  % Replace with your file
    time = meas.Time;
    current = -meas.Current;    % A
    voltage = meas.Voltage;     % V
    soc_true = (meas.Ah - meas.Ah(end)) / 2.9; % true SoC
    
    dt = mean(diff(time));
    N = length(time);
    
    %% Estimate OCV–SOC relationship using 5th-degree polynomial
    p = polyfit(soc_true, voltage, 5);       % Fit
    dp = polyder(p);                         % Derivative of OCV(SOC)
    
    %% Battery capacity
    Q = 2.9 * 3600;  % [As]
    
    %% Initial values
    % Initial state estimates
    x = [soc_true(1); 0];             % [SoC; V_RC]
    theta = [0.005; 0.005; 3000];     % [R0; R1; C1]
    
    % Covariances
    Px = diag([1e-5, 1e-5]);          % Smaller initial state uncertainty
    Ptheta = diag([0.001, 0.005, 100]); % Parameter uncertainty
    
    Qx = diag([1e-7, 1e-6]);          % Process noise for x
    Qtheta = diag([1e-7, 1e-6, 1]);   % Parameter process noise
    R = 1e-3;                         % Measurement noise
    
    
    %% Storage
    x_store = zeros(2, N);
    theta_store = zeros(3, N);
    V_store = zeros(N, 1);
    
    %% DEKF Loop
    for k = 2:N
        I = current(k-1);
        V = voltage(k);
    
        % --- Prediction ---
        SoC_pred = x(1) - dt * I / Q;
        V_RC_pred = x(2) + dt * (-x(2)/(theta(2)*theta(3)) + I/theta(3));
        R1 = theta(2); C1 = theta(3);
        x_pred = [SoC_pred; V_RC_pred];
    
        A = [1, 0;
             0, 1 - dt/(theta(2)*theta(3))];
    
        W = zeros(2,3);
        W(2,2) = dt * x(2) / (theta(2)^2 * theta(3));
        W(2,3) = dt * x(2) / (theta(2) * theta(3)^2) - dt * I / (theta(3)^2);
    
        Px = A * Px * A' + Qx;
        Ptheta = Ptheta + Qtheta;
    
        % --- Measurement ---
        OCV = polyval(p, x_pred(1));
        dOCV_dSoC = polyval(dp, x_pred(1));
    
        V_pred = OCV - theta(1)*I - x_pred(2);
        y_pred = V_pred;
        V_store(k) = y_pred;
    
        Hx = [dOCV_dSoC, -1];         % dy/dx
        Htheta = [-I, -1e-5*x_pred(2), -1*dt * I / (theta(3)^2)];          % dy/dtheta
    
        % Kalman gain
        Sx = Hx * Px * Hx' + R;
        Kx = Px * Hx' / Sx;
    
        Stheta = Htheta * Ptheta * Htheta' + R;
        Ktheta = Ptheta * Htheta' / Stheta;
    
        % Update
        y_res = V - y_pred;
        x = x_pred + Kx * y_res;
        Px = (eye(2) - Kx * Hx) * Px;
    
        theta = theta + Ktheta * y_res;
        Ptheta = (eye(3) - Ktheta * Htheta) * Ptheta;
    
        % Store
        x_store(:,k) = x;
        theta_store(:,k) = theta;
    end
    
    %% Compute Errors (exclude first point)
    soc_est = x_store(1,:)';
    soc_true = soc_true(:);
    voltage = voltage(:);
    V_store = V_store(:);
    
    soc_err = soc_true(2:end) - soc_est(2:end);
    voltage_err = voltage(2:end) - V_store(2:end);
    
    RMSE_soc = sqrt(mean(soc_err.^2));
    MAE_soc = mean(abs(soc_err));
    MAPE_soc = MAE_soc * 100;
    
    RMSE_V = sqrt(mean(voltage_err.^2));
    MAPE_V = mean(abs(voltage_err ./ voltage(2:end))) * 100;
    
    %% Display Metrics
    fprintf('--- SoC Estimation Errors ---\n');
    fprintf('RMSE : %.4f\n', RMSE_soc);
    fprintf('MAE  : %.4f\n', MAE_soc);
    fprintf('MAPE : %.2f %%\n\n', MAPE_soc);
    
    fprintf('--- Voltage Tracking Errors ---\n');
    fprintf('RMSE : %.4f V\n', RMSE_V);
    fprintf('MAPE : %.2f %%\n\n', MAPE_V);
    
    fprintf('--- Final Estimated Parameters ---\n');
    fprintf('R_0 = %.5f Ohm\n', theta(1));
    fprintf('R_1 = %.5f Ohm\n', theta(2));
    fprintf('C_1 = %.2f F\n', theta(3));
    
    %% Plot SoC
    figure;
    plot(time, soc_true, 'k', time, soc_est, 'r--');
    xlabel('Time (s)'); ylabel('SoC');
    legend('True','Estimated');
    title('SoC Estimation');
    grid on;
    
    %% Plot Parameter Estimates
    figure;
    subplot(3,1,1); plot(time, theta_store(1,:)); ylabel('R_0 (\Omega)'); title('Estimated R_0'); grid on;
    subplot(3,1,2); plot(time, theta_store(2,:)); ylabel('R_1 (\Omega)'); title('Estimated R_1'); grid on;
    subplot(3,1,3); plot(time, theta_store(3,:)); ylabel('C_1 (F)'); xlabel('Time (s)'); title('Estimated C_1'); grid on;
    
    %% Plot SoC Error
    figure;
    plot(time(2:end), soc_err);
    xlabel('Time (s)'); ylabel('SoC Error');
    title(sprintf('SoC Error | RMSE: %.4f | MAPE: %.2f %%', RMSE_soc, MAPE_soc));
    grid on;
    
    %% Plot Voltage Tracking
    figure;
    plot(time, voltage, 'k', time, V_store, 'b--');
    xlabel('Time (s)'); ylabel('Voltage (V)');
    legend('Measured','Predicted');
    title(sprintf('Voltage Tracking | RMSE: %.4f V | MAPE: %.2f %%', RMSE_V, MAPE_V));
    grid on;
%% --- OCV–SOC Fit Error Analysis from HPPC Data (Rest Periods) ---
% Ensure these variables are available: soc_all, ocv_all

% Fit polynomial to OCV vs SoC
ocv_fit = polyval(p, soc_true);
ocv_fit_error = voltage - ocv_fit;

% Error metrics
RMSE_ocv = sqrt(mean(ocv_fit_error.^2));
MAE_ocv = mean(abs(ocv_fit_error));
MAPE_ocv = mean(abs(ocv_fit_error ./ voltage)) * 100;

fprintf('--- OCV–SOC Polynomial Fit Errors (Using HPPC rest data) ---\n');
fprintf('RMSE : %.4f V\n', RMSE_ocv);
fprintf('MAE  : %.4f V\n', MAE_ocv);
fprintf('MAPE : %.2f %%\n\n', MAPE_ocv);

% OCV–SOC curve plot
figure;
plot(soc_true, voltage, 'ko', 'DisplayName', 'Actual OCV (rest data)');
hold on;
plot(soc_true, ocv_fit, 'r-', 'LineWidth', 1.5, 'DisplayName', 'Polynomial Fit');
xlabel('SoC'); ylabel('OCV (V)');
legend(); grid on;
title(sprintf('OCV–SOC Fit | RMSE: %.4f V | MAPE: %.2f %%', RMSE_ocv, MAPE_ocv));

% Residuals plot
figure;
plot(soc_true, ocv_fit_error, 'b');
xlabel('SoC'); ylabel('OCV Error (V)');
title('OCV Fit Residuals (Actual - Fitted)');
grid on;
