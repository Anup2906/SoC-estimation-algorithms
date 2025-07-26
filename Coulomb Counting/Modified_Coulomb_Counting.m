clc; clear; close all;

% Load CSV file
load("03-18-17_02.17 25degC_Cycle_1_Pan18650PF.mat")

% Extract required columns
time = meas.Time;                         % Time in seconds
voltage = meas.Voltage;                  % Voltage in Volts
current = meas.Current;                  % Current in Amperes
nominal_capacity = 2.9;
reamining_capacity = meas.Ah;
reamining_capacity = reamining_capacity - reamining_capacity(end);
soc = reamining_capacity/nominal_capacity;                          % Depth of Discharge

% Actual SOC from file
actual_soc = soc * 100;
   

% Battery Parameters
current_cycle = 5000;
total_cycle = 10000;
factor = current_cycle / total_cycle;
%Q_nominal = 3474.369004*2.8992/3600;     % Nominal battery capacity in Ah
Q_nominal=2.9;
eta_c = 0.97;                            % Charging efficiency
eta_d = 0.97;                            % Discharging efficiency
Vlim = 2.499;                            % Voltage limit for discharge
time_step = mean(diff(time));                          % Time step in seconds
initial_SOC = actual_soc(1);                % Initial State of Charge

% Initialize calculated SOC array
calculated_soc = zeros(length(time), 1);
calculated_soc(1) = initial_SOC;         % Set initial SOC

% Initialize capacity_used based on initial SOC
%capacity_used = Q_nominal - (initial_SOC * Q_nominal / 100);
initial_capacity=initial_SOC * Q_nominal / 100;
% Main loop for Coulomb Counting
for i = 2:length(time)
    Ib = current(i);
    Vb = voltage(i);
    prev_SOC = calculated_soc(i-1);      % SOC from previous time step

    if Ib == 0
        SOC = prev_SOC - 0.000001 * time_step * 100; % Self-discharge
    elseif Ib < 0
        if Vb > Vlim

            initial_capacity= initial_capacity + (Ib * time_step / 3600) * eta_d;
           
            SOC = (initial_capacity / Q_nominal) * 100;     
        else
            SOC = 0;
            initial_capacity = 0;
        end
    elseif Ib > 0
            initial_capacity= initial_capacity + (Ib * time_step / 3600) * eta_c;
            SOC = (initial_capacity / Q_nominal) * 100;      
     %else
        %if prev_SOC >= 100
            %SOC = 100;
         %end
    end

    % Clamp SOC between 0 and 100
    SOC = max(0, min(100, SOC));
    actual_soc = max(0, min(100, actual_soc));

    % Store calculated SOC
    calculated_soc(i) = SOC;
end

% Error between actual and calculated SOC
error_soc = actual_soc - calculated_soc;

% Plotting results
figure;
plot(time,current,'k', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Current(A)');
title('current vs time');

figure;
%subplot(2,1,1);
plot(time, actual_soc, 'b', 'LineWidth', 1.5); hold on;
plot(time, calculated_soc, 'r--', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('SOC (%)');
legend('Actual SOC', 'Calculated SOC');
title('Actual vs Calculated SOC');
figure;
%subplot(2,1,2);
plot(time, error_soc, 'k', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('SOC Error (%)');
title('SOC Estimation Error (Actual - Calculated)');

figure;
plot(time, current, 'g', 'LineWidth', 1.5); hold on;
plot(time, voltage, 'm', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Current (A) / Voltage (V)');
legend('Current', 'Voltage');
title('Battery Current and Voltage');

% Calculate error statistics
max_error = max(error_soc);
min_error = min(error_soc);
mean_error = mean(error_soc);

% Display results
fprintf('Maximum SOC Error: %.4f%%\n', max_error);
fprintf('Minimum SOC Error: %.4f%%\n', min_error);
fprintf('Mean SOC Error: %.4f%%\n', mean_error);


grid on;
box off;
set (gca, 'LineWidth', 1.2);