clear all; close all; 

load('COVID_data.mat'); 

% A time vector based on data length [cite: 25, 33]
t = 1:length(Australia.daily_new_cases);

% Set the split date for the two waves (around May 25th) [cite: 47]
split_idx = 150; 

% Organize data for the subplot loop 
% datasets = {Australia.daily_new_cases, Australia.daily_new_death, Australia.currently_infected};
titles = {'Daily New Cases', 'Daily New Deaths', 'Currently Infected'};

% Peak Fitting 
figure('Name', 'Australia COVID-19 Data Fitting');

for i = 1:3
    % Create a 3x1 grid of plots 
    subplot(3, 1, i); 
    y_data = datasets{i};

    
    % Splits data into two separate peaks 
    t1 = t(1:split_idx);         d1 = y_data(1:split_idx);
    t2 = t(split_idx+1:end);     d2 = y_data(split_idx+1:end);
    
    % Define Gaussian equation: Amplitude * exp(-((t - Peak)/Width)^2) 
    model = @(p, t) p(1) * exp(-((t - p(2)) / p(3)).^2);
    
    % Fit the first peak
    fit1 = lsqcurvefit(model, [max(d1), 50, 20], t1, d1);
    
    % Fit the second peak using the same equation 
    fit2 = lsqcurvefit(model, [max(d2), 180, 30], t2, d2);
    
    % Plot raw data points and the two fitting curves
    plot(t, y_data, 'k.'); hold on;
    plot(t1, model(fit1, t1), 'r-', 'LineWidth', 2);
    plot(t2, model(fit2, t2), 'b-', 'LineWidth', 2);
    
    title(titles{i});
    legend('Data', 'Peak 1 Fit', 'Peak 2 Fit');
    grid on;
end

% SIR Epidemic Model Fitting
% Fit only the second peak using differential equations 
t_sir = t(split_idx+1:end) - split_idx; 
y_sir = Australia.currently_infected(split_idx+1:end);

% Population constant for Australia 
N = 25e6; 

% Optimize Beta (transmission) and Gamma (recovery) rates
p_start = [0.2, 0.1]; 
p_final = lsqcurvefit(@(p, t) solve_sir(p, t, N), p_start, t_sir, y_sir);

% Plot the SIR model results 
[~, sir_curve] = solve_sir(p_final, t_sir, N);
figure('Name', 'SIR Model Bonus');
plot(t_sir, y_sir, 'go', t_sir, sir_curve, 'm-', 'LineWidth', 2);
title('Bonus: SIR Model Fit for Peak 2');
legend('Actual', 'SIR Fit');

% ODEs
function [t_out, I_out] = solve_sir(p, t, N)
    % Initial infected count at the start of the second wave
    I0 = 500; 
    % Define the SIR differential equations [cite: 115]
    sir_eq = @(t, y) [
        -p(1)*y(1)*y(2)/N;               % dS/dt (Susceptible)
        (p(1)*y(1)*y(2)/N) - p(2)*y(2);  % dI/dt (Infected)
        p(2)*y(2)];                      % dR/dt (Recovered)
    
    % Solve the ODE system over time 
    [t_out, soln] = ode45(sir_eq, t, [N-I0, I0, 0]);
    I_out = soln(:, 2); % Return only the Infected group
end
