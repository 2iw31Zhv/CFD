%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Shock Tube
% Author: Ziwei Zhu
% Email: ziweizhu95@gmail.com
% Date: 20171206
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% One-Dimensional Euler Equation
% For ideal gas
% rho E = rho e + 0.5 * rho * u^2 = p / (gamma - 1) + 0.5 * rho * u^2
% p = rho E - 0.5 * rho * u^2 * (gamma - 1);
% Euler Equation
%        rho            f1
% D/Dt [  m  ] + D/Dx [ f2 ] = 0
%      epsilon          f3
% where:
% m = rho E 
% epsilon = rho E
% f1 = m
% f2 = m^2 / rho + (gamma - 1) (epsilon - m^2 / 2 / rho)
% f3 = m / rho * (epsilon + (gamma -1)(epsilon - m^2 / 2 / rho))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Problem Settings
x0 = -0.5;
x1 = 0.5;
t_end = 0.25;
gamma = 1.4;
mode = 'Rusanov'; % 'Rusanov', 'Jameson', 'FVS-order-1', 'FVS-order-2'


% Discretization
nx = 1000;
nt = 1000;

% Intermediate variables
x_step = (x1 - x0) / nx;
t_step = t_end / nt;

% Initialize
U = zeros(3, nx + 1);
U_begin = [1; 0.75; 1];
U_end = [0.125; 0; 0.1];

for i = 0 : nx
    x = x0 + x_step * i;
    if x <= 0
        U(:,i+1) = U_begin;
    else
        U(:,i+1) = U_end;
    end
end

t = 0;
for i = 0 : nt
    t = nt * t_step;
    lambda_max = evaluate_lambda_max(U, gamma);
    lambda_begin = evaluate_lambda_max(U_begin, gamma);
    lambda_end = evaluate_lambda_max(U_end, gamma);
    
    F = evaluate_f(U, gamma);
    F_begin = evaluate_f(U_begin, gamma);
    F_end = evaluate_f(U_end, gamma);
    
    % Rusanov
    lambda_positive = 0.5 * (lambda_max + circshift(lambda_max, [0,-1])); % 1 x 101
    lambda_positive(end) = 0.5 * (lambda_max(end) + lambda_end);
    lambda_negative = 0.5 * (lambda_max + circshift(lambda_max, [0, 1])); % 1 x 101
    lambda_negative(1) = 0.5 * (lambda_max(1) + lambda_begin);
    
    
    F_positive = 0.5 * (circshift(F, [0, -1]) + F)...
        -  0.5 * lambda_positive .* (circshift(U, [0, -1]) - U);
    F_positive(:, end) = 0.5 * (F_end + F(:, end))...
        -  0.5 * lambda_positive(:, end) .* (U_end - U(:, end));
    F_negative = 0.5 * (circshift(F, [0, 1]) + F)...
        -  0.5 * lambda_negative .* (U - circshift(U, [0, 1]));
    F_negative(:, 1) = 0.5 * (F_begin + F(:, 1))...
        -  0.5 * lambda_negative(:, 1) .* (U(:, 1) - U_begin);
    
    % time: forward diff
    U = U - t_step / x_step * (F_positive - F_negative);
    
    rho = U(1, :);
    m = U(2, :);
    epsilon = U(3, :);
    
    u = m ./ rho;
    E = epsilon ./ rho;
    p = (rho .* E - 0.5 * rho .* u.^2) * (gamma - 1);
    
    x = -0.5:x_step:0.5;
    cla;
    plot(x, rho);
    hold on
    plot(x, u);
    hold on
    plot(x, p);
    hold on
    pause(0.01);
end

