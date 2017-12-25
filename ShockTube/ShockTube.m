%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Shock Tube
% Author: Ziwei Zhu
% Email: ziweizhu95@gmail.com
% Date: 20171225
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
% m = rho u 
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
mode = 'Roe-order-1';
% The following modes are supported:
% 'Rusanov', 'Jameson', 'FVS-order-1', 'FVS-order-2'
% 'Roe-order-1', 'FVD-order-2'

% Discretization
nx = 400;
nt = 2000;

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
    rho = U(1, :);
    m = U(2, :);
    epsilon = U(3, :);
    
    u = m ./ rho;
    E = epsilon ./ rho;
    p = (rho .* E - 0.5 * rho .* u.^2) * (gamma - 1);

    p_begin = (U_begin(3,:) - 0.5 * U_begin(2,:) .^2 ./ U_begin(1,:) )...
        * (gamma - 1);
    p_end = (U_end(3,:) - 0.5 * U_end(2,:).^2 ./ U_end(1,:) )...
        * (gamma - 1);
    
    p_p1 = circshift(p, [0,-1]); p_p1(:, end) = p_end;
    p_m1 = circshift(p, [0, 1]); p_m1(:, 1) = p_begin;
    
    t = i * t_step;
    
    lambda_max = evaluate_lambda_max(U, gamma);
    lambda_begin = evaluate_lambda_max(U_begin, gamma);
    lambda_end = evaluate_lambda_max(U_end, gamma);
    
    lambda_p1 = circshift(lambda_max, [0,-1]); lambda_p1(end) = lambda_end;
    lambda_m1 = circshift(lambda_max, [0,1]); lambda_m1(1) = lambda_begin;
    
    
    F = evaluate_f(U, gamma);
    F_begin = evaluate_f(U_begin, gamma);
    F_end = evaluate_f(U_end, gamma);
    
    F_p1 = circshift(F, [0,-1]); F_p1(:, end) = F_end;
    F_m1 = circshift(F, [0,1]); F_m1(:, 1) = F_begin;
    
    U_p1 = circshift(U, [0,-1]); U_p1(:, end) = U_end;
    U_p2 = circshift(U, [0,-2]); U_p2(:, end) = U_end; U_p2(:, end-1) = U_end;
    
    U_m1 = circshift(U, [0,1]); U_m1(:, 1) = U_begin;
    U_m2 = circshift(U, [0,2]); U_m2(:, 1) = U_begin; U_m2(:, 2) = U_begin;
    
    F_positive = zeros(3, nx+1);
    F_negative = zeros(3, nx+1);
    if strcmp(mode, 'Rusanov')
        lambda_positive = 0.5 * (lambda_max + lambda_p1); % 1 x 101
        lambda_negative = 0.5 * (lambda_max + lambda_m1); % 1 x 101
    
    
        F_positive = 0.5 * (F_p1 + F)...
            -  0.5 * lambda_positive .* (U_p1 - U);
        F_negative = 0.5 * (F_m1 + F)...
            -  0.5 * lambda_negative .* (U - U_m1);
    elseif strcmp(mode, 'Jameson')
        lambda_positive = 0.5 * (lambda_max + lambda_p1); % 1 x 101
        lambda_negative = 0.5 * (lambda_max + lambda_m1); % 1 x 101
        
        nu = abs(p_p1 - 2 * p + p_m1) ./ abs(p_p1 + 2 * p + p_m1);
        nu_begin = abs(p(1) - p_begin) ./ abs(p(1) + 3 * p_begin);
        nu_end = abs(p(end) - p_end) ./ abs(3 * p_end + p(end - 1));
        
        nu_p2 = circshift(nu, [0,-2]); nu_p2(end) = 0;
        nu_p2(end-1) = nu_end;
        nu_p1 = circshift(nu, [0,-1]); nu_p1(end) = nu_end;
        
        nu_m1 = circshift(nu, [0,1]); nu_m1(1) = nu_begin;
        nu_m2 = circshift(nu, [0,2]); nu_m2(1) = 0;
        nu_m2(2) = nu_begin;
        
        e2_positive = 4.1 * max([nu_p2;nu_p1;nu;nu_m1]);
        %e2_positive = 0.5 * ones(1, nx + 1);
        e4_positive = max([zeros(1,nx+1);
        1 / 64.0 * ones(1,nx+1) - e2_positive]);
    
        e2_negative = 4.1 * max([nu_p1;nu;nu_m1;nu_m2]);
        %e2_positive = 0.5 * ones(1, nx+1);
        e4_negative = max([zeros(1,nx+1);
        1 / 64.0 * ones(1,nx+1) - e2_negative]);
        
        F_positive = 0.5 * (F_p1 + F)...
            - lambda_positive .* e2_positive .* (U_p1 - U)...
            + lambda_positive .* e4_positive .* (U_p2 - 3 * U_p1 + 3 * U - U_m1);
        F_negative = 0.5 * (F_m1 + F)...
            - lambda_negative .* e2_negative .* (U - U_m1)...
            + lambda_negative .* e4_negative .* (U_p1 - 3 * U + 3 * U_m1 - U_m2);
    elseif strcmp(mode, 'FVS-order-1')
        A = evaluate_a(U, gamma);
        A_begin = evaluate_a(U_begin, gamma);
        A_end = evaluate_a(U_end, gamma);
        
        [A_positive, A_negative] = splitting_a(A);
        
        F_pos = zeros(3, nx+1);
        F_neg = zeros(3, nx+1);
        
        for j = 1 : nx+1
            F_pos(:, j) = A_positive(:, :, j) * U(:, j);
            F_neg(:, j) = A_negative(:, :, j) * U(:, j);
        end
        
        [A_positive_begin, A_negative_begin] = splitting_a(A_begin);
        [A_positive_end, A_negative_end] = splitting_a(A_end);
        
        F_pos_m1 = circshift(F_pos, [0, 1]); 
        F_pos_m1(:, 1) = A_positive_begin * U_begin;
        
        F_neg_p1 = circshift(F_neg, [0, -1]);
        F_neg_p1(:, end) = A_negative_end * U_end;
        
        F_positive = F_pos + F_neg_p1;
        F_negative = F_neg + F_pos_m1;
    elseif strcmp(mode, 'FVS-order-2')
        A = evaluate_a(U, gamma);
        A_begin = evaluate_a(U_begin, gamma);
        A_end = evaluate_a(U_end, gamma);
        
        [A_positive, A_negative] = splitting_a(A);
        
        F_pos = zeros(3, nx+1);
        F_neg = zeros(3, nx+1);
        
        for j = 1 : nx+1
            F_pos(:, j) = A_positive(:, :, j) * U(:, j);
            F_neg(:, j) = A_negative(:, :, j) * U(:, j);
        end
        
        [A_positive_begin, A_negative_begin] = splitting_a(A_begin);
        [A_positive_end, A_negative_end] = splitting_a(A_end);
        
        F_pos_m1 = circshift(F_pos, [0, 1]); 
        F_pos_m1(:, 1) = A_positive_begin * U_begin;
        
        F_pos_m2 = circshift(F_pos, [0, 2]);
        F_pos_m2(:, 1) = A_positive_begin * U_begin;
        F_pos_m2(:, 2) = A_positive_begin * U_begin;
        
        F_neg_p1 = circshift(F_neg, [0, -1]);
        F_neg_p1(:, end) = A_negative_end * U_end;
        
        F_neg_p2 = circshift(F_neg, [0, -2]);
        F_neg_p2(:, end) = A_negative_end * U_end;
        F_neg_p2(:, end-1) = A_negative_end *  U_end;
        
        F_positive = (1.5 * F_pos - 0.5 * F_pos_m1)...
            + (1.5 * F_neg_p1 - 0.5 * F_neg_p2);
        F_negative = (1.5 * F_pos_m1 - 0.5 * F_pos_m2)...
            + (1.5 * F_neg - 0.5 * F_neg_p1);
    elseif strcmp(mode, 'Roe-order-1')
        U_positive = roe_average(U, U_p1);
        U_negative = roe_average(U_m1, U);
        A_positive = evaluate_a(U_positive, gamma);
        A_negative = evaluate_a(U_negative, gamma);
        
        for j = 1 : nx + 1
            F_positive(:, j) = 0.5 * (evaluate_f(U(:, j), gamma)...
                + evaluate_f(U_p1(:, j), gamma))...
                - 0.5 * abs_eig(A_positive(:,:, j)) * (U_p1(:, j) - U(:, j));
            F_negative(:, j) = 0.5 * (evaluate_f(U_m1(:, j), gamma)...
                + evaluate_f(U(:, j), gamma))...
                - 0.5 * abs_eig(A_negative(:,:, j)) * (U(:, j) - U_m1(:, j));
        end
    elseif strcmp(mode, 'FVD-order-2')
        
    end
    
    % time: forward diff
    U = U - t_step / x_step * (F_positive - F_negative);   
    x = -0.5:x_step:0.5;
    
    
    % draw figure
    cla;
    plot(x, rho);
    hold on
    plot(x, u);
    hold on
    plot(x, p);
    hold on   
    legend('\rho', 'u', 'p');
    title_str = sprintf('Shock Tube, mode = %s, nx = %d, nt = %d, t = %f s', mode, nx, nt, t);
    title(title_str);
    axis([-0.5, 0.5, 0.0, 2.0]);
    pause(t_step * 0.01);
end

