function F = evaluate_f(U, gamma)
% f1 = m
% f2 = m^2 / rho + (gamma - 1) (epsilon - m^2 / 2 / rho)
% f3 = m / rho * (epsilon + (gamma -1)(epsilon - m^2 / 2 / rho))
[~,length_U] = size(U);
F = zeros(3, length_U);
rho = U(1, :);
m = U(2, :);
epsilon = U(3, :);

F(1, :) = m;
F(2, :) = m.^2 ./ rho + (gamma - 1) * (epsilon - m.^2 ./2 ./ rho);
F(3, :) = m ./ rho .* (epsilon + (gamma - 1) * (epsilon - m.^2 ./ 2 ./ rho));
end