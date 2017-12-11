function A = evaluate_a(U, gamma)
rho = U(1, :);
m = U(2, :);
epsilon = U(3, :);

[~, length_U] = size(U);
u = m ./ rho;
E = epsilon ./ rho;

A = zeros(3, 3, length_U);
A(1, 1, :) = zeros(1, length_U);
A(1, 2, :) = ones(1, length_U);
A(1, 3, :) = zeros(1, length_U);
A(2, 1, :) = 0.5 * (gamma - 3) * u.^2;
A(2, 2, :) = (3 - gamma) .* u;
A(2, 3, :) = (gamma - 1) .* ones(1, length_U);
A(3, 1, :) = (gamma - 1) * u.^3 - gamma * u .* E;
A(3, 2, :) = -1.5 * (gamma - 1) * u.^2 + gamma * E;
A(3, 3, :) = gamma * u;

end