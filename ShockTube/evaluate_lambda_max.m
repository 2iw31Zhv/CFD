function lambda_max = evaluate_lambda_max(U, gamma)
rho = U(1, :);
m = U(2, :);
epsilon = U(3, :);

u = m ./ rho;
E = epsilon ./ rho;
[~, length_U] = size(U);
lambda_max = zeros(1, length_U);

for i = 1 : length_U
    A = [0, 1, 0;
        0.5 * (gamma - 3) * u(i).^2, (3 - gamma) .* u(i), gamma - 1;
        (gamma - 1)*u(i)^3 - gamma * u(i) * E(i), -1.5 * (gamma - 1) * u(i)^2 + gamma * E(i), gamma * u(i)];
    eigen_value = eig(A);
    lambda_max(i) = max(abs(eigen_value));
end

end