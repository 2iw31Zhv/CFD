function lambda_max = evaluate_lambda_max(U, gamma)
rho = U(1, :);
m = U(2, :);
epsilon = U(3, :);

u = m ./ rho;
E = epsilon ./ rho;
[~, length_U] = size(U);
lambda_max = zeros(1, length_U);

A = evaluate_a(U, gamma);
for i = 1 : length_U
    eigen_value = eig(A(:, :, i));
    lambda_max(i) = max(abs(eigen_value));
end

end