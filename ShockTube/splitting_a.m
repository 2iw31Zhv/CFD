function [A_positive, A_negative] = splitting_a(A)

[~,~,length] = size(A);
A_positive = zeros(3, 3, length);
A_negative = zeros(3, 3, length);

for i = 1 : length
    [R, Lambda] = eig(A(:, :, i));
    abs_A = R * abs(Lambda) * R^-1;
    A_positive(:, :, i) = 0.5 * (A(:, :, i) + abs_A);
    A_negative(:, :, i) = 0.5 * (A(:, :, i) - abs_A);
end
end