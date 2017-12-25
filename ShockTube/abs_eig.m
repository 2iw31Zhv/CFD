function [abs_A, L] = abs_eig(A)

[~,~,length] = size(A);
abs_A = zeros(3, 3, length);
L = zeros(3, 3, length);

for i = 1 : length
    [lmatrix, Lambda] = eig(A(:, :, i));
    abs_A(:, :, i) = lmatrix * abs(Lambda) * lmatrix^-1;
    L(:, :, i) = lmatrix;
end

end
