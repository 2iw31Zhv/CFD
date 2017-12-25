function abs_A = abs_eig(A)

[~,~,length] = size(A);
abs_A = zeros(3, 3, length);

for i = 1 : length
    [R, Lambda] = eig(A(:, :, i));
    abs_A(:, :, i) = R * abs(Lambda) * R^-1;
end

end
