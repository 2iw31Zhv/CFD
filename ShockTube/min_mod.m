function slope = min_mod(L, U_m1, U, U_p1)

[~, length] = size(U);

slope = zeros(3, length);
for i = 1 : length
    eigen_vec_left = L(:, :, i) * (U_p1(:, i) - U(:, i));
    eigen_vec_right = L(:, :, i) * (U(:, i) - U_m1(:, i));
    for j = 1 : 3
        a = eigen_vec_left(j);
        b = eigen_vec_right(j);
        if a * b < 0
            slope(j, i) = 0;
        else
            slope(j, i) = sign(a) * min(abs(a), abs(b));
        end
    end
end

end