function flux = finite_volume_flux(UL, UR, abs_A, gamma)

[~, length] = size(UL);
flux = zeros(3, length);
for j = 1 : length
    flux(:, j) = 0.5 * (evaluate_f(UL(:, j), gamma)...
        + evaluate_f(UR(:, j), gamma))...
        - 0.5 * abs_A(:, :, j) * (UR(:, j) - UL(:, j));
end

end