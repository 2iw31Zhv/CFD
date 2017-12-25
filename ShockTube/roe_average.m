function U = roe_average(UL, UR)
rhoL = UL(1, :);
rhoR = UR(1, :);

mL = UL(2, :);
mR = UR(2, :);

epsilonL = UL(3, :);
epsilonR = UR(3, :);

[~, length_U] = size(UL);

U = zeros(3, length_U);
U(1, :) = sqrt(rhoL .* rhoR);
U(2, :) = (sqrt(rhoR) .* mL + sqrt(rhoL) .* mR) ./ (sqrt(rhoL) + sqrt(rhoR));
U(3, :) = (sqrt(rhoR) .* epsilonL + sqrt(rhoL) .* epsilonR)...
    ./ (sqrt(rhoL) + sqrt(rhoR));
end