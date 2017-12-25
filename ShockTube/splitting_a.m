function [A_positive, A_negative] = splitting_a(A)

abs_A = abs_eig(A);
A_positive= 0.5 * (A + abs_A);
A_negative = 0.5 * (A - abs_A);

end