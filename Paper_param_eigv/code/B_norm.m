function norm_B = B_norm(x, B)
% Calculates the B-norm of x with B a symmetric, positive definite matrix:
% It is defined as sqrt(x' B x)
% INPUT:
%   (*) x: vector from which we calculate the B-norm
%   (*) B: positive definite, symmetric matrix
% OUTPUT:
%   (*) norm_B: float representing the B-norm of x
norm_B = sqrt(x'*B*x);
end