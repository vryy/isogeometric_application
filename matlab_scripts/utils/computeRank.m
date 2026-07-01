function true_rank = computeRank(A)

% 1. Compute the sparse QR factorization
[Q, R, P] = qr(A);

% 2. Set a strict numerical tolerance
tol = max(size(A)) * eps(max(abs(diag(R))));

% 3. Count how many diagonal values are actually non-zero
true_rank = sum(abs(diag(R)) > tol);

end
