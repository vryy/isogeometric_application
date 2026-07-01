function [X, independent_cols, dependent_cols] = extractDependentColumns(A)

% 1. Run QR with vector permutation output 'p'
[Q, R, p] = qr(A, 'vector');

% 2. Recalculate your tolerance and rank
tol = max(size(A)) * eps(max(abs(diag(R))));
r = sum(abs(diag(R)) > tol); % This will be 96

% The first 96 columns are linearly independent base columns
independent_cols = p(1:r);

% The last 9 columns are the dependent columns
dependent_cols = p(r+1:end);

disp('The 9 dependent columns in original A are:');
disp(dependent_cols);

% Extract the blocks from R
R11 = R(1:r, 1:r);
R12 = R(1:r, r+1:end);

% Solve R11 * X = R12 to get the combinations
X = R11 \ R12;

% Clean up tiny rounding noises
X(abs(X) < tol) = 0;

end
