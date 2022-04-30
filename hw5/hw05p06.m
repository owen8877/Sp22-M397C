clear; 


k = 2; n = 5;

At_12 = randn(2*n, k) * randn(k, 2*n);
At_21 = randn(2*n, k) * randn(k, 2*n);

A_12 = randn(n, k) * randn(k, n);
A_21 = randn(n, k) * randn(k, n);
A_34 = randn(n, k) * randn(k, n);
A_43 = randn(n, k) * randn(k, n);

A_11 = randn(n);
A_22 = randn(n);
A_33 = randn(n);
A_44 = randn(n);

fprintf('Rank of A_12: %d\n', rank(A_12))
fprintf('Rank of A_21: %d\n', rank(A_21))
fprintf('Rank of A_34: %d\n', rank(A_34))
fprintf('Rank of A_43: %d\n', rank(A_43))
fprintf('Rank of A_12 (intermediate): %d\n', rank(At_12))
fprintf('Rank of A_21 (intermediate): %d\n', rank(At_21))
A = [[A_11, A_12; A_21, A_22], At_12; At_21, [A_33, A_34; A_43, A_44]];
B = inv(A);
fprintf('Rank of B_12: %d\n', rank(B(1:n, n+1:2*n)))
fprintf('Rank of B_21: %d\n', rank(B(n+1:2*n, 1:n)))
