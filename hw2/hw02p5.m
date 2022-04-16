m=12; n=20; b=3; k=5;
A = rand(m,n);
Gc = randn(n,k);
Gr = randn(m,k);
Yc = zeros(m,k);
Yr = zeros(n,k);
for i = 1:ceil(n/b)
    ind = (i-1)*b + (1:min(b, n-(i-1)*b));
    A_slice = A(:,ind);
    Yc = Yc + A_slice*Gc(ind,:);
    Yr(ind, :) = A_slice' * Gr;
end
fprintf(1,'Error in Yc = %12.3e\n',max(max(abs(Yc - A *Gc))))
fprintf(1,'Error in Yr = %12.3e\n',max(max(abs(Yr - A'*Gr))))