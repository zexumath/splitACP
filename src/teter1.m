function M = teter1(kpt, NsGrid)
[m,~] = size(kpt);
prec = zeros(m,1);
for i =1:m
    a = kpt(i,1)^2 ;
    b = 27 + a*(18+a*(12+a*8));
    prec(i) = b/(b+16*a^4);
end

PrecMat = spdiags(prec,0,NsGrid, NsGrid);

M = @(x) real(ifft(PrecMat*fft(x)));
