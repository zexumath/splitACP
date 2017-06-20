% This script test on the error when Npole and Ntot is changed.
% Ncheb = 4 fixed.

opt1D = initOptnlpp4m();

% opt1D.temperature = 20000; Natom = 60;

metaltest2;
metaltest6;


Nptest = 20:4:40;
Ntottest = 0.6:0.1:1.1;
errNpNtot = zeros(6,6);

for inp = 1:6
    for intot = 1:6
        opt1D.NeExtra = opt1D.Ne * Ntottest(intot);
        opt1D.Npole   = Nptest(inp);
        metaltest5;
        geterror;
        errNpNtot(inp,intot) = max(abs(errorACPsp(:,3)));
    end
end
