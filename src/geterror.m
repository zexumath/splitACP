% This script get the error of ACP and FD.
% Call metaltest2, metaltest5, metaltest6 before using.
%
% Ze Xu
% 06/06/17

eigFD   = eig(HessFD);
eigDFPT = eig(HessDFPT);
eigACP  = eig(HessACP);
errorFD  = zeros(opt1D.atom.NsCell,3);
errorACP = zeros(opt1D.atom.NsCell,3);
for I = 1:opt1D.atom.NsCell
    ref = norm(drhodR(:,I));
    errorFD(I,2)  = norm((rhoPerturbNew(:,I) - rho)/delta - drhodR(:,I))/ref;
    errorACP(I,1) = norm(drhodRbACP(:,I) - drhodRb(:,I))/ref;
    errorACP(I,2) = norm(drhodRACP(:,I)  - drhodR(:,I))/ref;
end
errorFD(:,3)    = eigFD  - eigDFPT;
errorACP(:,3)   = eigACP - eigDFPT;

