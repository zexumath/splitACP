% This script get the error of ACP and FD.
% Call metaltest2, metaltest5, metaltest6 before using.
%
% Ze Xu
% 06/06/17

eigFD   = eig(HessFD);
eigDFPT = eig(HessDFPT);
eigACPsp  = eig(HessACPsp);
errorFD  = zeros(opt1D.atom.NsCell,3);
errorACPsp = zeros(opt1D.atom.NsCell,3);
for I = 1:opt1D.atom.NsCell
    ref = norm(drhodR(:,I));
    errorFD(I,2)  = norm((rhoPerturbNew(:,I) - rho)/delta - drhodR(:,I))/ref;
    errorACPsp(I,1) = norm(drhodRbACPsp(:,I) - drhodRb(:,I))/ref;
    errorACPsp(I,2) = norm(drhodRACPsp(:,I)  - drhodR(:,I))/ref;
end
errorFD(:,3)    = eigFD  - eigDFPT;
errorACPsp(:,3)   = eigACPsp - eigDFPT;

