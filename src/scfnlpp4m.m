function [Energy, force, rho] = scfnlpp4m(opt1D)
% This code wrap up scf procedure, returning total energy(scalar), force on
% each atom(NsCell2D vector), and electron density(NsGrid vector).
% opt1D = initOptnlpp(opt1D);

[rhoNew, ~, DNew,  VNew,occ,~,~, ~, ~] = metaltest1(opt1D);
id      = find(occ>1e-6);
Nocc    = length(id);
atom    = opt1D.atom;
hs      = opt1D.hs;
NsCell  = atom.NsCell;
Ne      = opt1D.Ne;

mTotal0 = opt1D.mTotal0;
mx      = opt1D.mx;
gpx     = cal_hartree1D(mx,opt1D);
b       = opt1D.nlpp;
bx      = opt1D.nlppx;

rho     = rhoNew;
Eband   = sum(DNew.*occ);

temperature = opt1D.temperature;
Tbeta       = temperature*3.166815d-6;

Entropy     = getEntropy(occ,Tbeta);
Ecorrection = -0.5*sum((rho-mTotal0).*cal_hartree1D(rho+mTotal0,...
                opt1D))*hs;
Energy      = Ecorrection + Eband + Entropy;
force       = zeros(NsCell,1);
% Force is negative gradient

for i = 1:NsCell
    for j = 1:Nocc
        force(i) = force(i) + 2 * occ(j) * sum( b(:,i) .* VNew(:,j) ) ...
            * sum(bx(:,i).*VNew(:,j)) * hs^2;
    end 
    force(i) = force(i) - sum(gpx(:,i).*(mTotal0+rho)) * hs;
end

end
