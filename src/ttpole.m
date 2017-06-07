% TTPOLE is an example to test that the pole expansion works for Fermi-Dirac
% distribution.
%
% This tests the electron density for complex Hermitian matrices
%
% Lin Lin
% Lastest revision: 09/09/2016

% The unit of energy is hartree!

opt1D = initOptnlpp4m();
opt1D.temperature = 0;
[rhoNew, VpotNew, DNew, VNew, occ, efermi, INDNew, HMatNew, iter,update] = metaltest1(opt1D);

Npole    = 40;
T        = opt1D.temperature;
if T>0
        Gap = 0.0;
    else
        Gap = abs(DNew(opt1D.Ne) - efermi);
    end
DeltaE  = max(2.0, DNew(end) - DNew(1));

mu       = efermi;

disp('Fermi-Dirac');

identity = eye(opt1D.NsGrid);
H = HMatNew(identity);
H = (H + H')/2;
Ev = DNew;

hs = opt1D.hs;
Hc      = VNew * diag(DNew) * VNew' * hs;
Hc      = 1/2 * (Hc + Hc');

[zshift, zweight] = getpole(Npole, T, Gap, DeltaE, mu);

RhoExact = VNew*diag(occ)*VNew'*opt1D.hs;

RhoPole = zeros(size(RhoExact));

for l = 1 : Npole
  Gl = VNew / diag(DNew - zshift(l)) *VNew' * hs;
  % Does not work for complex Hermitian matrices.
  % The number of electrons is correct, but the electron density is not
  % RhoPole = RhoPole + imag(zweight(l)*Gl);
  
  % Correct version: sum and correct
  RhoPole = RhoPole + zweight(l)*Gl;
end
% Correction step
RhoPole = 1/(2i)*(RhoPole-RhoPole');

fprintf('Number of electrons (exact) = %25.15f\n', sum(diag(RhoExact)) );
fprintf('Number of electrons (pole)  = %25.15f\n', sum(diag(RhoPole)) );
fprintf('||RhoExact - RhoPole||_2    = %25.15f\n', ...
  norm(RhoExact-RhoPole));
