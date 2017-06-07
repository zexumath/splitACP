function [zshift, zweight] = getpole0( Npole, Gap, DeltaE, mu )
% getpole0.m generates the poles and weights for the calculation of
% electronic density when zero temperature is calculated.
% 
% Only lower half of the circle is included.

% T = 0;

% K2au = 3.166815d-6;
% au2K = 315774.67;
% beta = au2K/(T);
beta = 0;

Npolehalf = Npole/2;
M = DeltaE;
% mshift = (pi/beta)^2;
mshift = 0;
m2 = mshift+(Gap)^2;
M2 = M^2;
k = (sqrt(M2/m2)-1)/(sqrt(M2/m2)+1);
L = -log(k)/pi;
[K,Kp] = ellipkkp(L);

t = .5i*Kp - K + (.5:Npolehalf)*2*K/Npolehalf;  %t = .9999999999i*Kp - K + (.5:Npolehalf)*2*K/Npolehalf;
[u cn dn] = ellipjc(t,L);
z = sqrt(m2*M2)*((1/k+u)./(1/k-u));
dzdt = cn.*dn./(1/k-u).^2;

zsqrt = sqrt(z-mshift);


zweight = zeros(Npole, 1);
zshift = zeros(Npole, 1);

% fd = @(z) 1.0./(1+exp(beta*z));  % NOTE: No spin degeneracy here!
fd = @(z) (1-sign(real(z)))/2;

% From Eq. (2.10) in 
%   Lin, Lu, Ying and E, " Pole-Based Approximation of the Fermi-Dirac
%   Function",  Chin. Ann. Math. 30B, 729, 2009
%
for j = 1 : Npolehalf
  zshift(j) = mu + zsqrt(j);
  zshift(j+Npolehalf) = mu - zsqrt(j);
  % Old version which does not include the identity in the formulation.
  %
  % zweight(j) = ...
    % 2*K*sqrt(m2*M2)/(k*pi*Npolehalf) / zsqrt(j) * ...
    % dzdt(j) * (-tanh(beta*zsqrt(j)/2));
  % zweight(j+Npolehalf) = ...
    % 2*K*sqrt(m2*M2)/(k*pi*Npolehalf) / (-zsqrt(j)) * ...
    % dzdt(j) * (-tanh(beta*(-zsqrt(j))/2));
  %
  % New version which takes into account the identity.
  zweight(j) = ...
    2*K*sqrt(m2*M2)/(k*pi*Npolehalf) / zsqrt(j) * ...
    dzdt(j) * fd(zsqrt(j));
  zweight(j+Npolehalf) = ...
    2*K*sqrt(m2*M2)/(k*pi*Npolehalf) / (-zsqrt(j)) * ...
    dzdt(j) * fd(-zsqrt(j));
end