function atom = initAtom1D(N)

if nargin ==0
    N = 30;
end

atom.NsCell = N;
atom.LsCell = 2.4;

atom.Z = ones(atom.NsCell,1);
atom.sigma = atom.LsCell * 0.125 * ones(atom.NsCell,1);

ion1 = atom.LsCell *[0.5];

iontmp = (0:atom.NsCell-1)'*atom.LsCell;

ion = zeros(atom.NsCell,1);
ion(1:1:end,:) = iontmp+repmat(ion1,atom.NsCell,1);
atom.ion = ion;
% Add some randomness
% atom.ion = atom.ion + 0.05*randn(size(atom.ion));

end
