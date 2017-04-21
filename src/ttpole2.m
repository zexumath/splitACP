
opt1D   = initOptnlpp4m();
[rhoNew, VpotNew, DNew, VNew, occ, ...
    efermi, INDNew, HMatNew, iter,update] = metaltest1(opt1D);

Npole   = 60;
T       = opt1D.temperature;
Gap     = 0.0;
DeltaE  = 2.0;
mu      = efermi;
hs      = opt1D.hs;
Ls      = opt1D.Ls;
NsGrid  = opt1D.NsGrid;
itertol = opt1D.itertol;
Ne      = opt1D.Ne;

[zshift, zweight] = getpole(Npole, T, Gap, DeltaE, mu);

Hc      = VNew * diag(DNew) * VNew' * hs;
Hc      = 1/2 * (Hc + Hc');
Nt      = length(DNew);
m0      = opt1D.m0;
mTotal0 = opt1D.mTotal0;
mx      = opt1D.mx;
mxx     = opt1D.mxx;


gpx     = cal_hartree1D(mx , opt1D );
gpxx    = cal_hartree1D(mxx, opt1D );
gp      = gpx;

b       = opt1D.nlpp;
bx      = opt1D.nlppx;
bxx     = opt1D.nlppxx;

gpfunc  = @(x,j) gp(:,j).*x ...
    - b(:,j) * (bx(:,j)'*x)*hs...
    - bx(:,j) * (b(:,j)'*x)*hs;



res = zeros(NsGrid, NsGrid);
identity = eye(NsGrid, NsGrid);
for l = 1:Npole
    Gl = VNew / diag(DNew - zshift(l)) *VNew' * hs;
    tmp1 = zeros(NsGrid, NsGrid);
    for i = 1:Nt
        tmp2 = gpfunc(VNew(:,i),1);
        tmp1 = tmp1 + Gl * tmp2 / (zshift(l) - DNew(i)) * VNew(:,i)';
    end
    res = res + zweight(l) * tmp1;
end
res = 1/(2i)*(res-res');


% perfect match with drhodRb1 from metalTest3.m


