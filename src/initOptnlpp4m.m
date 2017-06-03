function opt1D = initOptnlpp4m(opt1D)

if nargin ==0
    opt1D = struct;
    atom = initAtom1D(); % {sigma Z ion LsCell NsCell}
    opt1D.atom = atom;
else
    atom = opt1D.atom;
end
hs = 0.1;
opt1D.Ne = sum(atom.Z);
opt1D.Ls = atom.LsCell * atom.NsCell;
Ls = opt1D.Ls;
opt1D.NsGrid = 2 * round(opt1D.Ls / hs /2);

opt1D.hs = opt1D.Ls / opt1D.NsGrid;

opt1D.epsl0 = 80;
opt1D.kappa = 0.1;
opt1D.itertol = 1e-8;
opt1D.itermaxite = 400;
opt1D.pcgtol = 1e-6;
opt1D.alpha = 0.5;

opt1D.temperature = 3000;

opt1D.select = 12;
opt1D.NeExtra = ceil(opt1D.Ne * 1); 

opt1D.GptGrid = (2*pi*[0:(opt1D.NsGrid)/2  -(opt1D.NsGrid)/2+1:-1]./opt1D.Ls)';

GptGrid = opt1D.GptGrid;
LapMult = 0.5 * GptGrid.^2 ;

LapMultMat = spdiags(LapMult,0,opt1D.NsGrid, opt1D.NsGrid);

opt1D.Lap = @(x) real(ifft(LapMultMat*fft(x)));

opt1D.rptGrid = (0:opt1D.NsGrid-1)'*opt1D.hs;
rptGrid = opt1D.rptGrid;

mFUN = @(x,Z,mu,sigma) -Z/sqrt(2*pi*sigma^2) * ...
    exp(-((x(:,1)-mu(:,1)).^2)/2/sigma^2);
m1FUN = @(x,Z,mu,sigma) -Z/sqrt(2*pi)/sigma ...
        * ((x(:,1)-mu(:,1))./sigma^2).*exp(-(x(:,1)-mu(:,1)).^2/2/sigma^2);

m0 = zeros(opt1D.NsGrid,opt1D.atom.NsCell);
mTotal0 = 0;
mx = zeros(opt1D.NsGrid, atom.NsCell);
mxx = zeros(opt1D.NsGrid, atom.NsCell);

for j = 1:atom.NsCell
    m0(:,j) = m0(:,j) + mFUN(rptGrid,atom.Z(j),atom.ion(j,:),atom.sigma(j));
    m0(:,j) = m0(:,j) + mFUN(rptGrid,atom.Z(j),atom.ion(j,:)+[ Ls],atom.sigma(j));
    m0(:,j) = m0(:,j) + mFUN(rptGrid,atom.Z(j),atom.ion(j,:)+[-Ls],atom.sigma(j));
    mTotal0 = mTotal0 + m0(:,j);
end


difftmp = spdiags(GptGrid,0,opt1D.NsGrid,opt1D.NsGrid);

mx  = -real(ifft(1i*difftmp*fft(m0)));
mxx = real(ifft(-1*difftmp.^2*fft(m0)));

opt1D.m0  = m0;
opt1D.mTotal0 = mTotal0;
opt1D.mx = mx;
opt1D.mxx = mxx;

opt1D.Npole = 10;
opt1D.NchebNodes = 4;
opt1D.tau = 1e-4;
opt1D.requiredeps = 1e-3;
opt1D.MTeter = teter1(opt1D.GptGrid,opt1D.NsGrid);

b = zeros(opt1D.NsGrid,atom.NsCell);
sigmab = 0.1;
for j = 1:atom.NsCell
    b(:,j) = b(:,j) + mFUN(rptGrid,-1,atom.ion(j,:)+ [Ls],sigmab);
    b(:,j) = b(:,j) + mFUN(rptGrid,-1,atom.ion(j,:)      ,sigmab);
    b(:,j) = b(:,j) + mFUN(rptGrid,-1,atom.ion(j,:)- [Ls],sigmab);
end
bx = zeros(opt1D.NsGrid,atom.NsCell);
for j = 1:atom.NsCell
    bx(:,j) = bx(:,j) + m1FUN(rptGrid,-1,atom.ion(j,:)+ [Ls],sigmab);
    bx(:,j) = bx(:,j) + m1FUN(rptGrid,-1,atom.ion(j,:)      ,sigmab);
    bx(:,j) = bx(:,j) + m1FUN(rptGrid,-1,atom.ion(j,:)- [Ls],sigmab);
end
fac = 0.1;
b = b* fac;
bx = bx* fac;
% bx1  = -real(ifft(1i*difftmp*fft(b)));
bxx = real(ifft(-1*difftmp.^2*fft(b)));


indzero = find(abs(b)<1e-8);
b(indzero) = 0;
bx(indzero) = 0;

opt1D.nlpp = sparse(b);
opt1D.nlppx = sparse(bx);
opt1D.nlppxx = bxx;

end

