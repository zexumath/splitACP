function [rhoNew, VpotNew, DNew, VNew, occ,efermi, INDNew, HMatNew, iter,update] = metaltest1(opt1D)
% This script tests adding nlpp term

atom = opt1D.atom;
hs = opt1D.hs;
Ls = opt1D.Ls;
Ne  = opt1D.Ne;

NsGrid = opt1D.NsGrid;
itertol = opt1D.itertol;
% NsCell = atom.NsCell;
% NsCell2D = atom.NsCell2D;
% Ne = sum(atom.Z);

mTotal0 = opt1D.mTotal0;


temperature = opt1D.temperature;

Tbeta       = temperature*3.166815d-6;

Vpot0 = cal_hartree1D(mTotal0,opt1D);
Vpot0 = Vpot0 - mean(Vpot0);

opts.issym = 1;
opts.isreal = 1;
opts.tol=1e-8;
opts.maxit=1000;
Neigs = Ne + opt1D.NeExtra;

b = opt1D.nlpp;

Vpot0Mat = spdiags(Vpot0,0,NsGrid,NsGrid);
HMat0 = @(x) opt1D.Lap(x) + Vpot0Mat*x - b* (b'*x)*hs;
% Compute teter preconditioner
MTeter = teter1(opt1D.GptGrid, NsGrid);
% FIXME  include in the opt
tolLOBPCG = 1e-8;
maxIterLOBPCG = 50;

% Replace by LOBPCG solver
% [V0,D0,flag] = eigs(HMat0, NsGrid, Neigs, 'sa',opts);
%
% assert(flag == 0);
%D0 = diag(D0);

% LOBPCG solver for solving the eigenvalue problem
if(1)
    V0 = randn(NsGrid, Neigs);
    [V0,D0,failureFlag,lambdaHistory, residualNormsHistory] = ...
        lobpcg(V0, HMat0, [], MTeter, ...
        [], tolLOBPCG, 20, 0);
end

[D0, IND0] = sort(D0); V0 = V0(:, IND0);

[occ,efermi] = getocc(D0,Ne,Tbeta);


V0 = V0/sqrt(hs);

rho0 = sum(repmat(occ',NsGrid,1).*V0.^2,2);

mixdim=10;
VpotNew = Vpot0;
rhoNew = rho0;
alpha = opt1D.alpha;
maxite = opt1D.itermaxite;
rStore = zeros(NsGrid,maxite+1);
yStore = zeros(NsGrid,maxite+1);
sStore = zeros(NsGrid,maxite+1);
sStore(:,1) = Vpot0;
update = [];
iter = 1;

while iter<=1
    VpotOld = VpotNew;
    rhoOld = rhoNew;
    FVpotOld = cal_hartree1D(mTotal0 + rhoOld,opt1D);
    rStore(:,iter) = VpotOld - FVpotOld;
    if iter>1
        yStore(:,iter) = rStore(:,iter) - rStore(:,iter-1);
    else
        yStore(:,iter) = rStore(:,iter);
    end
    
    VpotNew = VpotOld - alpha * rStore(:,iter);
    
    sStore(:,iter+1) = VpotNew - VpotOld;
    
    VpotNewMat = spdiags(VpotNew, 0, NsGrid, NsGrid);
    HMatNew = @(x) opt1D.Lap(x) + VpotNewMat*x - b* (b'*x)*hs;
    
    % [VNew, DNew,flagNew] = eigs(HMatNew,NsGrid,Neigs,'sa',opts);
    % assert(flagNew==0);
    % DNew = diag(DNew);
    
    
    % Use LOBPCG
    if(1)
        [VNew,DNew,failureFlag,lambdaHistory, residualNormsHistory] = ...
            lobpcg(V0, HMatNew, [], MTeter, ...
            [], tolLOBPCG, maxIterLOBPCG, 0);
    end
    
    [DNew, INDNew] = sort(DNew,'ascend'); VNew = VNew(:,INDNew);
    
    [occ,efermi] = getocc(DNew,Ne,Tbeta);
    
    VNew = VNew/sqrt(hs);
    rhoNew = sum(repmat(occ',NsGrid,1).*VNew.^2,2);
    
    update(end+1) = norm(VpotNew - VpotOld)./norm(VpotOld);
    
    iter = iter + 1;
%     fprintf('iter = %d, update = %2.15f\n',iter-1,update);
end


while update(end) >itertol && iter<=maxite
    VpotOld = VpotNew;
    rhoOld = rhoNew;
    FVpotOld = cal_hartree1D(mTotal0 + rhoOld,opt1D);
    
    rStore(:,iter) = VpotOld - FVpotOld;
    yStore(:,iter) = rStore(:,iter) - rStore(:,iter-1);
    
    tmp = min(iter-1,mixdim);
    
    Sk = sStore(:,iter:-1:iter-tmp+1);
    Yk = yStore(:,iter:-1:iter-tmp+1);
    tmpYkrS = Yk\rStore(:,iter);
    sStore(:,iter+1) = -alpha*rStore(:,iter) + (alpha*Yk-Sk)*tmpYkrS;
    VpotNew = VpotOld + sStore(:,iter+1);
    
    VpotNewMat = spdiags(VpotNew, 0, NsGrid, NsGrid);
    HMatNew = @(x) opt1D.Lap(x) + VpotNewMat*x - b* (b'*x)*hs;
    
    % [VNew, DNew,flagNew] = eigs(HMatNew,NsGrid,Neigs,'sa',opts);
    % assert(flagNew==0);
    % DNew = diag(DNew);
    
    
    % Use LOBPCG
    if(1)
        [VNew,DNew,failureFlag,lambdaHistory, residualNormsHistory] = ...
            lobpcg(VNew, HMatNew, [], MTeter, ...
            [], tolLOBPCG, maxIterLOBPCG, 0);
    end
    
    
    % DNew(Ne-1:Ne+1)
    [DNew, INDNew] = sort(DNew,'ascend');VNew = VNew(:,INDNew);
    [occ,efermi] = getocc(DNew,Ne,Tbeta);
    
    VNew = VNew/sqrt(hs);
    rhoNew = sum(repmat(occ',NsGrid,1).*VNew.^2,2);
    
    update(end+1) = norm(VpotNew - VpotOld)./norm(VpotOld);
    
    iter = iter + 1;
%     fprintf('iter = %d, update = %2.15f\n',iter-1,update);
    
%         figure(3)
%         plot(DNew,'-o'); grid on;
%         pause
end


% fprintf('norm of update = %3.20f\n',update);
% fprintf('energy = %3.20f\n',sum(DNew(1:Ne)));
% figure(2)
% imagesc(reshape(rhoNew,NsGridx,NsGridy));colorbar;
% pause
end
