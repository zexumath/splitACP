function Hess = compHess(drho, gp, HMat, DNew, Vocc, QProj,opt1D)

Ne = opt1D.Ne;
N = opt1D.NsGrid;
MTeter = opt1D.MTeter;

f1 = drho;
f2 = opt1D.nlpp;
f3 = opt1D.nlppx;

dV = cal_hartree1D(drho,opt1D);

dVtot = repmat(dV,1,opt1D.atom.NsCell) + gp;

Hess = zeros(opt1D.atom.NsCell);

y=zeros(N,N);
relrestmp = zeros(Ne,1);

for i=1:Ne
    f2(:,i) = QProj(f2(:,i));
    AFun = @(x) DNew(i) * x - HMat(x);
    [tmp,~,relrestmp(i),ite,~] = minres(AFun,f1(:,i),1e-6,50,MTeter,[],0);
    tmp = QProj(tmp);
    fprintf('ite = %d, relrestmp(%d) = %f\n',ite, i,relrestmp(i));
end


end
