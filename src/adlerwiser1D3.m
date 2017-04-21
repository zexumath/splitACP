function y = adlerwiser1D3(gfunc, HMat, DNew, Vocc, QProj, opt1D,index)
global tmp0;
global RELRES;
Ne = opt1D.Ne;
N = opt1D.NsGrid;
MTeter = opt1D.MTeter;

b=Vocc;
y=zeros(N,1);
relrestmp = zeros(Ne,1);
for i=1:Ne
    b(:,i) = gfunc(b(:,i));
    b(:,i) = QProj(b(:,i));
    AFun = @(x) QProj( DNew(i) * x - HMat(x) ) ;
    [tmp,~,relrestmp(i),ite,~] = minres(AFun,b(:,i),1e-6,50,MTeter,[],tmp0(:,i,index));
    fprintf('ite = %d, relrestmp(%d) = %f\n',ite, i,relrestmp(i));
    tmp = QProj(tmp);
    tmp0(:,i,index) = tmp;
    y = y+2*tmp.*Vocc(:,i);
end
RELRES(:,end+1) = relrestmp;