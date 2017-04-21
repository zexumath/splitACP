function y = adlerwiser1D4m(gfunc, HMat, DNew, Vocc,occ, QtProj, opt1D,index)
global tmp0;
global RELRES;
Nocc = size(Vocc,2);
N = opt1D.NsGrid;
MTeter = opt1D.MTeter;

b=Vocc;
y=zeros(N,1);
relrestmp = zeros(Nocc,1);
for i=1:Nocc
    if occ(i)>=1e-6
        b(:,i) = gfunc(b(:,i));
        b(:,i) = QtProj(b(:,i));
        AFun = @(x) QtProj( DNew(i) * x - HMat(x) ) ;
        [tmp,~,relrestmp(i),ite,~] = minres(AFun,b(:,i),1e-6,50,MTeter,[],tmp0(:,i,index));
        fprintf('ite = %d, relrestmp(%d) = %f\n',ite, i,relrestmp(i));
        tmp = QtProj(tmp);
        tmp0(:,i,index) = tmp;
        y = y+2*tmp.*Vocc(:,i)*occ(i);
    else
        continue;
    end
end
RELRES(:,end+1) = relrestmp;