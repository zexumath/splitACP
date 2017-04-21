function [RHSsel,sel,CHI0gCS, CHI0gCStmp,W,PI,CStime,Pivotsel,indrand,rphase,eqnres,dpsib,coeffMatsel]...
    = compresschi01Dtest(opt1D,gpfunc,Vocc,occ,DNew,HMatNew,QProj,NchebNodes,CStime,...
    Pivotsel,indrand,rphase,select,eqnres)
% This function apply the compressed CHI0 to the right-hand-side matrix g.
% Script is used to build the adaptively compressed CHI.

atom = opt1D.atom;
NsGrid = opt1D.NsGrid;
NsCell = atom.NsCell;
Nocc = size(Vocc,2);

% disp('Computing fft on each column of G...');
% Generate the columns of G' and compute fft.
% Sample with a preselected indrand and store into Gsel, which is
% already transposed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% time of column selection
tic 
preselect = min(opt1D.select * Nocc, Nocc*NsCell);
if nargin <=10
    indrand = randperm(Nocc*NsCell,preselect);
end
if nargin<=11
    rphase = exp(1i*2*pi*rand(Nocc*NsCell,1));
end
Gsel = zeros(preselect,NsGrid);

G = zeros(NsGrid,Nocc*NsCell);
for i = 1:Nocc
    for j = 1:NsCell
        G(:,(i-1)*NsCell+j) = gpfunc(Vocc(:,i),j);
    end
end


for rowindex = 1:NsGrid
    Gtmp = fft(G(rowindex,:)'.*rphase);
    Gsel(:,rowindex) = Gtmp(indrand);
end
tau = opt1D.tau;



% disp('Compute QR factorization of Gsel...');
% tic
if nargin <=9
    [Qsel,Rsel,Pivotsel] =qr(Gsel,'vector');
else
    [Qsel, Rsel] = qr(Gsel(:,Pivotsel),'vector');
end
if nargin <=12
    ind = find(abs(diag(Rsel))>tau*abs(Rsel(1,1)));
    select = length(ind);
end
% select = 6*Ne;
% select = opt1D.forcedselect*Ne;
% toc

% Compute CHI0gCS
% disp('Compute CHI0gCS with results of QR given...');
% tic
RHSsel = [eye(select),...
    Rsel(1:select,1:select)\Rsel(1:select,select+1:end)];
% norm(RHSsel,'fro')
RHSsel(:,Pivotsel) = RHSsel;
RHSsel = real(RHSsel');
RHSselP = QProj(RHSsel);

coeffMatsel = G(Pivotsel(1:select),:);
sel = Pivotsel(1:select);
if nargout <=2
    return
end

CStime.colsel =  CStime.colsel+ toc; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% time of solving equations
tic

% 
% cl = DNew(1) - 1;
% cr = (DNew(Ne) + DNew(Ne+1))/2;

cl = DNew(1);
cr = DNew(Nocc);

% NchebNodes = 20;

nodes = cheb_nodes(NchebNodes,[cl,cr]);

if nargin <=12
    eqnres = zeros(NsGrid,select,NchebNodes);
end

MTeter = opt1D.MTeter;
% global RELRESCS;
relrestmp = zeros(select,NchebNodes);
for j = 1:NchebNodes
    %     HMatInvCheb = pinv(nodes(j) * IMat - HMat);
    AFUN = @(u) QProj( nodes(j)*u - HMatNew(u) );
    %     vChebtmp(:,:,j) = HMatInvCheb*RHSP;
    for i = 1:select
        [tmp,~,relrestmp(i,j),ite,RESVEC] = minres(AFUN,RHSselP(:,i),1e-6,50,...
            MTeter,[],eqnres(:,i,j));
        eqnres(:,i,j) = QProj(tmp);
        fprintf('ite = %d, relrestmp(%d,%d) = %f\n',ite, i,j,relrestmp(i,j));
%         for i=1:length(RESVEC)
%             fprintf('%d : RESVEC(i)= %g\n',i,RESVEC(i));
%         end
    end
end
% RELRESCS(:,:,end+1) = relrestmp;

CHI0gCStmp = zeros(NsGrid,select,Nocc);
CStime.soleqn = CStime.soleqn + toc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% time of reconstructing solutions
tic

for i = 1:Nocc
    CHI0gCStmp(:,:,i) = lagpoly(DNew(i),nodes,eqnres);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% time of getting chi0g.
%     stop

W = zeros(NsGrid,select);

for i = 1:Nocc
    W = W + 2*occ(i)*spdiags(conj(Vocc(:,i)),0,NsGrid,NsGrid)...
        * CHI0gCStmp(:,:,i)...
        * spdiags(Vocc(Pivotsel(1:select),i),0,select,select);
end



PI = @(A) A(Pivotsel(1:select),:);
% CHI0 = W*PI;
%     stop
CHI0gCS = zeros(NsGrid,NsCell);
dpsib = zeros(NsGrid,Nocc,NsCell);
for i = 1:Nocc
    tmp = CHI0gCStmp(:,:,i)* coeffMatsel(:,(i-1)*NsCell+1 : i*NsCell);
    dpsib(:,i,:) = tmp;
    CHI0gCS = CHI0gCS + 2*diag(Vocc(:,i))*occ(i) * tmp;
end
% toc
CStime.reconstruct = CStime.reconstruct + toc;
end
