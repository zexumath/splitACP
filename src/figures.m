%% eigenvalues / gap
Ne = opt1D.Ne;
numocc = 1:Ne;
close all
plot(numocc,DNew(numocc),'bo')
set(gca,'FontSize',20)
hold on
numunocc = Ne+1:Ne+opt1D.NeExtra;
plot(numunocc,DNew(numunocc),'rd')
axis tight
legend('occupied','unoccupied')
xlabel('index($i$)','Interpreter','latex','FontSize', 20)
ylabel('$\varepsilon_i$','Interpreter','latex','FontSize', 20)
wysiwyg;

%% electron density

xgrid = (0:opt1D.NsGrid-1)*opt1D.hs;
plotregion = 10/opt1D.atom.NsCell * opt1D.NsGrid;
close all
plot(xgrid(1:plotregion),rhoNew(1:plotregion),'b-','Linewidth',1);
set(gca,'FontSize',20)
xlabel('$x$','Interpreter','latex','FontSize', 20)
ylabel('$\rho(x)$','Interpreter','latex','FontSize', 20)
wysiwyg;
% axis tight

%% \delta \rho
close all
plot(rptGrid,gp(:,10),'k--d')
set(gca,'FontSize',20)
hold on;
plot(rptGrid,drhodR(:,10),'b-*')
plot(rptGrid,drhodRCSsp(:,10),'r-o')
plot(rptGrid,drhodRCS(:,10),'g-o')
% plot(rptGrid,(rhoPerturbNew(:,5)-rho)/delta,'k--')
legend({'$\frac{\partial V_I(x-R_I)}{\partial R_I}$',...
'$\frac{\delta \rho}{\delta R_I}$ DFPT ',...
'$\frac{\delta \rho}{\delta R_I}$ split-ACP'},'Interpreter','latex')
axis tight
xlabel('$x$','Interpreter','latex','FontSize', 20);
ylabel('Perturbation and Response','Interpreter','latex','FontSize', 20);
wysiwyg

%% \delta \rho diff
close all
plot(rptGrid,-drhodR(:,10) + drhodRCSsp(:,10),'k-o')
set(gca,'FontSize',20)
hold on;
plot(rptGrid,-drhodR(:,10) + drhodRCS(:,10),'r-o')
% plot(rptGrid,(rhoPerturbNew(:,5)-rho)/delta,'k--')
legend({'Error split-ACP',...
'Error ACP'},'Interpreter','latex')
axis tight
xlabel('$x$','Interpreter','latex','FontSize', 20);
ylabel('Difference on Response','Interpreter','latex','FontSize', 20);
wysiwyg

%% phonon spectrum
% eigACPsp = eig(HessACPsp);
% eigFD  = eig(HessFD);
% eigDFPT = eig(HessDFPT);

close all;sig = 0.002;
specplot(eigDFPT,sig,[-10*sig,eigDFPT(end)+10*sig],'b-',1.5);
set(gca,'FontSize',20)
hold on
specplot(eigACPsp,sig,[-10*sig,eigDFPT(end)+10*sig],'ro',1);
specplot(eigFD,sig,[-10*sig,eigDFPT(end)+10*sig],'k--',1.5);
legend('DFPT','split ACP','FD');
xlabel('phonon frequency','FontSize', 20)
ylabel('$\varrho_D$','Interpreter','latex','FontSize', 20);
wysiwyg;

%% Ntot ratio
close all;
errscale = zeros(8,3);
DFPTtimescale = zeros(8,1);
ACPtimescale = zeros(8,1);

ACPtimeregularscale  = zeros(8,1);
ACPtimesingularscale = zeros(8,1);

idx = 0;
for testratio = 0.9:-0.1:0.2
    idx = idx +1;
    resfilename = sprintf('../res//res052017/metal1Dratio%1.1s.mat',testratio);
    load(resfilename);
    
    errscale(idx,:)     = max(abs(errorACP));
    DFPTtimescale(idx)  = DFPTtime.total;
    ACPtimescale(idx)   = sACPtime.total;
    ACPtimeregularscale(idx)   = sACPtime.regular;
    ACPtimesingularscale(idx)   = sACPtime.singular;
end

figure(1)
ratio = 1.9:-0.1:1.2 ;
semilogy(ratio,errscale(:,1),'-o');
set(gca,'FontSize',20)
hold on
semilogy(ratio,errscale(:,2),'-o');
semilogy(ratio,errscale(:,3),'-o');
legend('bare response','total response','phonon spectrum');
xlabel('$N^{(s)}_{\mathrm{cut}} / N_e$','Interpreter','latex','FontSize', 20)
ylabel('Error','FontSize', 20);
wysiwyg;

% figure(2)
% semilogy(ratio,DFPTtimescale,'-o');
% set(gca,'FontSize',20)
% hold on
% semilogy(ratio,ACPtimescale,'-o');
% semilogy(ratio,ACPtimesingularscale,'-o');
% semilogy(ratio,ACPtimeregularscale,'-o');
% legend('DFPTtime','ACPtime','ACPsingular','ACPregular');
% xlabel('Nextra / Ne ratio','FontSize', 20)
% ylabel('Time','FontSize', 20);
% wysiwyg;

%% Time scale
FD = [];
DFPT = [];
ACP = [];
spACP = [];
for Natom = 10:10:80
    resfilename = sprintf('../res/res060717/metal1DT0-%d.mat',Natom);
    load(resfilename);
    FD(end+1) = FDtime;
    DFPT(end+1) = DFPTtime.total;
    spACP(end+1) = spACPtime.total;
    ACP(end+1) = ACPtime.total;
end

Natom = 10:10:80;
close all
figure(1)
loglog(Natom,DFPT,'b-o');
set(gca,'FontSize',20)
hold on;
loglog(Natom,FD,'k-^');
loglog(Natom,spACP,'m-*');
loglog(Natom,ACP,'r-d');
legend('DFPT','FD','split ACP','ACP','Location','NorthWest')
axis tight
xlabel('System size: # of Atom','FontSize', 20)
ylabel('Time (s) ','FontSize', 20);
wysiwyg