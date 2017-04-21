function response1 = getResponse1(gfunc,DNew,Vocc,occ,opt1D,index)

Nt = length(occ);
NsGrid = opt1D.NsGrid;
DIM = 1;
na = opt1D.Ne; % na = Ne;

response1 = zeros(NsGrid,1);

for i = 1:Nt
    for j = 1:Nt
        if( abs(occ(i) - occ(j)) > 1e-8 )
           response1 = response1 +  gfunc(Vocc(:,i)).*conj(Vocc(:,j))
    end
end
