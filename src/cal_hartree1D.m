function V = cal_hartree1D(rho, opt1D)

Ls = opt1D.Ls;
kappa = opt1D.kappa;
epsl0 = opt1D.epsl0;
GptGrid = opt1D.GptGrid;

mult = GptGrid.^2 + kappa^2;
mult = [0;1./mult(2:end)];
mult = spdiags(mult,0,opt1D.NsGrid,opt1D.NsGrid);
V  = real(ifft(4*pi*mult*fft(rho)));
V  = V / epsl0;
end
