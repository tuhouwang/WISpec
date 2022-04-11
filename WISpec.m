clc
clear
close all
tic;

[casename, Layers, Ns, kmax, M, freq, zs, dz, rmax, dr, tlmin, tlmax,...
 dep, c, rho, alpha, Lb, ch, rhoh, alphah] = ReadEnvParameter('input_pekeris.txt');

%在声源深度上增加一个虚拟的界面
[dep, c, rho, alpha, Layers, Ns, R] = VirtualInterface(dep, c, rho, alpha, zs, Layers, Ns);

[c, rho, alpha] = ChebInterpolation(dep, c, rho, alpha, Layers, Ns);

[r, k, kh] = ChebInitialization(Layers, freq, rmax, dr, c, alpha, ch, alphah);

k0 = max(real(k{1}));
kr = linspace(0, kmax * k0, M);
eps= 3 * k0 / pi / (M - 1) / log10(exp(1.0));
kr = kr - 1i * eps;
%------------------------------Depth equation------------------------------
z = 0 : dz : dep{end}(end);
psi = zeros(length(z),M);
for m = 1 : M

    Vec = ChebDepthSolution(Lb, Ns, Layers, dep, k, rho, kh, rhoh, kr(m), R);

    psi(:, m) = KernelFunc(Vec, dz, dep, Layers);

end

Plot(kr, psi(73,:));
toc;
%--------------------------Wavenumber Integration--------------------------
phi = zeros(length(z),length(r));
for ir = 1 : length(r)
    for iz = 1 : length(z)
        kernel = psi(iz, :) .* besselj(0, real(kr) * r(ir)) .* real(kr);
        phi(iz, ir) = trapz(kr, kernel);
    end
end

phi0 = exp(1i * k0) / 4 / pi; 
tl   = - 20 * log10(abs(phi / phi0));

Pcolor(r,z,tl,casename,tlmin,tlmax);
toc;
