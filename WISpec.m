clc
clear
% close all
tic;

[casename, Layers, Ns, kmax, M, freq, zs, dz, rmax, dr, tlmin, tlmax,...
 dep, c, rho, alpha, Lb, ch, rhoh, alphah] = ReadEnvParameter('input_pekeris.txt');

%Add a virtual interface to the depth of the sound source
[dep, c, rho, alpha, Layers, Ns, R] = VirtualInterface(dep, c, rho, alpha, zs, Layers, Ns);

[c, rho, alpha] = ChebInterpolation(dep, c, rho, alpha, Layers, Ns);

[r, k, kh] = ChebInitialization(Layers, freq, rmax, dr, c, alpha, ch, alphah);

%Set up complex contour
k0 = max(real(k{1}));
kr = linspace(0, kmax * k0, M);
eps= 1.5 * kmax * k0 / pi / (M - 1) / log10(exp(1.0));
kr = kr - 1i * eps;

%------------------------------Depth equation------------------------------
z   = 0 : dz : dep{end}(end);
psi = zeros(length(z), M);
for m = 1 : M
    m
    Vec = ChebDepthSolution(Lb, Ns, Layers, dep, k, rho, kh, rhoh, kr(m), R);

    psi(:, m) = KernelFunc(Vec, dz, dep, Layers);

end

% Plot(kr, psi(73,:));
toc;
%--------------------------Wavenumber Integration--------------------------
%Ponit source
phi = zeros(length(z),length(r));
for ir = 1 : length(r)
    r(ir)
    for iz = 1 : length(z)
        kernel = psi(iz, :) .* besselj(0, kr * r(ir)) .* real(kr);
        phi(iz, ir) = trapz(kr, kernel);
    end
end

%Line Source

%Sound pressure from displacement potential function
phi0 = exp(1i * k0) / 4 / pi;
tl   = - 20 * log10(abs(phi / phi0));

% Plot(r, tl(73,:));
Pcolor(r,z,tl,casename,tlmin,tlmax);

toc;
