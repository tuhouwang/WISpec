clc
clear
% close all
tic;

[casename, Src, Layers, Ns, kmax, M, freq, zs, dz, rmax, dr, tlmin, tlmax,...
 dep, c, rho, alpha, Lb, ch, rhoh, alphah] = ReadEnvParameter('input_line_pekeris.txt');

%Get the z and rho of the final resolution
[z, rhoi] = FinalResolute(dep, dz, rho, Layers);

%Add a virtual interface to the depth of the sound source
[dep, c, rho, alpha, Layers, Ns, R, s] = VirtualInterface(dep, c, rho, alpha, zs, Layers, Ns);

[c, rho, alpha] = ChebInterpolation(dep, c, rho, alpha, Layers, Ns);

[r, k, kh, w] = ChebInitialization(Layers, freq, rmax, dr, c, alpha, ch, alphah);

%Set up complex contour
k0  = max(real(k{1}));
kr  = linspace(0, kmax * k0, M);
eps = 1.5 * kmax * k0 / pi / (M - 1) / log10(exp(1.0));
kr  = kr - 1i * eps;

%------------------------------Depth equation------------------------------
psi = zeros(length(z), M);
parfor m = 1 : M
    m
    Vec = ChebDepthSolution(Lb, Ns, Layers, dep, k, rho, kh, rhoh, kr(m), R);

    psi(:, m) = KernelFunc(Vec, dz, dep, Layers);

end

% Plot(kr, psi(201,:));
toc;
%--------------------------Wavenumber Integration--------------------------
if(Src == 'P')
    % Ponit source
    phi  = kmax * k0 / (M - 1) * psi * diag(kr) * besselj(0, kr.' * r); 
    phi0 = 0.25 * exp(1i * k{s}(end)) / pi;
else    
    % Line Source
%     phi  = kmax * k0 / (M - 1) * psi * exp(1i * kr.' * r);
    phi  = 2 * kmax * k0 / (M - 1) * psi * cos(kr.' * r);
    phi0 = 0.25 * 1i * besselh(0, 1, 1);
%     phi0 = 0.25 * 1i * besselh(0, 1, k{s}(end));
end

%Sound pressure from displacement potential function
phi  = w ^ 2 * diag(rhoi)  * phi;
phi0 = w ^ 2 * rho{s}(end) * phi0;
tl   = - 20  * log10(abs(phi / phi0));

% Plot(r, tl(185,:));
Pcolor(r, z, tl, casename, tlmin, tlmax);

toc;
