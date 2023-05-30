% Ocean acoustics.

% Copyright (C) 2023 Houwang Tu
% ------------------------------------------------------------------------|
% This program is free software: you can redistribute it and/or modify it |
% under the terms of the GNU General Public License as published by the   |
% Free Software Foundation, either version 3 of the License, or (at your  |
% option) any later version.                                              |
%                                                                         |
% This code is distributed in the hope that it will be useful, but without|
% any warranty; without even the implied warranty of merchantability or   |
% fitness for a particular purpose. See the GNU General Public License for|
% more details.                                                           |
%                                                                         |
% You should have received a copy of the GNU General Public License along |
% with this program. If not, see <http://www.gnu.org/licenses/>.          |
%                                                                         |
% The program WISpec.m computes the acoustic field in arbitrary           |
% horizontally stratified media using the Chebyshev--Tau spectral method. |
% The method is described in the article (H. Tu, Y. Wang, W. Liu et al.,  |
% A spectral method for the depth-separated solution of a wavenumber      |
% integration model for horizontally stratified fluid acoustic waveguides |
% Physics of Fluids, 2023(35), 057127, https://doi.org/10.1063/5.0150221).|
% under the supervision of Prof. Yongxian Wang, National University of    |
% Defense Technology, China.                                              |                                                  
% ------------------------------------------------------------------------|

clc
clear
close all
tic;

[casename, Src, Layers, Nl, kmax, M, freq, zs, dz, rmax, dr, tlmin, tlmax,...
 dep, c, rho, alpha, Lb, ch, rhoh, alphah] = ReadEnvParameter('input.txt');

% Get the z and rho of the final resolution.
[z, rhoi] = FinalResolute(dep, dz, rho, Layers);

% Add a virtual interface to the depth of the sound source.
[dep, c, rho, alpha, Layers, Nl, R, s] = VirtualInterface(dep, c, rho, ...
 alpha, zs, Layers, Nl);

[c, rho, alpha] = ChebInterpolation(dep, c, rho, alpha, Layers, Nl);

[r, k, kh, w] = ChebInitialization(Layers, freq, rmax, dr, c, alpha, ch, alphah);

% Set up complex contour.
k0  = max(real(k{1}));
kr  = linspace(0, kmax * k0, M);
eps = 1.5 * kmax * k0 / pi / (M - 1) / log10(exp(1.0));
kr  = kr - 1i * eps;

%------------------------------Depth equation------------------------------
Vec  = cell(Layers, M);
Vec2 = cell(Layers, 1);
for m = 1 : M
    m
    Vec(:, m) = ChebDepthSolution(Lb, Nl, Layers, dep, k, rho, kh, rhoh, kr(m), R);
end

for i = 1 : Layers
    Vec2(i) = {cell2mat(Vec(i,:))};
end

psi = KernelFunc(Vec2, dz, dep, Layers);

% Plot(kr, psi(201,:));
% Pcolor(real(kr), z, abs(psi), casename, 0, 5);colormap(jet);

toc;
%--------------------------Wavenumber Integration--------------------------
if(Src == 'P')
    % Ponit source
    phi  = kmax * k0 / (M - 1) * psi * diag(kr) * besselj(0, kr.' * r); 
    phi0 = 0.25 * exp(1i * k{s}(end)) / pi;
else    
    % Line Source
    phi  = 2 * kmax * k0 / (M - 1) * psi * cos(kr.' * r);
    phi0 = 0.25 * 1i * besselh(0, 1, 1);
%     phi0 = 0.25 * 1i * besselh(0, 1, k{s}(end));
end

% Sound pressure from displacement potential function.
phi  = w ^ 2 * diag(rhoi)  * phi;
phi0 = w ^ 2 * rho{s}(end) * phi0;
tl   = - 20  * log10(abs(phi / phi0));
Pcolor(r, z, tl, casename, tlmin, tlmax);

toc;
