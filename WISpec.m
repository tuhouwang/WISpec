clc
clear
close all
tic;
N  = 30;  % Number of truncated order
M  = 3000;% Number of discrete wavenumber
f  = 20;  % Frequency of sound source
zs = 36;  % Depth of sound source
H  = 100; % Depth of the water column
dz = 0.5;
dr = 1;
rmax = 3000;
c    = 1500; % Speed of sound in water coulmn
rho  = 1.0;  % Density of water coulmn

% Differentiation matrix and Chebyshev-Gauss-Lobatto points
D   = DerivationMatrix(N+1);
x   = cos((0 : N) * pi / N)';

z1  = 0  : dz : zs;
z2  = zs : dz : H;
z   = 0  : dz : H;

x1  = -2 / zs * z1 + 1;
x2  = -2 /(H - zs) * z2 + (H + zs) / (H - zs);

k0 = 2 * pi * f / c;
kr = linspace(0, 2 * k0, M);
eps= 3 * k0 / pi / (M - 1) / log10(exp(1.0));
kr = kr - 1i * eps;

psi = zeros(length(z),M);
for m = 1 : M

    A = 4.0 / zs ^ 2 * D * D + ...
        ConvolutionMatrix(ChebTransFFT(N,(k0^2-kr(m)^2)*ones(N+1,1)));

    B = 4.0 / (H - zs) ^ 2 * D * D + ...
        ConvolutionMatrix(ChebTransFFT(N,(k0^2-kr(m)^2)*ones(N+1,1)));
    
    U = zeros(2*N+2);
    
    U(1 :N-1,  1    :N-1)   = A(1:N-1, 1:N-1);
    U(1 :N-1,  2*N-1:2*N)   = A(1:N-1, N:N+1);
    U(N:2*N-2, N    :2*N-2) = B(1:N-1, 1:N-1);
    U(N:2*N-2, 2*N+1:2*N+2) = B(1:N-1, N:N+1);

    %upper boundary
    U(2*N-1, 1    :N-1) = 1.0;
    U(2*N-1, 2*N-1:2*N) = 1.0;
    
    %lower boundary
    %perfectly free/rigid
    Lb = (-1.0) .^ (0 : N);% * D;
    Lb = Lb - Lb * D * 2.0 / (H - zs) * 1.5 / sqrt(kr(m)^2 - (2 * pi * f / 2000)^2);
    
    U(2*N+2, N    :2*N-2) = Lb(1:N-1);
    U(2*N+2, 2*N+1:2*N+2) = Lb(N:N+1);
    
    
    
    
    %sound pressure is continuous
    U(2*N, 1    :N-1  ) = (-1.0).^(0:N-2);
    U(2*N, N    :2*N-2) =  -1.0;
    U(2*N, 2*N-1:2*N  ) = (-1.0).^(N-1:N);
    U(2*N, 2*N+1:2*N+2) =  -1.0;
    
    %displacement potential function is continuous 
    Pu = -2 / zs * ((-1.0).^(0 : N)) * D;
    Pd =  2 /(H - zs) * ones(1, N+1) * D;

    %second interface boundary
    U(2*N+1, 1    :N-1  ) = Pu(1 : N-1);
    U(2*N+1, N    :2*N-2) = Pd(1 : N-1);
    U(2*N+1, 2*N-1:2*N  ) = Pu(N : N+1);
    U(2*N+1, 2*N+1:2*N+2) = Pd(N : N+1);
    
    R = zeros(2*N+2,1);
    R(2*N+1) = - 0.5 / pi;
    
    vec  = U \ R;
    vec1 = [vec(1 :   N-1); vec(2*N-1:2*N  )];
    vec2 = [vec(N : 2*N-2); vec(2*N+1:2*N+2)];
    
    psi1 = InvChebTrans(vec1, x1);
    psi2 = InvChebTrans(vec2, x2);
    
    psi(:, m) = [psi1(1:end-1); psi2];

end

r   = dr : dr : rmax;
phi = zeros(length(z),length(r));

for ir = 1 : length(r)
    for iz = 1 : length(z)
        kenel = psi(iz, :) .* besselj(0, real(kr) * r(ir)) .* real(kr);
        phi(iz, ir) = trapz(kr, kenel);
    end
end

phi0 = exp(1i * k0) / 4 / pi; 
tl   = - 20 * log10(abs(phi / phi0));

%-------------------plot------------------------
% figure; hold on;
% plot(kr, abs(psi(37,:)), 'r-',  'LineWidth', 3);
% plot(kr, abs(psi(47,:)), 'b--', 'LineWidth', 2);
% xlabel('kr (1/m)'); ylabel('Magnitude');
% legend('z=36 m','z=46 m','box','off');
% set(gca,'Position',[0.1,0.15,0.8,0.75]);
% set(gcf,'Position',[400,200,1200,600]);
% set(gca,'FontSize',16,'LineWidth',2,'box','on');
% axis([0 0.2 0 100]);


pcolor(r, z, tl);
colormap(flipud(jet));caxis([40 70]); 
shading flat; view(0, -90);
xlabel('Range (m)'); ylabel('Depth (m)');
colorbar('YDir','Reverse','FontSize', 16);
set(gca,'Position',[0.1,0.15,0.75,0.75],'FontSize',16);
set(gcf,'Position',[100,100,800,800]); 

toc;
