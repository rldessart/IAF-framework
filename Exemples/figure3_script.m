clear all; close all; clc;
graphics_toolkit("qt");

disp('Inclined Contact Model (30 degrees) - IAF')

% -------------------------
% GRID
% -------------------------
dx = 10;
dz = 10;

x = 0:dx:1000;
z = 0:dz:400;

[X,Z] = meshgrid(x,z);

% -------------------------
% PARÂMETROS DO CONTATO
% -------------------------
x0 = 300;
z0 = 50;
theta = 30 * pi/180;

chi1 = 0.01;
chi2 = 0.05;

% -------------------------
% MODELO GEOLÓGICO
% -------------------------
z_interface = z0 + (X - x0)*tan(theta);

model = chi1 * ones(size(X));
model(Z >= z_interface) = chi2;

% -------------------------
% CAMPO SINTÉTICO
% (mantido simples, como no original)
% -------------------------
T = atan((X-x0)./Z);

% -------------------------
% DERIVADAS
% -------------------------
[dTdz, dTdx] = gradient(T,dz,dx);

ASA = sqrt(dTdx.^2 + dTdz.^2);
epsilon = 1e-6;

% -------------------------
% ATRIBUTOS
% -------------------------
Theta = acos( dTdz ./ (ASA + epsilon) );
ISA   = atan( dTdz ./ (ASA + epsilon) );

% ✔ IAF ORIGINAL (CORRETO)
IAF = atan(ISA .* Theta);

% -------------------------
% NORMALIZAÇÃO
% -------------------------
Tn     = T / max(abs(T(:)));
ASAn   = ASA / max(abs(ASA(:)));
Thetan = Theta / max(abs(Theta(:)));
ISAn   = ISA / max(abs(ISA(:)));
IAFn   = IAF / max(abs(IAF(:)));

% -------------------------
% PLOT
% -------------------------
figure(1); clf;
set(gcf,'color','w','position',[100 100 900 800]);

subplot(2,3,1)
imagesc(x,z,model)
set(gca,'YDir','reverse')
title('(a) Geological Model')
colorbar

subplot(2,3,2)
imagesc(x,z,Tn)
set(gca,'YDir','reverse')
title('(b) Total Field')
colorbar

subplot(2,3,3)
imagesc(x,z,ASAn)
set(gca,'YDir','reverse')
title('(c) ASA')
colorbar

subplot(2,3,4)
imagesc(x,z,ISAn)
set(gca,'YDir','reverse')
title('(d) ISA')
colorbar

subplot(2,3,5)
imagesc(x,z,Thetan)
set(gca,'YDir','reverse')
title('(e) Theta')
colorbar

subplot(2,3,6)
imagesc(x,z,IAFn)
set(gca,'YDir','reverse')
title('(f) IAF')
colorbar

colormap(jet)
drawnow
