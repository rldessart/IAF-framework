clear; clc; close all;

addpath('../scripts')
addpath('../models')

disp('Running Figure 1: Magnetic Contact')

% ============================================
% MODEL 1 — CONTACT (sem ruído)
% ============================================
[T, dx] = synthetic_models(3);

% ============================================
% COMPUTE ATTRIBUTES
% ============================================
[ISA, pTheta, IAF] = compute_all_attributes(T, dx);

% ============================================
% NORMALIZATION
% ============================================
ISA    = ISA    / max(abs(ISA(:)));
pTheta = pTheta / max(abs(pTheta(:)));
IAF    = IAF    / max(abs(IAF(:)));

% ============================================
% PLOT
% ============================================
figure

subplot(2,2,1)
imagesc(T)
caxis([-1 1])
title('(a) Magnetic Field')
axis equal tight
colorbar

subplot(2,2,2)
imagesc(ISA)
caxis([-1 1])
title('(b) ISA')
axis equal tight
colorbar

subplot(2,2,3)
imagesc(pTheta)
caxis([-1 1])
title('(c) p-Theta')
axis equal tight
colorbar

subplot(2,2,4)
imagesc(IAF)
caxis([-1 1])
title('(d) IAF')
axis equal tight
colorbar
