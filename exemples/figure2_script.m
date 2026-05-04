clear; clc; close all;

addpath('../scripts')
addpath('../models')

disp('Running Figure 2: Multi-source with noise')

% ============================================
% MODEL 2 — MULTI-SOURCE + NOISE
% ============================================
[T, dx] = synthetic_models(2);

figure
imagesc(T)
title('(a) Synthetic Magnetic Field (Noisy)')
axis equal tight
colorbar

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

subplot(1,3,1)
imagesc(ISA)
caxis([-1 1])
title('(b) ISA')
axis equal tight
colorbar

subplot(1,3,2)
imagesc(pTheta)
caxis([-1 1])
title('(c) p-Theta')
axis equal tight
colorbar

subplot(1,3,3)
imagesc(IAF)
caxis([-1 1])
title('(d) IAF')
axis equal tight
colorbar
