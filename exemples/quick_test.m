% QUICK TEST SCRIPT (IAF FRAMEWORK)
% This script generates synthetic models and computes
% ISA, p-Theta, and IAF attributes.

addpath('../scripts')
addpath('../models')

clear; close all; clc;

model_type = 2;   % 1=body, 2=dike, 3=contact

disp('Running synthetic model...')

% ============================================
% GENERATE MODEL
% ============================================
[T, dx] = synthetic_models(model_type);

figure
imagesc(T)
title('Synthetic Model')
axis equal tight
colorbar

% ============================================
% COMPUTE ATTRIBUTES
% ============================================
[ISA, pTheta, IAF] = compute_all_attributes(T, dx);

% ============================================
% NORMALIZATION (visual consistency)
% ============================================
ISA    = ISA    / max(abs(ISA(:)));
pTheta = pTheta / max(abs(pTheta(:)));
IAF    = IAF    / max(abs(IAF(:)));

% ============================================
% DISPLAY RESULTS
% ============================================
figure

subplot(1,3,1)
imagesc(ISA)
title('ISA'); axis equal tight; colorbar

subplot(1,3,2)
imagesc(pTheta)
title('p-Theta'); axis equal tight; colorbar

subplot(1,3,3)
imagesc(IAF)
title('IAF'); axis equal tight; colorbar
