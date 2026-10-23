% Figure2_Multiple_Sources_With_Noise_FINAL.m
% Modelo sintetico de duas fontes para a Figura 2 do IAF.
% Dois prismas retangulares em profundidades diferentes.
% O campo total dos dois corpos e calculado antes da adicao de ruido.
%
% IMPORTANTE:
% - Script independente; nao altera os arquivos da Figura 1.
% - Os parametros geometricos do modelo sintetico estao explicitados abaixo.
% - As equacoes dos operadores foram mantidas conforme a implementacao validada.
% - Confirmar a correspondencia do protocolo de ruido e do modelo com o manuscrito.
%
% Compatibilidade: GNU Octave / MATLAB (sem toolboxes especiais).

clear; close all; clc;

%% 1. Grade de observacao
dxObs = 10;
dyObs = 10;
xObs = -500:dxObs:500;
yObs = -200:dyObs:200;
[Xobs, Yobs] = meshgrid(xObs, yObs);

%% 2. Parametros fisicos
M = 1;                       % A/m (mesmo modulo ilustrativo da Figura 1)
inclination = 45;            % graus
declination = 0;             % graus; eixo x como direcao horizontal
I = inclination*pi/180;
D = declination*pi/180;
mHat = [cos(I)*cos(D), cos(I)*sin(D), sin(I)];
fHat = mHat;                 % magnetizacao paralela ao campo indutor
mu0_over_4pi = 1e-7;

% Corpo raso: menor e a esquerda, perto da superficie
shallow.x = [-300, -120];
shallow.y = [-60, 60];
shallow.z = [70, 150];       % profundidade positiva para baixo
shallow.n = [9, 6, 6];

% Corpo profundo: mais largo e a direita, com separacao lateral clara
deep.x = [220, 500];
deep.y = [-90, 90];
deep.z = [220, 360];
deep.n = [14, 9, 8];

%% 3. Discretizar os dois corpos em celulas magnetizadas
[srcX1, srcY1, srcZ1, vol1] = prism_cells(shallow.x, shallow.y, shallow.z, shallow.n);
[srcX2, srcY2, srcZ2, vol2] = prism_cells(deep.x, deep.y, deep.z, deep.n);

srcX = [srcX1; srcX2];
srcY = [srcY1; srcY2];
srcZ = [srcZ1; srcZ2];
cellVol = [repmat(vol1, numel(srcX1), 1); repmat(vol2, numel(srcX2), 1)];

%% 4. Campo magnetico total limpo, contribucao dos dois corpos
Tclean = total_field_variable_cells(Xobs, Yobs, 0, ...
    srcX, srcY, srcZ, cellVol, M, mHat, fHat, mu0_over_4pi);

% Campos limpos nos niveis usados para a derivada vertical.
zTop = -5;
zBottom = 5;
TcleanTop = total_field_variable_cells(Xobs, Yobs, zTop, ...
    srcX, srcY, srcZ, cellVol, M, mHat, fHat, mu0_over_4pi);
TcleanBottom = total_field_variable_cells(Xobs, Yobs, zBottom, ...
    srcX, srcY, srcZ, cellVol, M, mHat, fHat, mu0_over_4pi);

%% 5. Ruido gaussiano de 15%, conforme o manuscrito
noiseLevel = 0.15;
rng(2026, 'twister'); % semente fixa para reproducibilidade
sigmaNoise = noiseLevel * std(Tclean(:), 0);
noise = sigmaNoise * randn(size(Tclean));
noiseTop = sigmaNoise * randn(size(TcleanTop));
noiseBottom = sigmaNoise * randn(size(TcleanBottom));
Tnoisy = Tclean + noise;
TnoisyTop = TcleanTop + noiseTop;
TnoisyBottom = TcleanBottom + noiseBottom;

fprintf('Campo limpo: std = %.6g nT\n', std(Tclean(:), 0));
fprintf('Ruido: nivel = %.1f%%; sigma = %.6g nT\n', ...
    100*noiseLevel, sigmaNoise);
fprintf('Celulas magnetizadas: %d\n', numel(srcX));
fprintf('Alvo raso: x=[%.0f, %.0f] m; z=[%.0f, %.0f] m\n', ...
    shallow.x(1), shallow.x(2), shallow.z(1), shallow.z(2));
fprintf('Alvo profundo: x=[%.0f, %.0f] m; z=[%.0f, %.0f] m\n', ...
    deep.x(1), deep.x(2), deep.z(1), deep.z(2));

% Campos individuais para diagnostico da contribuicao de cada alvo.
Tshallow = total_field_variable_cells(Xobs, Yobs, 0, ...
    srcX1, srcY1, srcZ1, repmat(vol1, numel(srcX1), 1), ...
    M, mHat, fHat, mu0_over_4pi);
Tdeep = total_field_variable_cells(Xobs, Yobs, 0, ...
    srcX2, srcY2, srcZ2, repmat(vol2, numel(srcX2), 1), ...
    M, mHat, fHat, mu0_over_4pi);
fprintf('Pico absoluto alvo raso: %.6g nT\n', max(abs(Tshallow(:))));
fprintf('Pico absoluto alvo profundo: %.6g nT\n', max(abs(Tdeep(:))));

%% 6. Calcular operadores no campo LIMPO e no campo RUIDOSO
% Mesmo procedimento de derivacao e mesmas equacoes nos dois casos.
[dTdx_clean, ~] = gradient(Tclean, dxObs, dyObs);
dTdz_clean = (TcleanBottom - TcleanTop) ./ (zBottom - zTop);

[dTdx_noisy, ~] = gradient(Tnoisy, dxObs, dyObs);
dTdz_noisy = (TnoisyBottom - TnoisyTop) ./ (zBottom - zTop);

%% 7. Equacoes do manuscrito aplicadas sem alteracao
epsilon = 1e-6;
[ISA_clean, Theta_clean, IAF_clean, ASA_clean] = ...
    operators_from_derivatives(dTdx_clean, dTdz_clean, epsilon);
[ISA_noisy, Theta_noisy, IAF_noisy, ASA_noisy] = ...
    operators_from_derivatives(dTdx_noisy, dTdz_noisy, epsilon);

fprintf('ASA minima campo limpo: %.6g\n', min(ASA_clean(:)));
fprintf('ASA minima campo ruidoso: %.6g\n', min(ASA_noisy(:)));
fprintf('Pontos ASA <= 10*epsilon (limpo): %d/%d\n', ...
    nnz(ASA_clean <= 10*epsilon), numel(ASA_clean));
fprintf('Pontos ASA <= 10*epsilon (ruidoso): %d/%d\n', ...
    nnz(ASA_noisy <= 10*epsilon), numel(ASA_noisy));

%% 8. Mapas: limpo e ruidoso com espacamento e escalas legiveis
figMaps = figure('Color','w','Position',[30 40 1600 900], ...
    'Name','IAF Figure 2 - Multiple sources with noise');
cleanMethods = {ISA_clean, Theta_clean, IAF_clean};
noisyMethods = {ISA_noisy, Theta_noisy, IAF_noisy};
methodNames = {'ISA','Theta','IAF'};

% Positions leave horizontal room for individual colorbars.
lefts = [0.055 0.375 0.695];
mapWidth = 0.165;
mapHeight = 0.30;
rowTop = 0.57;
rowBottom = 0.12;

for k = 1:3
    commonMin = min([cleanMethods{k}(:); noisyMethods{k}(:)]);
    commonMax = max([cleanMethods{k}(:); noisyMethods{k}(:)]);
    if k == 1
        % Symmetric range is appropriate for the signed ISA response.
        lim = max(abs([cleanMethods{k}(:); noisyMethods{k}(:)]));
        commonMin = -lim; commonMax = lim;
    end
    if commonMax <= commonMin
        commonMax = commonMin + 1;
    end

    ax1 = axes('Parent',figMaps,'Position',[lefts(k) rowTop mapWidth mapHeight]);
    imagesc(ax1,xObs,yObs,cleanMethods{k});
    set(ax1,'YDir','normal','FontName','Arial','FontSize',9);
    axis(ax1,'image'); axis(ax1,'tight'); caxis(ax1,[commonMin commonMax]);
    xlabel(ax1,'x (m)'); ylabel(ax1,'y (m)');
    title(ax1,sprintf('(%c) %s - clean field',char('a'+k-1),methodNames{k}), ...
        'FontName','Arial','FontSize',11,'FontWeight','bold');
    cb1 = colorbar(ax1,'Position',[lefts(k)+mapWidth+0.022 rowTop 0.010 mapHeight]);
    set(cb1,'FontName','Arial','FontSize',9);

    ax2 = axes('Parent',figMaps,'Position',[lefts(k) rowBottom mapWidth mapHeight]);
    imagesc(ax2,xObs,yObs,noisyMethods{k});
    set(ax2,'YDir','normal','FontName','Arial','FontSize',9);
    axis(ax2,'image'); axis(ax2,'tight'); caxis(ax2,[commonMin commonMax]);
    xlabel(ax2,'x (m)'); ylabel(ax2,'y (m)');
    title(ax2,sprintf('(%c) %s - 15%% noise',char('d'+k-1),methodNames{k}), ...
        'FontName','Arial','FontSize',11,'FontWeight','bold');
    cb2 = colorbar(ax2,'Position',[lefts(k)+mapWidth+0.022 rowBottom 0.010 mapHeight]);
    set(cb2,'FontName','Arial','FontSize',9);
end
annotation(figMaps,'textbox',[0.12 0.94 0.76 0.04], ...
    'String','Angular operators: clean field versus 15% Gaussian noise', ...
    'EdgeColor','none','HorizontalAlignment','center', ...
    'FontName','Arial','FontSize',14,'FontWeight','bold');

%% 9. Perfis: legenda fora dos eixos e limites dos corpos bem identificados
[~,iy0] = min(abs(yObs-0));
figProfiles = figure('Color','w','Position',[60 40 1400 950], ...
    'Name','IAF Figure 2 - Multiple sources with noise profiles');
names = {'ISA','Theta','IAF'};
cleanProfiles = {ISA_clean(iy0,:),Theta_clean(iy0,:),IAF_clean(iy0,:)};
noisyProfiles = {ISA_noisy(iy0,:),Theta_noisy(iy0,:),IAF_noisy(iy0,:)};
xPlotLimits = [min(xObs), max(xObs)];
% Draw all four boundaries; use a shared legend below the panels.
for k = 1:3
    ax = subplot(3,1,k); hold(ax,'on');
    hClean = plot(ax,xObs,cleanProfiles{k},'k-','LineWidth',1.8);
    hNoisy = plot(ax,xObs,noisyProfiles{k},'Color',[0.85 0.33 0.10], ...
        'LineWidth',1.0);
    yLimits = ylim(ax);
    hShallow = plot(ax,[shallow.x(1) shallow.x(1)],yLimits, ...
        '--','Color',[0.00 0.45 0.74],'LineWidth',1.2);
    plot(ax,[shallow.x(2) shallow.x(2)],yLimits, ...
        '--','Color',[0.00 0.45 0.74],'LineWidth',1.2,'HandleVisibility','off');
    hDeep = plot(ax,[deep.x(1) deep.x(1)],yLimits, ...
        ':','Color',[0.85 0.33 0.10],'LineWidth',1.5);
    plot(ax,[deep.x(2) deep.x(2)],yLimits, ...
        ':','Color',[0.85 0.33 0.10],'LineWidth',1.5,'HandleVisibility','off');
    ylim(ax,yLimits);
    grid(ax,'on'); box(ax,'on');
    xlim(ax,xPlotLimits);
    set(ax,'FontName','Arial','FontSize',10,'Layer','top');
    ylabel(ax,[names{k} ' (rad)']);
    title(ax,sprintf('(%c) %s - central profile, y = 0 m', ...
        char('a'+k-1),names{k}), ...
        'FontName','Arial','FontSize',11,'FontWeight','bold');
    if k < 3
        set(ax,'XTickLabel',[]);
    else
        xlabel(ax,'x (m)');
    end
    if k == 1
        handlesForLegend = [hClean hNoisy hShallow hDeep];
    end
end
% A single legend below the plots avoids hiding data.
legend(handlesForLegend, ...
    {'Clean field','15% noise','Shallow body limits','Deep body limits'}, ...
    'Orientation','horizontal','Location','southoutside','FontSize',9);

%% 10. Quantificar efeito do ruido nos operadores (diagnostico)
fprintf('\nRMS difference between noisy and clean operators:\n');
fprintf('ISA:   %.6g rad\n', sqrt(mean((ISA_noisy(:)-ISA_clean(:)).^2)));
fprintf('Theta: %.6g rad\n', sqrt(mean((Theta_noisy(:)-Theta_clean(:)).^2)));
fprintf('IAF:   %.6g rad\n', sqrt(mean((IAF_noisy(:)-IAF_clean(:)).^2)));

%% 11. Salvar resultados organizados em arquivos separados
save('Figure2_Multiple_Sources_With_Noise_Data.mat', ...
    'xObs','yObs','Tclean','Tnoisy','TcleanTop','TcleanBottom', ...
    'TnoisyTop','TnoisyBottom','noise','noiseTop','noiseBottom', ...
    'sigmaNoise','noiseLevel','ISA_clean','Theta_clean','IAF_clean', ...
    'ISA_noisy','Theta_noisy','IAF_noisy','ASA_clean','ASA_noisy', ...
    'shallow','deep','inclination','declination');

try
    print(figMaps,'Figure2_Multiple_Sources_With_Noise_Maps.png','-dpng','-r300');
    print(figProfiles,'Figure2_Multiple_Sources_With_Noise_Profiles.png','-dpng','-r300');
    disp('Mapas e perfis da Figura 2 exportados com sucesso.');
catch err
    fprintf('Aviso ao exportar figuras: %s\n',err.message);
end

disp('Figura 2 concluida; parametros e equacoes foram mantidos.');

%% Funcoes locais
function [ISA,Theta,IAF,ASA] = operators_from_derivatives(dTdx,dTdz,epsilon)
    ASA = hypot(dTdx,dTdz);
    ratio = dTdz ./ (ASA + epsilon);
    ratio = max(-1,min(1,ratio));
    ISA = atan(ratio);
    Theta = acos(ratio);
    IAF = (ISA + Theta) ./ 2;
end

function [srcX, srcY, srcZ, cellVol] = prism_cells(xLim, yLim, zLim, n)
    xEdges = linspace(xLim(1), xLim(2), n(1)+1);
    yEdges = linspace(yLim(1), yLim(2), n(2)+1);
    zEdges = linspace(zLim(1), zLim(2), n(3)+1);
    xc = (xEdges(1:end-1) + xEdges(2:end))/2;
    yc = (yEdges(1:end-1) + yEdges(2:end))/2;
    zc = (zEdges(1:end-1) + zEdges(2:end))/2;
    [VX,VY,VZ] = ndgrid(xc,yc,zc);
    srcX = VX(:); srcY = VY(:); srcZ = VZ(:);
    cellVol = (xLim(2)-xLim(1))/n(1) * ...
              (yLim(2)-yLim(1))/n(2) * ...
              (zLim(2)-zLim(1))/n(3);
end

function TnT = total_field_variable_cells(X,Y,zObs,srcX,srcY,srcZ, ...
        cellVol,M,mHat,fHat,mu0_over_4pi)
    Bx = zeros(size(X)); By = zeros(size(X)); Bz = zeros(size(X));
    for k = 1:numel(srcX)
        moment = M*cellVol(k)*mHat;
        rx = X-srcX(k);
        ry = Y-srcY(k);
        rz = zObs-srcZ(k);
        r2 = rx.^2 + ry.^2 + rz.^2;
        r = sqrt(r2);
        r3 = r2.*r;
        r5 = r3.*r2;
        mdotr = moment(1).*rx + moment(2).*ry + moment(3).*rz;
        Bx = Bx + mu0_over_4pi*(3.*rx.*mdotr./r5 - moment(1)./r3);
        By = By + mu0_over_4pi*(3.*ry.*mdotr./r5 - moment(2)./r3);
        Bz = Bz + mu0_over_4pi*(3.*rz.*mdotr./r5 - moment(3)./r3);
    end
    TnT = (fHat(1).*Bx + fHat(2).*By + fHat(3).*Bz)*1e9;
end

function xline_compat(x1, x2)
    yl = ylim;
    plot([x1 x1], yl, 'k--', 'LineWidth', 0.9, 'HandleVisibility', 'off');
    plot([x2 x2], yl, 'k--', 'LineWidth', 0.9, 'HandleVisibility', 'off');
end
