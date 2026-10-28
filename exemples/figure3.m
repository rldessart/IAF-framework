% IAF_Figura3_V5.m
% Figura 3 - Avaliacao inicial da preservacao dos contatos geologicos sob ruido.
% Modelo sintetico independente: dois prismas magnetizados.
% Objetivo: comparar ISA, Theta e IAF em campos limpo e ruidoso.
%
% IMPORTANTE:
% Esta e a V1 de desenvolvimento. Antes de interpretar os resultados como
% evidencia cientifica, confirmar que o modelo e o protocolo de ruido estao
% alinhados com o manuscrito. O criterio de deteccao de contato aqui usado
% e exploratorio e deve ser validado antes da versao final.

clear; close all; clc;

%% 1. Grade de observacao
dxObs = 10;
dyObs = 10;
xObs = -500:dxObs:500;
yObs = -200:dyObs:200;
[Xobs, Yobs] = meshgrid(xObs, yObs);

%% 2. Parametros fisicos do modelo sintetico
M = 1;                         % A/m, valor ilustrativo
inclination = 45;               % graus
declination = 0;                % graus
I = inclination*pi/180;
D = declination*pi/180;
mHat = [cos(I)*cos(D), cos(I)*sin(D), sin(I)];
fHat = mHat;
mu0_over_4pi = 1e-7;

% Corpos sinteticos: limites laterais conhecidos para avaliar os contatos.
shallow.x = [-300, -120];
shallow.y = [-60, 60];
shallow.z = [70, 150];
shallow.n = [9, 6, 6];

deep.x = [220, 500];
deep.y = [-90, 90];
deep.z = [220, 360];
deep.n = [14, 9, 8];

%% 3. Discretizacao dos prismas
[srcX1, srcY1, srcZ1, vol1] = prism_cells(shallow.x, shallow.y, shallow.z, shallow.n);
[srcX2, srcY2, srcZ2, vol2] = prism_cells(deep.x, deep.y, deep.z, deep.n);

srcX = [srcX1; srcX2];
srcY = [srcY1; srcY2];
srcZ = [srcZ1; srcZ2];
cellVol = [repmat(vol1, numel(srcX1), 1); repmat(vol2, numel(srcX2), 1)];

%% 4. Campo total limpo em tres niveis de observacao
zObs = 0;
zTop = -5;
zBottom = 5;

Tclean = total_field_variable_cells(Xobs,Yobs,zObs,srcX,srcY,srcZ,cellVol,M,mHat,fHat,mu0_over_4pi);
TcleanTop = total_field_variable_cells(Xobs,Yobs,zTop,srcX,srcY,srcZ,cellVol,M,mHat,fHat,mu0_over_4pi);
TcleanBottom = total_field_variable_cells(Xobs,Yobs,zBottom,srcX,srcY,srcZ,cellVol,M,mHat,fHat,mu0_over_4pi);

%% 5. Ruido gaussiano reproducivel
noiseLevel = 0.15;
rng(2026,'twister');
sigmaNoise = noiseLevel * std(Tclean(:),0);
noise = sigmaNoise * randn(size(Tclean));
noiseTop = sigmaNoise * randn(size(TcleanTop));
noiseBottom = sigmaNoise * randn(size(TcleanBottom));

Tnoisy = Tclean + noise;
TnoisyTop = TcleanTop + noiseTop;
TnoisyBottom = TcleanBottom + noiseBottom;

%% 6. Derivadas e operadores
[dTdx_clean, ~] = gradient(Tclean,dxObs,dyObs);
dTdz_clean = (TcleanBottom-TcleanTop)/(zBottom-zTop);

[dTdx_noisy, ~] = gradient(Tnoisy,dxObs,dyObs);
dTdz_noisy = (TnoisyBottom-TnoisyTop)/(zBottom-zTop);

epsilon = 1e-6;
[ISA_clean,Theta_clean,IAF_clean,ASA_clean] = operators_from_derivatives(dTdx_clean,dTdz_clean,epsilon);
[ISA_noisy,Theta_noisy,IAF_noisy,ASA_noisy] = operators_from_derivatives(dTdx_noisy,dTdz_noisy,epsilon);

%% 7. Perfis centrais para comparacao visual
[~,iy0] = min(abs(yObs-0));
cleanProfiles = {ISA_clean(iy0,:),Theta_clean(iy0,:),IAF_clean(iy0,:)};
noisyProfiles = {ISA_noisy(iy0,:),Theta_noisy(iy0,:),IAF_noisy(iy0,:)};
methodNames = {'ISA','Theta','IAF'};

figProfiles = figure('Color','w','Position',[80 40 1300 900], ...
    'Name','IAF Figure 3 V5 - Contact preservation');
for k = 1:3
    ax = subplot(3,1,k);
    hold(ax,'on');
    hClean = plot(ax,xObs,cleanProfiles{k},'k-','LineWidth',1.7);
    hNoisy = plot(ax,xObs,noisyProfiles{k},'Color',[0.85 0.33 0.10],'LineWidth',1.0);

    yl = ylim(ax);
    hShallow = plot(ax,[shallow.x(1) shallow.x(1)],yl,'--','Color',[0.00 0.45 0.74],'LineWidth',1.1);
    plot(ax,[shallow.x(2) shallow.x(2)],yl,'--','Color',[0.00 0.45 0.74],'LineWidth',1.1,'HandleVisibility','off');
    hDeep = plot(ax,[deep.x(1) deep.x(1)],yl,':','Color',[0.49 0.18 0.56],'LineWidth',1.3);
    plot(ax,[deep.x(2) deep.x(2)],yl,':','Color',[0.49 0.18 0.56],'LineWidth',1.3,'HandleVisibility','off');
    ylim(ax,yl);
    xlim(ax,[min(xObs),max(xObs)]);
    grid(ax,'on'); box(ax,'on');
    ylabel(ax,[methodNames{k} ' (rad)']);
    title(ax,sprintf('(%c) %s - perfil central, y = 0 m','a'+k-1,methodNames{k}));
    if k == 1
        handlesForLegend = [hClean hNoisy hShallow hDeep];
    end
    if k == 3
        xlabel(ax,'x (m)');
    else
        set(ax,'XTickLabel',[]);
    end
end

% Legenda horizontal compartilhada, abaixo dos tres paineis.
legend(handlesForLegend, ...
    {'Campo limpo','Campo com 15% de ruido','Limites do alvo raso','Limites do alvo profundo'}, ...
    'Orientation','horizontal','Location','southoutside','FontSize',12);

% Solicita ao sistema grafico que atualize as janelas antes de salvar/terminar.
drawnow;

%% 8. Metricas exploratorias de perturbacao e rugosidade
% Estas metricas quantificam alteracao do perfil, mas nao equivalem sozinhas
% a uma deteccao automatica e validada de contato geologico.
fprintf('\nFIGURA 3 V1 - diagnosticos exploratorios\n');
fprintf('Ruido aplicado: %.1f%%; sigma = %.6g nT\n',100*noiseLevel,sigmaNoise);
fprintf('Numero de celulas magnetizadas: %d\n',numel(srcX));

for k = 1:3
    cp = cleanProfiles{k};
    np = noisyProfiles{k};
    rmsPert = sqrt(mean((np-cp).^2));
    roughClean = sqrt(mean(diff(cp,2).^2));
    roughNoisy = sqrt(mean(diff(np,2).^2));
    fprintf('%s: RMS perturbacao = %.6g rad; rugosidade limpa = %.6g; rugosidade ruidosa = %.6g\n', ...
        methodNames{k},rmsPert,roughClean,roughNoisy);
end

%% 9. Avaliacao dos contatos: criterio local com diagnostico de borda
% Para cada contato, procura-se o maximo da magnitude da derivada do perfil
% dentro de uma janela local. O criterio e exploratorio: operadores angulares
% podem ter extremos deslocados em relacao ao contato fisico.
%
% O contato em x=500 m coincide com a borda da grade. Por isso, nao e justo
% usa-lo como medida equivalente aos contatos internos: reportamos esse caso,
% mas calculamos tambem as medias somente com contatos internos.
trueContacts = [shallow.x(1), shallow.x(2), deep.x(1), deep.x(2)];
contactLabels = {'Raso esquerdo','Raso direito','Profundo esquerdo','Profundo direito'};
windowHalfWidth = 100; % metros em torno do contato verdadeiro
profileDx = mean(diff(xObs));
isInterior = (trueContacts > min(xObs)) & (trueContacts < max(xObs));

fprintf('\nAVALIACAO DOS CONTATOS - V5\n');
fprintf('Criterio exploratorio: maximo de |derivada do perfil| em janela de +/- %.0f m\n',windowHalfWidth);
fprintf('Contato                 Verdadeiro   ISA limpo   ISA ruido   Theta limpo Theta ruido IAF limpo  IAF ruido\n');

edgeResults = nan(numel(trueContacts),7);
for j = 1:numel(trueContacts)
    xTrue = trueContacts(j);
    idxWin = find(abs(xObs-xTrue) <= windowHalfWidth);
    edgeResults(j,1) = xTrue;

    for k = 1:3
        cp = cleanProfiles{k};
        np = noisyProfiles{k};
        dcp = abs(gradient(cp,profileDx));
        dnp = abs(gradient(np,profileDx));
        [~,ic] = max(dcp(idxWin));
        [~,in] = max(dnp(idxWin));
        edgeResults(j,2+(k-1)*2) = xObs(idxWin(ic));
        edgeResults(j,3+(k-1)*2) = xObs(idxWin(in));
    end

    fprintf('%-20s %8.0f %11.0f %11.0f %11.0f %11.0f %10.0f %10.0f\n', ...
        contactLabels{j},xTrue,edgeResults(j,2),edgeResults(j,3), ...
        edgeResults(j,4),edgeResults(j,5),edgeResults(j,6),edgeResults(j,7));
end

% Erro absoluto medio calculado em dois modos:
% (1) todos os contatos, para rastreabilidade;
% (2) apenas contatos internos, excluindo x=500 m, que coincide com a borda.
fprintf('\nERRO ABSOLUTO MEDIO DE LOCALIZACAO (m)\n');
fprintf('A media interna exclui o contato em x = %.0f m, na borda da grade.\n',max(xObs));
for k = 1:3
    cleanCol = 2+(k-1)*2;
    noisyCol = cleanCol+1;
    errClean = abs(edgeResults(:,cleanCol)-trueContacts(:));
    errNoisy = abs(edgeResults(:,noisyCol)-trueContacts(:));
    fprintf('%s - todos os contatos: limpo = %.2f m; ruido = %.2f m\n', ...
        methodNames{k},mean(errClean),mean(errNoisy));
    fprintf('%s - contatos internos: limpo = %.2f m; ruido = %.2f m\n', ...
        methodNames{k},mean(errClean(isInterior(:))),mean(errNoisy(isInterior(:))));
end

%% 10. Grafico comparativo do erro medio
% Barras agrupadas para contatos internos, evitando que o contato na borda
% domine a comparacao. O criterio continua exploratorio e deve ser interpretado
% junto da tabela de posicoes estimadas.
meanErrorClean = zeros(1,3);
meanErrorNoisy = zeros(1,3);
for k = 1:3
    cleanCol = 2+(k-1)*2;
    noisyCol = cleanCol+1;
    meanErrorClean(k) = mean(abs(edgeResults(isInterior,cleanCol)-trueContacts(isInterior)'));
    meanErrorNoisy(k) = mean(abs(edgeResults(isInterior,noisyCol)-trueContacts(isInterior)'));
end

figError = figure('Color','w','Position',[120 120 1000 600], ...
    'Name','IAF Figure 3 V5 - Mean contact localization error');
barData = [meanErrorClean(:), meanErrorNoisy(:)];
bar(barData);
set(gca,'XTick',1:3,'XTickLabel',methodNames,'FontName','Arial','FontSize',12);
ylabel('Erro absoluto medio (m)');
title('Erro medio de localizacao dos contatos internos');
legend({'Campo limpo','Campo com 15% de ruido'}, ...
    'Location','northoutside','Orientation','horizontal','FontSize',12);
grid on; box on;
drawnow;

%% 10. Exportar resultados
save('IAF_Figura3_V5_Data.mat', ...
    'xObs','yObs','Tclean','Tnoisy','TcleanTop','TcleanBottom', ...
    'TnoisyTop','TnoisyBottom','noise','noiseTop','noiseBottom', ...
    'sigmaNoise','noiseLevel','ISA_clean','Theta_clean','IAF_clean', ...
    'ISA_noisy','Theta_noisy','IAF_noisy','ASA_clean','ASA_noisy', ...
    'shallow','deep','inclination','declination','trueContacts','edgeResults','windowHalfWidth','meanErrorClean','meanErrorNoisy','isInterior');

try
    print(figProfiles,'IAF_Figura3_V5_Perfis.png','-dpng','-r300');
    print(figError,'IAF_Figura3_V5_Erro_Medio.png','-dpng','-r300');
    disp('Perfis exportados: IAF_Figura3_V5_Perfis.png');
    disp('Grafico de erro exportado: IAF_Figura3_V5_Erro_Medio.png');
catch err
    fprintf('Aviso: nao foi possivel exportar o perfil: %s\n',err.message);
end

drawnow;
disp('IAF_Figura3_V5 concluida.');

%% Funcoes locais
function [ISA,Theta,IAF,ASA] = operators_from_derivatives(dTdx,dTdz,epsilon)
    ASA = hypot(dTdx,dTdz);
    ratio = dTdz ./ (ASA + epsilon);
    ratio = max(-1,min(1,ratio));
    ISA = atan(ratio);
    Theta = acos(ratio);
    IAF = (ISA + Theta)./2;
end

function [srcX,srcY,srcZ,cellVol] = prism_cells(xLim,yLim,zLim,n)
    xEdges = linspace(xLim(1),xLim(2),n(1)+1);
    yEdges = linspace(yLim(1),yLim(2),n(2)+1);
    zEdges = linspace(zLim(1),zLim(2),n(3)+1);
    xc = (xEdges(1:end-1)+xEdges(2:end))/2;
    yc = (yEdges(1:end-1)+yEdges(2:end))/2;
    zc = (zEdges(1:end-1)+zEdges(2:end))/2;
    [VX,VY,VZ] = ndgrid(xc,yc,zc);
    srcX = VX(:); srcY = VY(:); srcZ = VZ(:);
    cellVol = ((xLim(2)-xLim(1))/n(1))* ...
              ((yLim(2)-yLim(1))/n(2))* ...
              ((zLim(2)-zLim(1))/n(3));
end

function TnT = total_field_variable_cells(X,Y,zObs,srcX,srcY,srcZ,cellVol,M,mHat,fHat,mu0_over_4pi)
    Bx = zeros(size(X)); By = zeros(size(X)); Bz = zeros(size(X));
    for k = 1:numel(srcX)
        moment = M*cellVol(k)*mHat;
        rx = X-srcX(k);
        ry = Y-srcY(k);
        rz = zObs-srcZ(k);
        r2 = rx.^2+ry.^2+rz.^2;
        r = sqrt(r2);
        r3 = r2.*r;
        r5 = r3.*r2;
        mdotr = moment(1).*rx+moment(2).*ry+moment(3).*rz;
        Bx = Bx+mu0_over_4pi*(3.*rx.*mdotr./r5-moment(1)./r3);
        By = By+mu0_over_4pi*(3.*ry.*mdotr./r5-moment(2)./r3);
        Bz = Bz+mu0_over_4pi*(3.*rz.*mdotr./r5-moment(3)./r3);
    end
    TnT = (fHat(1).*Bx+fHat(2).*By+fHat(3).*Bz)*1e9;
end
