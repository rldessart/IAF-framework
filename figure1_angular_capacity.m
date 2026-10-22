% IAF Figure 1 - Angular capacity under magnetization inclinations
% Self-contained GNU Octave / MATLAB script.
%
% PURPOSE
%   Build a synthetic magnetic anomaly from a fixed 3-D rectangular prism
%   represented by volume-cell dipoles; vary the magnetization inclination
%   (30, 45, 60 degrees); compute ISA, Theta, and IAF; export maps, profiles,
%   and numerical data.
%
% MODEL LIMITATIONS
%   The field is a numerical volume-dipole approximation, not a closed-form
%   analytical prism solution. The magnetization magnitude is fixed and the
%   magnetization/inducing-field directions are assumed coincident.
%   The operator equations below are the current project definitions and
%   must be checked against the final methodological derivation before
%   publication. No Gaussian smoothing or noise is applied.
%
% COORDINATES
%   x/y horizontal; z positive downward; observation surface z = 0.
%
% Run from this file's directory, or set data/output paths below.

clear; close all; clc;

%% -------------------- Configuration --------------------
outDir = 'IAF_Figure1_output';
if exist(outDir,'dir') ~= 7
    mkdir(outDir);
end

dxObs = 10; dyObs = 10;
xObs = 0:dxObs:1000;
yObs = -200:dyObs:200;
[Xobs,Yobs] = meshgrid(xObs,yObs);

% Fixed prism geometry (m)
xLim = [400 600];
yLim = [-75 75];
zLim = [100 250];

% Volume-cell discretization (production grid validated against a finer
% numerical quadrature reference in the project's separate convergence test).
nx = 12; ny = 9; nz = 9;
xCell = linspace(xLim(1),xLim(2),nx+1);
yCell = linspace(yLim(1),yLim(2),ny+1);
zCell = linspace(zLim(1),zLim(2),nz+1);
xc = (xCell(1:end-1)+xCell(2:end))/2;
yc = (yCell(1:end-1)+yCell(2:end))/2;
zc = (zCell(1:end-1)+zCell(2:end))/2;
[VX,VY,VZ] = ndgrid(xc,yc,zc);
cellVol = ((xLim(2)-xLim(1))/nx) * ...
          ((yLim(2)-yLim(1))/ny) * ...
          ((zLim(2)-zLim(1))/nz);
srcX = VX(:); srcY = VY(:); srcZ = VZ(:);

M = 1.0;                        % A/m, illustrative fixed magnitude
inclinations = [30 45 60];      % degrees
mu0_over_4pi = 1e-7;
epsGrad = 1e-12;
zAbove = -5;                    % m; above z=0 surface
zBelow =  5;                    % m; positive downward
nAngles = numel(inclinations);
results = struct();

fprintf('IAF Figure 1: fixed prism, volume-dipole forward model\n');
fprintf('Grid: %d x %d observation points; prism cells: %d\n', ...
    numel(xObs),numel(yObs),numel(srcX));

%% -------------------- Forward model and operators --------------------
for ia = 1:nAngles
    inc = inclinations(ia);
    I = inc*pi/180;

    % Inclination measured from horizontal, positive downward.
    mHat = [cos(I), 0, sin(I)];
    fHat = mHat;  % induced magnetization assumption

    T0 = totalFieldDipoles(Xobs,Yobs,0,srcX,srcY,srcZ, ...
        cellVol,M,mHat,fHat,mu0_over_4pi);
    Tbelow = totalFieldDipoles(Xobs,Yobs,zBelow,srcX,srcY,srcZ, ...
        cellVol,M,mHat,fHat,mu0_over_4pi);
    Tabove = totalFieldDipoles(Xobs,Yobs,zAbove,srcX,srcY,srcZ, ...
        cellVol,M,mHat,fHat,mu0_over_4pi);

    % Horizontal derivative along x. Vertical derivative uses a centered
    % finite difference across the observation surface.
    [dTdx,~] = gradient(T0,dxObs,dyObs);
    dTdz = (Tbelow-Tabove)/(zBelow-zAbove);

    ASA = hypot(dTdx,dTdz);
    ratio = dTdz ./ (ASA + epsGrad);
    ratio = max(-1,min(1,ratio));

    % Current manuscript/project definitions. Do not silently reinterpret.
    ISA = atan(ratio);
    Theta = acos(ratio);
    IAF = (ISA + Theta)/2;

    iy0 = find(yObs==0,1);
    results(ia).inclination = inc;
    results(ia).T = T0;
    results(ia).dTdx = dTdx;
    results(ia).dTdz = dTdz;
    results(ia).ASA = ASA;
    results(ia).ISA = ISA;
    results(ia).Theta = Theta;
    results(ia).IAF = IAF;
    results(ia).profileT = T0(iy0,:);
    results(ia).profileISA = ISA(iy0,:);
    results(ia).profileTheta = Theta(iy0,:);
    results(ia).profileIAF = IAF(iy0,:);

    fprintf('Inclination %2d deg: max|T| = %.6f nT\n', ...
        inc,max(abs(T0(:))));
end

%% -------------------- Figure 1A: maps --------------------
figMaps = figure('Color','w','Position',[40 40 1550 1000], ...
    'Name','Figure 1 - Angular response maps','NumberTitle','off');
operators = {'ISA','Theta','IAF'};
for ia = 1:nAngles
    for io = 1:3
        ax = subplot(3,3,(ia-1)*3+io);
        A = results(ia).(operators{io});
        sc = max(abs(A(:)));
        if sc > 0, Aplot = A/sc; else, Aplot = A; end
        imagesc(xObs,yObs,Aplot);
        set(ax,'YDir','normal','FontName','Arial','FontSize',11, ...
            'LineWidth',0.8,'Layer','top');
        axis tight;
        caxis([-1 1]);
        xlabel('x (m)','FontSize',11);
        ylabel('y (m)','FontSize',11);
        letter = char('a'+(ia-1)*3+io-1);
        title(sprintf('(%c) %s | I = %d deg',letter, ...
            operators{io},inclinations(ia)), ...
            'FontSize',12,'FontWeight','bold');
        cb = colorbar;
        set(cb,'FontName','Arial','FontSize',10);
    end
end
colormap(jet);
print(figMaps,fullfile(outDir,'IAF_Figure1_Maps.png'),'-dpng','-r300');

%% -------------------- Figure 1B: central profiles --------------------
% Manual normalized-figure annotations are used for the legend to prevent
% Octave versions from moving it over a title or data curve.
figProfiles = figure('Color','w','Position',[80 40 1400 1000], ...
    'Name','Figure 1 - Central profiles','NumberTitle','off');

% Legend items in a dedicated header band (line + text annotations).
legendY = 0.955;
legendX = [0.39 0.49 0.61];
legendColors = {[0.000 0.447 0.741],[0.850 0.325 0.098],[0.929 0.694 0.125]};
legendLabels = {'ISA','Theta','IAF'};
for k = 1:3
    annotation(figProfiles,'line',[legendX(k) legendX(k)+0.035], ...
        [legendY legendY],'Color',legendColors{k},'LineWidth',2.5);
    annotation(figProfiles,'textbox',[legendX(k)+0.04 legendY-0.018 0.06 0.035], ...
        'String',legendLabels{k},'EdgeColor','none','FontName','Arial', ...
        'FontSize',12,'VerticalAlignment','middle');
end

% Footer annotation avoids collision with legend and profile titles.
annotation(figProfiles,'textbox',[0.18 0.018 0.72 0.03], ...
    'String','Dashed vertical lines indicate the source footprint (x = 400-600 m).', ...
    'EdgeColor','none','HorizontalAlignment','center', ...
    'FontName','Arial','FontSize',10);

panelBottom = [0.665 0.385 0.105]; % top-to-bottom: 30, 45, 60 deg
panelLetters = {'(j)','(k)','(l)'};
for ia = 1:nAngles
    ax = axes('Parent',figProfiles,'Position',[0.17 panelBottom(ia) 0.79 0.19]);
    hold(ax,'on');

    pISA = results(ia).profileISA(:)';
    pTheta = results(ia).profileTheta(:)';
    pIAF = results(ia).profileIAF(:)';
    pISA = pISA/(max(abs(pISA))+eps);
    pTheta = pTheta/(max(abs(pTheta))+eps);
    pIAF = pIAF/(max(abs(pIAF))+eps);

    plot(ax,xObs,pISA,'LineWidth',2.0,'Color',legendColors{1});
    plot(ax,xObs,pTheta,'LineWidth',2.0,'Color',legendColors{2});
    plot(ax,xObs,pIAF,'LineWidth',2.2,'Color',legendColors{3});
    yl = [-1.08 1.08];
    for k = 1:2
        plot(ax,[xLim(k) xLim(k)],yl,'k--','LineWidth',1.0);
    end
    set(ax,'FontName','Arial','FontSize',11,'LineWidth',0.8, ...
        'Layer','top','Box','on');
    grid(ax,'on');
    xlim(ax,[0 1000]); ylim(ax,yl);
    set(ax,'YTick',[-1 -0.5 0 0.5 1]);
    ylabel(ax,'Normalized response','FontName','Arial','FontSize',11);
    title(ax,sprintf('%s Central profile - Magnetization inclination = %d deg', ...
        panelLetters{ia},inclinations(ia)), ...
        'FontName','Arial','FontSize',12,'FontWeight','bold');
    if ia < nAngles
        set(ax,'XTickLabel',[]);
    else
        xlabel(ax,'x (m)','FontName','Arial','FontSize',11);
    end
end
print(figProfiles,fullfile(outDir,'IAF_Figure1_Profiles.png'),'-dpng','-r300');

%% -------------------- Save numerical data and metadata --------------------
save(fullfile(outDir,'IAF_Figure1_Data.mat'), ...
    'xObs','yObs','inclinations','results','xLim','yLim','zLim', ...
    'nx','ny','nz','M','dxObs','dyObs','zAbove','zBelow');

fid = fopen(fullfile(outDir,'README_results.txt'),'w');
fprintf(fid,'IAF Figure 1 outputs\n');
fprintf(fid,'Forward model: volume-cell dipole approximation (not closed-form analytical prism).\n');
fprintf(fid,'Prism x=[%g,%g] m; y=[%g,%g] m; z=[%g,%g] m.\n', ...
    xLim(1),xLim(2),yLim(1),yLim(2),zLim(1),zLim(2));
fprintf(fid,'Discretization: %d x %d x %d cells.\n',nx,ny,nz);
fprintf(fid,'Inclinations: 30, 45, 60 degrees; no Gaussian smoothing or noise.\n');
fprintf(fid,'Operators follow the current project definitions: ISA=atan(dTdz/(ASA+eps)); Theta=acos(dTdz/(ASA+eps)); IAF=(ISA+Theta)/2.\n');
fprintf(fid,'Map panels and profile curves are normalized individually for shape comparison.\n');
fclose(fid);

fprintf('\nOutputs saved in %s:\n',outDir);
fprintf('  IAF_Figure1_Maps.png\n  IAF_Figure1_Profiles.png\n  IAF_Figure1_Data.mat\n  README_results.txt\n');
fprintf('Reminder: validate operator definitions before scientific interpretation/publication.\n');

%% -------------------- Local function --------------------
function TnT = totalFieldDipoles(X,Y,zObs,srcX,srcY,srcZ, ...
        cellVol,M,mHat,fHat,mu0_over_4pi)
    % Point-dipole field:
    % B = mu0/(4*pi) * (3*r*(m dot r)/|r|^5 - m/|r|^3)
    % Total-field anomaly approximation: dT ~= dot(F_hat,B).
    Bx = zeros(size(X)); By = zeros(size(X)); Bz = zeros(size(X));
    moment = M*cellVol*mHat;
    for k = 1:numel(srcX)
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
    T_Tesla = fHat(1).*Bx + fHat(2).*By + fHat(3).*Bz;
    TnT = T_Tesla*1e9;
end
