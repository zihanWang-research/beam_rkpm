%% Quarter-ring RKPM beam: geometrically nonlinear learning program
% This file follows the section-based style of the linear quarter-ring
% programs.  Every control point has [u_x,u_y], so nDof = 2*nNode.
% The current tangent is obtained by central differences of the complete
% internal force; it can later be replaced by the analytic Kmat + Kgeo form.

clear; clc; close all;

%% 1. Model parameters and reference quarter-circle

R = 1.0;
E = 1.0e6;
sectionWidth = 0.2;
slenderness = 100.0;                 % R / effective thickness
tEffective = R / slenderness;
sectionHeight = 2.0 * sqrt(3.0) * tEffective;
sectionArea = sectionWidth * sectionHeight;
sectionInertia = sectionWidth * sectionHeight^3 / 12.0;
axialRigidity = E * sectionArea;
bendingRigidity = E * sectionInertia;

totalForce = 5.0;
quarterLoad = totalForce / 2.0;
nNode = 17;
rkpmOrder = 2;
alphaRKPM = 3.2;                    % support radius = alphaRKPM * h
pointsPerCell = 2;

arcLength = 0.5 * pi * R;
sNode = linspace(0.0, arcLength, nNode).';
thetaNode = sNode / R;
XNode = R * [cos(thetaNode), sin(thetaNode)];
dX = XNode(2:nNode,:) - XNode(1:nNode-1,:);
segmentLength = vecnorm(dX,2,2);
h = mean(segmentLength);
supportRKPM = alphaRKPM * h * ones(nNode, 1);
nDof = 2 * nNode;

nLoadStep = 20;
maximumNewtonIteration = 25;
newtonTolerance = 1.0e-9;
finiteDifferenceScale = 1.0e-7;

%% 2. Gauss-point data and RKPM shape functions

data = buildQuarterRingGaussData(R, sNode, supportRKPM, ...
    rkpmOrder, pointsPerCell);
fprintf('Quarter ring: %d nodes, %d Gauss points, EA = %.6e, EI = %.6e\n', ...
    nNode, numel(data.weight), axialRigidity, bendingRigidity);

%% 3. External load

[loadActive, NLoad, ~, ~] = rkpm_shape_1d(arcLength, sNode, supportRKPM, rkpmOrder);
Fext = zeros(nDof, 1);
Fext(2 * loadActive) = -quarterLoad * NLoad;

%% 4. Symmetry constraints

[C, g] = buildQuarterRingConstraints(arcLength, sNode, supportRKPM, ...
    rkpmOrder, nDof);
if rank(C) < size(C, 1)
    error('Quarter-ring constraint matrix is rank deficient.');
end
Z = null(C);

%% 5. Load stepping and Newton-Raphson iteration

d = zeros(nDof, 1);
history = zeros(nLoadStep * maximumNewtonIteration, 5);
historyRow = 0;

for loadStep = 1:nLoadStep
    loadFactor = loadStep / nLoadStep;
    converged = false;

    for iteration = 1:maximumNewtonIteration
        Fint = assembleInternalForce(d, data, ...
            axialRigidity, bendingRigidity, nDof);
        residual = Fint - loadFactor * Fext;
        projectedResidual = Z.' * residual;
        residualNorm = norm(projectedResidual);
        referenceNorm = max(1.0, norm(loadFactor * Z.' * Fext));
        relativeResidual = residualNorm / referenceNorm;
        constraintResidual = norm(C * d - g);

        historyRow = historyRow + 1;
        history(historyRow, :) = [loadStep, iteration, loadFactor, ...
            relativeResidual, constraintResidual];
        fprintf('step %02d/%02d, iter %02d: projected residual %.3e\n', ...
            loadStep, nLoadStep, iteration, relativeResidual);

        if relativeResidual < newtonTolerance && constraintResidual < newtonTolerance
            converged = true;
            break;
        end

        Kt = finiteDifferenceTangent(d, data, axialRigidity, ...
            bendingRigidity, nDof, finiteDifferenceScale);
        Kaug = [Kt, C.'; C, zeros(size(C, 1), size(C, 1))];
        rhs = [-residual; g - C * d];
        increment = Kaug \ rhs;
        deltaD = increment(1:nDof);
        if any(~isfinite(deltaD))
            error('Newton increment contains NaN or Inf.');
        end
        d = d + deltaD;
    end

    if ~converged
        error('Load step %d did not converge.', loadStep);
    end
end
history = history(1:historyRow, :);

%% 6. Recover final displacement and plot the deformed centerline

FintFinal = assembleInternalForce(d, data, ...
    axialRigidity, bendingRigidity, nDof);
finalResidual = FintFinal - Fext;
finalProjectedResidual = Z.' * finalResidual;
dNode = reshape(d, 2, []).';
[beginActive, NBegin, ~, ~] = rkpm_shape_1d(0.0, sNode, supportRKPM, rkpmOrder);
uA = NBegin.' * dNode(beginActive, 1);
fprintf('Final projected residual: %.3e\n', norm(finalProjectedResidual));
fprintf('Displacement at s = 0 symmetry point: ux = %.6e, uy = %.6e\n', uA, ...
    NBegin.' * dNode(beginActive, 2));

SPlot = linspace(0.0, arcLength, 240).';
XRefPlot = R * [cos(SPlot / R), sin(SPlot / R)];
XDefPlot = zeros(size(XRefPlot));
for k = 1:numel(SPlot)
    [active, N, ~, ~] = rkpm_shape_1d(SPlot(k), sNode, supportRKPM, rkpmOrder);
    XDefPlot(k, :) = XRefPlot(k, :) + N.' * dNode(active, :);
end
maximumDisplacement = max(vecnorm(XDefPlot - XRefPlot, 2, 2));
displayScale = 0.18 * R / max(maximumDisplacement, eps);
XDisplay = XRefPlot + displayScale * (XDefPlot - XRefPlot);

figure('Color', 'w');
plot(XRefPlot(:, 1), XRefPlot(:, 2), 'k--', 'LineWidth', 1.2); hold on;
plot(XDisplay(:, 1), XDisplay(:, 2), 'b-', 'LineWidth', 1.8);
plot(XNode(:, 1), XNode(:, 2), 'ko', 'MarkerFaceColor', 'w');
plot(XNode(:, 1) + displayScale * dNode(:, 1), ...
    XNode(:, 2) + displayScale * dNode(:, 2), 'ro', 'MarkerFaceColor', 'r');
axis equal; grid on; xlabel('x'); ylabel('y');
legend('reference', 'deformed (scaled)', 'reference nodes', 'deformed nodes', ...
    'Location', 'best');
title('Geometrically nonlinear quarter-ring');

figure('Color', 'w');
semilogy(history(:, 3), max(history(:, 4), eps), 'o-');
grid on; xlabel('load factor'); ylabel('projected relative residual');
title('Newton convergence history');

function data = buildQuarterRingGaussData(R, SNode, support, order, nGaussPerCell)
[parentPoint, parentWeight] = gaussRule(nGaussPerCell);
nCell = numel(SNode) - 1;
nPoint = nCell * nGaussPerCell;
data = struct();
data.S = zeros(nPoint, 1);
data.weight = zeros(nPoint, 1);
data.active = cell(nPoint, 1);
data.N = cell(nPoint, 1);
data.N_S = cell(nPoint, 1);
data.N_SS = cell(nPoint, 1);
data.A = zeros(nPoint, 2);
data.C = zeros(nPoint, 2);
data.kappa0 = zeros(nPoint, 1);
q = 0;
for cellIndex = 1:nCell
    SLeft = SNode(cellIndex);
    SRight = SNode(cellIndex + 1);
    cellJacobian = 0.5 * (SRight - SLeft);
    for gaussIndex = 1:nGaussPerCell
        q = q + 1;
        S = 0.5 * (SLeft + SRight) + cellJacobian * parentPoint(gaussIndex);
        [active, N, N_S, N_SS] = rkpm_shape_1d(S, SNode, support, order);
        theta = S / R;
        data.S(q) = S;
        data.weight(q) = cellJacobian * parentWeight(gaussIndex);
        data.active{q} = active;
        data.N{q} = N;
        data.N_S{q} = N_S;
        data.N_SS{q} = N_SS;
        data.A(q, :) = [-sin(theta), cos(theta)];
        data.C(q, :) = [-cos(theta) / R, -sin(theta) / R];
        data.kappa0(q) = 1.0 / R;
    end
end
end

function Fint = assembleInternalForce(d, data, EA, EI, nDof)
Fint = zeros(nDof, 1);
for q = 1:numel(data.weight)
    active = data.active{q};
    N_S = data.N_S{q};
    N_SS = data.N_SS{q};
    localDofs = reshape([2 * active - 1, 2 * active].', [], 1);
    dLocal = d(localDofs);
    nActive = numel(active);
    H1 = zeros(2, 2 * nActive);
    H2 = zeros(2, 2 * nActive);
    for a = 1:nActive
        columns = 2 * a - 1 : 2 * a;
        H1(:, columns) = N_S(a) * eye(2);
        H2(:, columns) = N_SS(a) * eye(2);
    end

    A = data.A(q, :).';
    C = data.C(q, :).';
    aCurrent = A + H1 * dLocal;
    cCurrent = C + H2 * dLocal;
    lambda = norm(aCurrent);
    if lambda < 1.0e-12
        error('Current tangent vector collapsed at Gauss point %d.', q);
    end

    membraneStrain = 0.5 * (aCurrent.' * aCurrent - 1.0);
    determinant = determinant2D(aCurrent, cCurrent);
    curvature = determinant / lambda^3;
    bendingStrain = curvature - data.kappa0(q);
    Bm = aCurrent.' * H1;

    J = [0.0, -1.0; 1.0, 0.0];
    gradientA = -J * cCurrent / lambda^3 ...
        - 3.0 * determinant * aCurrent / lambda^5;
    gradientC = J * aCurrent / lambda^3;
    Bb = zeros(1, 2 * nActive);
    for a = 1:nActive
        columns = 2 * a - 1 : 2 * a;
        Bb(columns) = N_S(a) * gradientA.' + N_SS(a) * gradientC.';
    end

    localForce = (Bm.' * (EA * membraneStrain) ...
        + Bb.' * (EI * bendingStrain)) * data.weight(q);
    Fint(localDofs) = Fint(localDofs) + localForce;
end
end

function Kt = finiteDifferenceTangent(d, data, EA, EI, nDof, relativeStep)
Kt = zeros(nDof, nDof);
for column = 1:nDof
    step = relativeStep * max(1.0, abs(d(column)));
    dPlus = d; dPlus(column) = dPlus(column) + step;
    dMinus = d; dMinus(column) = dMinus(column) - step;
    Kt(:, column) = (assembleInternalForce(dPlus, data, EA, EI, nDof) ...
        - assembleInternalForce(dMinus, data, EA, EI, nDof)) / (2.0 * step);
end
Kt = 0.5 * (Kt + Kt.');
end

function [C, g] = buildQuarterRingConstraints(L, SNode, support, order, nDof)
[activeA, ~, dNA, ~] = rkpm_shape_1d(0.0, SNode, support, order);
[activeB, ~, dNB, ~] = rkpm_shape_1d(L, SNode, support, order);
C = zeros(4, nDof);
for a = 1:numel(activeA)
    node = activeA(a);
    C(1, 2 * node) = C(1, 2 * node) + 1.0;       % uy(0) = 0
    C(2, 2 * node - 1) = C(2, 2 * node - 1) + dNA(a); % ux,S(0) = 0
end
for a = 1:numel(activeB)
    node = activeB(a);
    C(3, 2 * node - 1) = C(3, 2 * node - 1) + 1.0; % ux(L) = 0
    C(4, 2 * node) = C(4, 2 * node) + dNB(a);       % uy,S(L) = 0
end
g = zeros(4, 1);
end

function value = determinant2D(a, b)
value = a(1) * b(2) - a(2) * b(1);
end

function [point, weight] = gaussRule(n)
switch n
    case 1
        point = 0.0; weight = 2.0;
    case 2
        point = [-1.0; 1.0] / sqrt(3.0); weight = [1.0; 1.0];
    case 3
        point = [-sqrt(3.0 / 5.0); 0.0; sqrt(3.0 / 5.0)];
        weight = [5.0; 8.0; 5.0] / 9.0;
    otherwise
        error('Only 1, 2, or 3 Gauss points are implemented.');
end
end
