%% Quarter-ring RKPM beam: learning scaffold
% Build this example one section at a time.  The model is deliberately
% planar: every control point has [u_x, u_y], so the global system has 2*nNode
% unknowns.  Do not use the 3D helix code as a direct template here.
%
% Target benchmark
%   - reference centerline: X(theta) = R*[cos(theta), sin(theta)]
%   - domain: theta in [0, pi/2]
%   - load: downward F/2 at theta = pi/2
%   - symmetry conditions:
%       theta = 0:    u_y = 0,  du_x/deta = 0
%       theta = pi/2: u_x = 0,  du_y/deta = 0
%
% Work order: finish and test one numbered section before moving to the next.
%
% -------------------------------------------------------------------------
% Current notation and array layout
% -------------------------------------------------------------------------
% nNode              number of control points on the reference centerline
% thetaNode          nNode-by-1 circular angle coordinates [rad]
% X                  nNode-by-3 reference positions; X(A,:)=[x_A,y_A,z_A]
% h                  mean chord length between adjacent reference nodes
% alphaPCA/RKPM      dimensionless support scales in rho = alpha*h
% supportPCA         1-by-nNode PCA support radii [length]
% pcaNeighborIds     nNode-by-7 global node indices in each PCA domain
% pcaNeighborCount   nNode-by-1 actual number of valid entries in each row
% tPCA               nNode-by-3 oriented PCA tangent vectors
% etaPCA             nNode-by-7 local projected coordinates [length]
% N, Neta, Netaeta   shape values and first/second eta derivatives; their
%                    columns currently correspond to local PCA-domain slots
% Xeta/Xetaeta       reconstructed X_,eta and X_,etaeta at each node
% metric             ||X_,eta|| = ds/deta at each node
%
% Important: mechanics remains planar (u_x,u_y), although X is stored with
% three coordinates so that the same geometry/PCA notation extends to 3D.

clear; clc; close all;

%% 1. Model parameters and quarter-circle control points
% TODO 1.1: Choose nNode, radius R, E, and rectangular section dimensions.
% TODO 1.2: Set quadratic RKPM order and the dynamic support rule
%            rho = supportScale*h, where h is the mean nodal chord length.
% TODO 1.3: Create thetaNode, then X = [R*cos(thetaNode), R*sin(thetaNode)].
% TODO 1.4: Build etaNode as cumulative chord length (or choose theta as eta,
%            but then consistently include the parameter-to-arc-length factor).

nNode = 129;             % Replace with a value such as 17 after TODO 1.1.
R = 1;
E = 1e6;
totalForce = 100;         % Point-load magnitude F on the complete circular ring.

sectionArea = 0.2;

slenderness = 100;    % R / t
tEffective = R / slenderness;
sectionHeight = 2*sqrt(3)*tEffective;
sectionWidth = sectionArea / sectionHeight;

sectionInertia = sectionWidth * sectionHeight^3 / 12;

nodalArcWeight = zeros(nNode,1);
rkpmOrder = 2;
alphaPCA = 3.5;
alphaRKPM = 2.5;  % Dimensionless alpha in rho = alpha*h.

% nonlinear iteration set up
nLoadStep = 300;
maximumNewtonIteration = 30;
newtonTolerance = 1.0e-9;

%% reference part
% reference geo
thetaNode = linspace(0,pi/2,nNode)';
X = zeros(nNode,3);
X(:,1) = R * cos(thetaNode);
X(:,2) = R * sin(thetaNode);

dX = X(2:nNode,:) - X(1:nNode - 1,:);
segmentLength = vecnorm(dX,2,2);
h = mean(segmentLength);

nDof = 2 * nNode;

supportPCA = alphaPCA * h * ones(1,nNode);

pcaNeighborIds = zeros(nNode,7);
pcaNeighborCount = zeros(nNode,1);
for i = 1:nNode
    pcaCount = 0;
    for j = 1:nNode
        distance = norm(X(i,:) - X(j,:));
        if distance <= supportPCA(j)
            pcaCount = pcaCount + 1;
            pcaNeighborIds(i,pcaCount) = j;
        end
    end
    pcaNeighborCount(i) = pcaCount;
end

tPCA = zeros(nNode,3);
for i = 1:nNode
    pcaIds = pcaNeighborIds(i,1:pcaNeighborCount(i));
    XbarLocal = mean(X(pcaIds,:),1);
    YLocal = X(pcaIds,:) - XbarLocal;
    CPCA = YLocal' * YLocal;
    [pcaVectors,pcaValues] = eig(CPCA);
    [lambdaMax,indexMax] = max(diag(pcaValues));
    tPCA(i,:) = pcaVectors(:,indexMax)';
end
% determine direction
for i = 2:nNode
    if tPCA(i-1,:) * tPCA(i, :)' < 0
        tPCA(i,:) = -tPCA(i,:);
    end
end
startEdge = X(2,:) - X(1,:);
if startEdge * tPCA(1,:)' < 0
    tPCA = -tPCA;
end

etaPCA = zeros(size(pcaNeighborIds));
for i = 1:nNode
    for j = 1:pcaNeighborCount(i)
        etaPCA(i,j) = (X(pcaNeighborIds(i,j),:) - X(i,:)) * tPCA(i,:)';
    end
end

% rkpm parameters

maxNeighbor = ceil(2*alphaRKPM) + 1;
N = zeros(nNode,maxNeighbor);
Neta = zeros(nNode,maxNeighbor);
Netaeta = zeros(nNode,maxNeighbor);
rkpmActiveLocalIds = zeros(nNode,maxNeighbor);
rkpmNeighborIds = zeros(nNode,maxNeighbor);
rkpmNeighborCount = zeros(nNode,1);
supportRKPMGlobal = alphaRKPM * h * ones(1,nNode);
etaRKPM = zeros(size(rkpmNeighborIds));

% Reconstructed reference geometry derivatives and their metric quantities.
Xeta = zeros(nNode,3);
Xetaeta = zeros(nNode,3);
metric = zeros(nNode,1);
arcSegmentLength = (thetaNode(2:nNode) - thetaNode(1:nNode-1)) * R;
nodalArcWeight(2:nNode-1) = arcSegmentLength(1:nNode-2)/2 ...
    + arcSegmentLength(2:nNode-1)/2;
nodalArcWeight(1) = arcSegmentLength(1)/2;
nodalArcWeight(nNode) = arcSegmentLength(end)/2;

for i = 1:nNode
    rkpmCount = 0;
    for j = 1:nNode
        distance = norm(X(i,:) - X(j,:));
        if distance <= supportRKPMGlobal(j)
            rkpmCount = rkpmCount + 1;
            rkpmNeighborIds(i,rkpmCount) = j;
        end
    end
    rkpmNeighborCount(i) = rkpmCount;
end


for i = 1:nNode
    for j = 1:rkpmNeighborCount(i)
        etaRKPM(i,j) = (X(rkpmNeighborIds(i,j),:) - X(i,:)) * tPCA(i,:)';
    end
end


for i = 1:nNode
    rkpmNodeLocal = etaRKPM(i,1:rkpmNeighborCount(i));
    supportRKPMLocal = alphaRKPM * h * ...
        ones(1,rkpmNeighborCount(i));
    [active,NLocal,NetaLocal,NetaetaLocal] = ...
        rkpm_shape_1d(0,rkpmNodeLocal,supportRKPMLocal,rkpmOrder);
    N(i,1:rkpmNeighborCount(i)) = NLocal';
    Neta(i,1:rkpmNeighborCount(i)) = NetaLocal';
    Netaeta(i,1:rkpmNeighborCount(i)) = NetaetaLocal';
    rkpmActiveLocalIds(i,1:rkpmNeighborCount(i)) = active';

    Xeta(i,:) = Neta(i,1:rkpmNeighborCount(i)) ...
        * X(rkpmNeighborIds(i,1:rkpmNeighborCount(i))',:);
    Xetaeta(i,:) = Netaeta(i,1:rkpmNeighborCount(i)) ...
        * X(rkpmNeighborIds(i,1:rkpmNeighborCount(i))',:);
    metric(i) = norm(Xeta(i,:));
end

%% Boundary conditions and external force
% At theta = 0:     u_y = 0,  u_x,eta = 0
% At theta = pi/2:  u_x = 0,  u_y,eta = 0
C = zeros(4,nDof);
g = zeros(4,1);

bcNodalBegin = 1;
NetaNodalBegin = Neta( ...
    bcNodalBegin,1:rkpmNeighborCount(bcNodalBegin));
NNodalBegin = N( ...
    bcNodalBegin,1:rkpmNeighborCount(bcNodalBegin));
beginNodeIds = rkpmNeighborIds( ...
    bcNodalBegin,1:rkpmNeighborCount(bcNodalBegin));

C(1,2*beginNodeIds) = NNodalBegin;
C(2,2*beginNodeIds-1) = NetaNodalBegin;

bcNodalEnd = nNode;
NetaNodalEnd = Neta( ...
    bcNodalEnd,1:rkpmNeighborCount(bcNodalEnd));
NNodalEnd = N( ...
    bcNodalEnd,1:rkpmNeighborCount(bcNodalEnd));
endNodeIds = rkpmNeighborIds( ...
    bcNodalEnd,1:rkpmNeighborCount(bcNodalEnd));

C(3,2*endNodeIds-1) = NNodalEnd;
C(4,2*endNodeIds) = NetaNodalEnd;

Fext = zeros(nDof,1);
quarterRingLoad = totalForce / 2;
Fext(2*nNode) = -quarterRingLoad;



%% nonlinear part force evaluated at the current displacement d


d = zeros(nNode*2,1);
Z = null(C);

for loadStep = 1:nLoadStep

    loadFactor = loadStep / nLoadStep;
    FextCurrent = loadFactor * Fext;

    converged = false;

    for iteration = 1:maximumNewtonIteration

        FmemNonlinear = zeros(nNode*2,1);
        FbendNonlinear = zeros(nNode*2,1);
        KmatMembrane = zeros(nNode*2,nNode*2);
        KgeoMembrane = zeros(nNode*2,nNode*2);
        KmatBending = zeros(nNode*2,nNode*2);
        KgeoBending = zeros(nNode*2,nNode*2);
        for i = 1:nNode

            % local nodes and local DOFs
            localNodeIds = rkpmNeighborIds( ...
                i,1:rkpmNeighborCount(i))';

            nLocalNode = numel(localNodeIds);

            localDofIds = [2*localNodeIds-1, 2*localNodeIds];
            localDofIds = reshape(localDofIds',[],1);
            % Current local displacement vector
            dLocal = d(localDofIds);

            NetaLocal = Neta( ...
                i,1:rkpmNeighborCount(i))';
            NetaetaLocal = Netaeta( ...
                i,1:rkpmNeighborCount(i))';

            % current state

            state = evaluateCurrentConfiguration2D( ...
                dLocal, Xeta(i,1:2)', Xetaeta(i,1:2)', metric(i), ...
                NetaLocal, NetaetaLocal);

            % Nodal integration weight in reference arc length
            weight = nodalArcWeight(i);

            [localMembraneForce, ...
                localKmatMembrane, ...
                localKgeoMembrane] = ...
                evaluateQuarterRingMembrane2D( ...
                    E, sectionArea, weight, state);

            FmemNonlinear(localDofIds) = ...
                FmemNonlinear(localDofIds) ...
                + localMembraneForce;

            KmatMembrane(localDofIds,localDofIds) = ...
                KmatMembrane(localDofIds,localDofIds) ...
                + localKmatMembrane;

            KgeoMembrane(localDofIds,localDofIds) = ...
                KgeoMembrane(localDofIds,localDofIds) ...
                + localKgeoMembrane;

            [localFb, localKmatB, localKgeoB] = ...
                evaluateQuarterRingBending2D( ...
                    dLocal, Xeta(i,1:2)', Xetaeta(i,1:2)', metric(i), ...
                    Neta(i,1:rkpmNeighborCount(i))', ...
                    Netaeta(i,1:rkpmNeighborCount(i))', ...
                    E, sectionInertia, nodalArcWeight(i));

            FbendNonlinear(localDofIds) = ...
                FbendNonlinear(localDofIds) + localFb;
            KmatBending(localDofIds,localDofIds) = ...
                KmatBending(localDofIds,localDofIds) + localKmatB;
            KgeoBending(localDofIds,localDofIds) = ...
                KgeoBending(localDofIds,localDofIds) + localKgeoB;
        end

        Fint = FmemNonlinear + FbendNonlinear;
        Kt = KmatMembrane + KgeoMembrane ...
            + KmatBending + KgeoBending;

        % Current residual
        residual = Fint - FextCurrent;

        % Check convergence
        projectedResidual = Z' * residual;
        relativeResidual = norm(projectedResidual) ...
            / max(1.0,norm(Z' * FextCurrent));
        constraintResidual = norm(C*d-g);

        fprintf('step %d, iteration %d, residual = %.3e\n', ...
            loadStep, iteration, relativeResidual);

        if relativeResidual < newtonTolerance ...
                && constraintResidual < newtonTolerance
            converged = true;
            break;
        end

        % Newton correction
        Kaug = [Kt, C';
                C,  zeros(size(C,1),size(C,1))];

        rhs = [-residual;
                g-C*d];

        increment = Kaug \ rhs;
        deltaD = increment(1:nNode*2);

        % Update displacement
        d = d + deltaD;

    end


    if ~converged
        error('Current load step did not converge.');
    end

end



KtMembrane = KmatMembrane + KgeoMembrane;

KtBending = KmatBending + KgeoBending;

FintTrial = FmemNonlinear + FbendNonlinear;
KtTrial = KtMembrane + KtBending;

fprintf('nonlinear membrane force norm = %.6e\n', ...
    norm(FmemNonlinear));

fprintf('membrane material tangent norm = %.6e\n', ...
    norm(KmatMembrane,'fro'));

fprintf('membrane geometric tangent norm = %.6e\n', ...
    norm(KgeoMembrane,'fro'));

fprintf('nonlinear bending force norm = %.6e\n', ...
    norm(FbendNonlinear));

fprintf('bending material tangent norm = %.6e\n', ...
    norm(KmatBending,'fro'));

fprintf('bending geometric tangent norm = %.6e\n', ...
    norm(KgeoBending,'fro'));

fprintf('total trial tangent symmetry error = %.6e\n', ...
    norm(KtTrial-KtTrial','fro') ...
    / max(norm(KtTrial,'fro'),eps));

% internel force written explicity
% Fint = Klinear * d;
%
% residual = Fint + C' * lambdaConstraint - Fext;

dx = d(1:2:end);
dy = d(2:2:end);
dz = zeros(size(dx));

fprintf('free-DOF equilibrium residual = %.4e\n' , ...
    norm(Z' * residual));
fprintf('constraint residual  = %.4e\n', norm(C*d-g));


DX = zeros(nNode,1);
DY = DX;
DZ = DX;

for i = 1:nNode
    ids_local = i;
    neighborHood = 1:rkpmNeighborCount(ids_local);
    NNodal_local = N(ids_local,neighborHood);
    globalIds = rkpmNeighborIds(i,neighborHood);
    DX(i) = NNodal_local * dx(globalIds);
    DY(i) = NNodal_local * dy(globalIds);
    DZ(i) = NNodal_local * dz(globalIds);
end


% plot figure
deformationScale = 1;

U = [DX,DY,DZ];
Xdeformed = X + deformationScale * U;

figure('Color','w');
hold on;
grid on;
box on;

plot(X(:,1), X(:,2), 'k--', ...
    'LineWidth', 1.5, ...
    'DisplayName', 'Reference');

plot(Xdeformed(:,1), Xdeformed(:,2), 'r-', ...
    'LineWidth', 2.0, ...
    'DisplayName', 'Deformed');
scatter(X(1,1), X(1,2), 45, 'k', 'filled', ...
    'DisplayName', 'Symmetry point');
scatter(X(end,1), X(end,2), 45, 'b', 'filled', ...
    'DisplayName', 'Load point A');

axis equal;
xlabel('x');
ylabel('y');
title(sprintf('Quarter ring deformation (scale = %g)', deformationScale));
legend('Location','best');


% Analytical solution check at load point A
% The paper's t is an equivalent thickness/radius-of-gyration measure.
% For the current rectangular section, tEffective = sqrt(I/A).

% The analytical benchmark reports the transverse displacement at side
% point A, theta = 0. The load is applied at the top point theta = pi/2,
% so compare the horizontal side-point displacement with the benchmark
% sign convention.
uA_RKPM = -DX(1);

uA_exact = -totalForce * R^3 / (E * sectionInertia) ...
    * ((4 - pi) / (4*pi) - 0.25*(tEffective/R)^2);

relativeError_uA = abs(abs(uA_RKPM) - abs(uA_exact)) / abs(uA_exact);

fprintf('u_A (RKPM, side point theta=0)  = %.8e\n', uA_RKPM);
fprintf('u_A (exact, side point theta=0) = %.8e\n', uA_exact);
fprintf('relative error = %.4f %%\n', 100*relativeError_uA);
