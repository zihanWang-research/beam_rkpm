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

nNode = 33;             % Replace with a value such as 17 after TODO 1.1.
R = 1;
E = 1e6;
totalForce = 20;         % Point-load magnitude F on the complete circular ring.

sectionArea = 0.2;
density = 1.0;

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
nLoadStep = 200;
maximumNewtonIteration = 30;
newtonTolerance = 1.0e-9;

% dynamic part set up
dt = 0.01;
nTimeStep = 100;
betat = 0.25;
gammat = 0.5;
movieEvery = 5;          % Write one video frame every this many time steps.
movieFrameRate = 15;     % Video frame rate [frames/s].
deformationScale = 1;

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

% mass matrix
M = zeros(nDof);
for i = 1:nNode
    localNodeIds = rkpmNeighborIds( ...
        i,1:rkpmNeighborCount(i))';

    nLocalNode = numel(localNodeIds);

    localDofIds = [2*localNodeIds-1, 2*localNodeIds];
    localDofIds = reshape(localDofIds',[],1);
    H0 = zeros(2,2*nLocalNode);
    NLocal = N(i,1:rkpmNeighborCount(i))';
    for localNode = 1:nLocalNode
        columns = 2*localNode-1 : 2*localNode;
        H0(:,columns) = NLocal(localNode) * eye(2);
    end
    Mlocal = density * sectionArea * H0' * H0 * nodalArcWeight(i);
    M(localDofIds,localDofIds) = M(localDofIds,localDofIds) + Mlocal;
end

%% nonlinear part force evaluated at the current displacement d


d0 = zeros(nNode*2,1);
v0 = zeros(nNode*2,1);
a0 = zeros(nNode*2,1);

d = d0;
v = v0;
a = a0;

Z = null(C);

% Dynamic deformation movie.  Motion JPEG AVI is broadly supported on Linux.
scriptDirectory = fileparts(mfilename('fullpath'));
resultsDirectory = fullfile(scriptDirectory,'..','..','results');
if ~exist(resultsDirectory,'dir')
    mkdir(resultsDirectory);
end

movieFile = fullfile(resultsDirectory,'quarter_ring_dynamic.avi');
videoWriter = VideoWriter(movieFile,'Motion JPEG AVI');
videoWriter.FrameRate = movieFrameRate;
open(videoWriter);

movieFigure = figure('Color','w');
movieAxes = axes('Parent',movieFigure);
hold(movieAxes,'on');
grid(movieAxes,'on');
box(movieAxes,'on');

plot(movieAxes,X(:,1),X(:,2),'k--', ...
    'LineWidth',1.5,'DisplayName','Reference');
hDeformed = plot(movieAxes,X(:,1),X(:,2),'r-', ...
    'LineWidth',2.0,'DisplayName','Deformed');
scatter(movieAxes,X(1,1),X(1,2),45,'k','filled', ...
    'DisplayName','Symmetry point');
scatter(movieAxes,X(end,1),X(end,2),45,'b','filled', ...
    'DisplayName','Load point A');

axis(movieAxes,'equal');
xlabel(movieAxes,'x');
ylabel(movieAxes,'y');
title(movieAxes,'Quarter ring dynamic response');
legend(movieAxes,'Location','best');

% Keep the camera fixed while the structure deforms.
axisMargin = 0.2 * R;
axis(movieAxes,[min(X(:,1))-axisMargin, max(X(:,1))+axisMargin, ...
    min(X(:,2))-axisMargin, max(X(:,2))+axisMargin]);
drawnow;

for timeStep = 1:nTimeStep

    loadFactor = timeStep / nTimeStep;
    FextCurrent = loadFactor * Fext;

    converged = false;
    
    % newmark predict
    dPredict = d + dt * v + (0.5 - betat) * dt^2 * a;
    vPredict = v + (1-gammat) * dt * a;
    d = dPredict;

    for iteration = 1:maximumNewtonIteration

        FmemNonlinear = zeros(nNode*2,1);
        FbendNonlinear = zeros(nNode*2,1);
        KmatMembrane = zeros(nNode*2,nNode*2);
        KgeoMembrane = zeros(nNode*2,nNode*2);
        KmatBending = zeros(nNode*2,nNode*2);
        KgeoBending = zeros(nNode*2,nNode*2);

        
        aTrial = (d - dPredict) / (betat * dt^2);
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

        Keff = M / (betat * dt^2) + Kt;

        % Current residual
        residual = Fint - FextCurrent + M * aTrial;

        % Check convergence
        projectedResidual = Z' * residual;
        relativeResidual = norm(projectedResidual) ...
            / max(1.0,norm(Z' * FextCurrent));
        constraintResidual = norm(C*d-g);

        fprintf('step %d, iteration %d, residual = %.3e\n', ...
            timeStep, iteration, relativeResidual);

        if relativeResidual < newtonTolerance ...
                && constraintResidual < newtonTolerance
            converged = true;
            break;
        end

        % Newton correction
        Kaug = [Keff, C';
                C,  zeros(size(C,1),size(C,1))];

        rhs = [-residual;
                g-C*d];

        increment = Kaug \ rhs;
        deltaD = increment(1:nNode*2);

        % Update displacement
        d = d + deltaD;

    end
    a = (d - dPredict) / (betat * dt^2);
    v = vPredict + gammat * dt * a;

    % Reconstruct nodal displacement and write a frame at the requested
    % interval.  The mechanics is still evaluated at the RKPM integration
    % points; this loop only prepares the visualization.
    if mod(timeStep,movieEvery) == 0 || timeStep == nTimeStep
        dxMovie = d(1:2:end);
        dyMovie = d(2:2:end);
        DXMovie = zeros(nNode,1);
        DYMovie = zeros(nNode,1);

        for movieNode = 1:nNode
            movieNeighbors = 1:rkpmNeighborCount(movieNode);
            movieShape = N(movieNode,movieNeighbors);
            movieGlobalIds = rkpmNeighborIds( ...
                movieNode,movieNeighbors);
            DXMovie(movieNode) = movieShape * dxMovie(movieGlobalIds);
            DYMovie(movieNode) = movieShape * dyMovie(movieGlobalIds);
        end

        XdeformedMovie = X + deformationScale * ...
            [DXMovie,DYMovie,zeros(nNode,1)];
        set(hDeformed,'XData',XdeformedMovie(:,1), ...
            'YData',XdeformedMovie(:,2));
        title(movieAxes,sprintf( ...
            'Quarter ring: t = %.4f, load factor = %.3f', ...
            timeStep*dt,loadFactor));
        drawnow;
        writeVideo(videoWriter,getframe(movieFigure));
    end


    if ~converged
        error('Current load step did not converge.');
    end

end

close(videoWriter);
fprintf('dynamic video written to: %s\n',movieFile);



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


U = [DX,DY,DZ];
Xdeformed = X + deformationScale * U;

set(hDeformed,'XData',Xdeformed(:,1),'YData',Xdeformed(:,2));
title(movieAxes,sprintf( ...
    'Quarter ring final deformation (scale = %g)',deformationScale));
drawnow;


% Analytical solution check at load point A
% The paper's t is an equivalent thickness/radius-of-gyration measure.
% For the current rectangular section, tEffective = sqrt(I/A).

% The load point is theta = pi/2, i.e. the last node, and the applied
% load is vertical. Compare its y-displacement with the benchmark value.
uA_RKPM = DY(end);

uA_exact = -totalForce * R^3 / (E * sectionInertia) ...
    * ((4 - pi) / (4*pi) - 0.25*(tEffective/R)^2);

relativeError_uA = abs(abs(uA_RKPM) - abs(uA_exact)) / abs(uA_exact);

fprintf('u_A (RKPM)  = %.8e\n', uA_RKPM);
fprintf('u_A (exact) = %.8e\n', uA_exact);
fprintf('relative error = %.4f %%\n', 100*relativeError_uA);

function vMatrix = crossMatrix(v)
if length(v) == 3
    vMatrix = [0, -v(3), v(2);
               v(3), 0, -v(1);
               -v(2), v(1), 0];
end

end



function state = evaluateCurrentConfiguration2D( ...
    dLocal, X_eta_ref, X_etaeta_ref, metric, NetaLocal, NetaetaLocal)

dLocal = dLocal(:);
NetaLocal = NetaLocal(:);
NetaetaLocal = NetaetaLocal(:);
nLocalNode = numel(NetaLocal);

Heta = zeros(2,2*nLocalNode);
Hetaeta = zeros(2,2*nLocalNode);

for localNode = 1:nLocalNode
    columns = 2*localNode-1 : 2*localNode;
    Heta(:,columns) = NetaLocal(localNode) * eye(2);
    Hetaeta(:,columns) = NetaetaLocal(localNode) * eye(2);
end

metricEta = dot(X_eta_ref, X_etaeta_ref) / metric;

H1 = Heta / metric;
H2 = Hetaeta / metric^2 ...
    - metricEta * Heta / metric^3;

% Current displacement derivatives with respect to S
uS = H1 * dLocal;
uSS = H2 * dLocal;

a_reference = X_eta_ref / metric;

c_reference = X_etaeta_ref / metric^2 ...
    - metricEta * X_eta_ref / metric^3;

a_current = a_reference + uS;
c_current = c_reference + uSS;

state.H1 = H1;
state.H2 = H2;
state.aReference = a_reference;
state.cReference = c_reference;
state.aCurrent = a_current;
state.cCurrent = c_current;

end


function [localForce, localKmat, localKgeo] = ...
    evaluateQuarterRingMembrane2D( ...
        E, sectionArea, weight, state)

% membrane contribution
H1 = state.H1;
a_current = state.aCurrent;

membraneStrain = 0.5 * ...
    (dot(a_current,a_current) - 1.0);

% Current axial force resultant
axialForce = E * sectionArea * membraneStrain;

% Nonlinear membrane B matrix
Bm = a_current' * H1;


% ---------------------------------------------------------
% 1. Internal membrane force
% ---------------------------------------------------------
localForce = ...
    Bm' * axialForce * weight;

% ---------------------------------------------------------
% 2. Material membrane tangent
% ---------------------------------------------------------
localKmat = ...
    Bm' * (E*sectionArea) * Bm * weight;

% ---------------------------------------------------------
% 3. Geometric membrane tangent
% ---------------------------------------------------------
localKgeo = ...
    axialForce * (H1' * H1) * weight;

end

function [localForce, localKmat, localKgeo] = ...
    evaluateQuarterRingBending2D( ...
        dLocal, X_eta_ref, X_etaeta_ref, metric, NetaLocal, NetaetaLocal, ...
        youngsModulus, sectionInertia, weight)

X_S = X_eta_ref / metric;
metricEta = dot(X_eta_ref, X_etaeta_ref) / metric;
X_SS = X_etaeta_ref / metric^2 ...
     - metricEta * X_eta_ref / metric^3;

% Reference normal needed by the bending strain
%   kappa = x_2' * n - X_2' * n0.
X_1_3d = [X_S; 0];
X_2_3d = [X_SS; 0];
bHat0 = cross(X_1_3d,X_2_3d);
h_b0 = norm(bHat0);
if h_b0 <= eps
    error('Degenerate reference centerline: ||bHat0|| is too small.');
end
b0 = bHat0 / h_b0;
nHat0 = cross(b0,X_1_3d);
h_n0 = norm(nHat0);
if h_n0 <= eps
    error('Degenerate reference centerline: ||nHat0|| is too small.');
end
n0 = nHat0 / h_n0;

nLocalNode = numel(NetaLocal);
Neta = zeros(2,2*nLocalNode);
Netaeta = zeros(2,2*nLocalNode);

for localNode = 1:nLocalNode
    columns = 2*localNode-1 : 2*localNode;
    Neta(:,columns) = NetaLocal(localNode) * eye(2);
    Netaeta(:,columns) = NetaetaLocal(localNode) * eye(2);
end

N1 = Neta / metric;
N2 = Netaeta / metric^2 ...
    - metricEta * Neta / metric^3;

% Current displacement derivatives with respect to S
uS = N1 * dLocal;
uSS = N2 * dLocal;

x_1 = X_S + uS;
x_2 = X_SS + uSS;

% Embed the planar centerline derivatives into R^3 for the
% handwritten binormal/frame formulas. The third displacement DOF is
% still absent; the added zero is only used for the geometric cross product.
x_1_3d = [x_1; 0];
x_2_3d = [x_2; 0];

% Unnormalized and normalized binormal.
bHat = cross(x_1_3d,x_2_3d);
h_b = norm(bHat);
if h_b <= eps
    error('Degenerate current centerline: ||bHat|| is too small.');
end
b = bHat / h_b;

% Handwritten normalization matrix A.
A = (eye(3) - b*b.') / h_b;
B1 = - A * crossMatrix(x_2_3d);
B2 = A * crossMatrix(x_1_3d);
N1_3d = [N1;zeros(1,nLocalNode*2)];
N2_3d = [N2;zeros(1,nLocalNode*2)];

% normal direction
nHat = cross(b,x_1_3d);
h_n = norm(nHat);
if h_n <= eps
    error('Degenerate current centerline: ||nHat|| is too small.');
end
n = nHat / h_n;
C = (eye(3) - n*n.') / h_n;

D1 = C * (crossMatrix(b) - crossMatrix(x_1_3d)*B1);
D2 = - C * crossMatrix(x_1_3d) * B2;

Bb = x_2_3d.' * D1 * N1_3d + x_2_3d.' * D2 * N2_3d + n.' * N2_3d;

% Bending strain, moment, and local internal force.
bendingStrain = x_2_3d.' * n - X_2_3d.' * n0;
bendingMoment = youngsModulus * sectionInertia * bendingStrain;
localForce = Bb.' * bendingMoment * weight;

% Material bending tangent.
localKmat = Bb.' * Bb * youngsModulus * sectionInertia * weight;

% Geometric bending tangent is added after Delta(Bb) is derived.
localKgeo = zeros(2*nLocalNode,2*nLocalNode);

% Derivative of A with respect to every local displacement DOF.
% delta_A(:,:,p)  = partial(A)/partial(dLocal(p))
% delta_B1(:,:,p) = partial(B1)/partial(dLocal(p))
% delta_B2(:,:,p) = partial(B2)/partial(dLocal(p))
nLocalDof = 2*nLocalNode;
delta_A = zeros(3,3,nLocalDof);
delta_B1 = zeros(3,3,nLocalDof);
delta_B2 = zeros(3,3,nLocalDof);
delta_C = zeros(3,3,nLocalDof);
delta_D1 = zeros(3,3,nLocalDof);
delta_D2 = zeros(3,3,nLocalDof);
delta_Bb = zeros(1,nLocalDof,nLocalDof);
I_eye = eye(3);

for p = 1:nLocalDof
    % Derivatives of x_1 and x_2 with respect to dLocal(p).
    N1_p = N1_3d(:,p);
    N2_p = N2_3d(:,p);

    % bHat_(i,p) = epsilon_(i,j,k) *
    %              (N1_(j,p)*x_2(k) + x_1(j)*N2_(k,p)).
    db_hat_p = zeros(3,1);
    for i = 1:3
        for j = 1:3
            for k = 1:3
                db_hat_p(i) = db_hat_p(i) ...
                    + leviCivita(i,j,k) ...
                    * (N1_p(j) * x_2_3d(k) ...
                    + x_1_3d(j) * N2_p(k));
            end
        end
    end
    % h_b,p = b_i*bHat_(i,p),  b_(i,p) = A_(i,j)*bHat_(j,p).
    hb_p = b.' * db_hat_p;
    b_p = A * db_hat_p;

    % A_(i,j,p) = partial(A_(i,j))/partial(dLocal(p)).
    for i = 1:3
        for j = 1:3
            delta_A(i,j,p) = ...
                - (b_p(i) * b(j) + b(i) * b_p(j)) / h_b ...
                - (I_eye(i,j) - b(i)*b(j)) / h_b^2 * hb_p;
        end
    end

    % B1_(i,j,p) = -A_(i,k,p)*[x_2]_(k,j)
    %               -A_(i,k)*[x_2,p]_(k,j).
    % B2_(i,j,p) =  A_(i,k,p)*[x_1]_(k,j)
    %               +A_(i,k)*[x_1,p]_(k,j).
    for i = 1:3
        for j = 1:3
            for k = 1:3
                for l = 1:3
                    % [x]_(k,j) = epsilon_(k,l,j) * x_l
                    delta_B1(i,j,p) = delta_B1(i,j,p) ...
                        - delta_A(i,k,p) * leviCivita(k,l,j) * x_2_3d(l) ...
                        - A(i,k) * leviCivita(k,l,j) * N2_p(l);

                    delta_B2(i,j,p) = delta_B2(i,j,p) ...
                        + delta_A(i,k,p) * leviCivita(k,l,j) * x_1_3d(l) ...
                        + A(i,k) * leviCivita(k,l,j) * N1_p(l);
                end
            end
        end
    end

    % nHat_(i,p) = epsilon_(i,j,k) *
    %              (b_(j,p)*x_1(k) + b(j)*N1_(k,p)).
    dn_hat_p = zeros(3,1);
    for i = 1:3
        for j = 1:3
            for k = 1:3
                dn_hat_p(i) = dn_hat_p(i) ...
                    + leviCivita(i,j,k) ...
                    * (b_p(j) * x_1_3d(k) ...
                    + b(j) * N1_p(k));
            end
        end
    end

    % h_n,p = n_i*nHat_(i,p),  n_(i,p) = C_(i,j)*nHat_(j,p).
    hn_p = n.' * dn_hat_p;
    n_p = C * dn_hat_p;

    % C_(i,j,p) = partial(C_(i,j))/partial(dLocal(p)).
    for i = 1:3
        for j = 1:3
            delta_C(i,j,p) = ...
                - (n_p(i) * n(j) + n(i) * n_p(j)) / h_n ...
                - (I_eye(i,j) - n(i)*n(j)) / h_n^2 * hn_p;
        end
    end

    % D1 = C * ( [b]_x - [x_1]_x * B1 )
    % D2 = -C * [x_1]_x * B2
    % Apply the product rule for each local DOF p.
    delta_D1(:,:,p) = ...
        delta_C(:,:,p) * ...
        (crossMatrix(b) - crossMatrix(x_1_3d) * B1) ...
        + C * ( ...
            crossMatrix(b_p) ...
            - crossMatrix(N1_p) * B1 ...
            - crossMatrix(x_1_3d) * delta_B1(:,:,p));

    delta_D2(:,:,p) = ...
        - delta_C(:,:,p) * crossMatrix(x_1_3d) * B2 ...
        - C * crossMatrix(N1_p) * B2 ...
        - C * crossMatrix(x_1_3d) * delta_B2(:,:,p);

    % Bb_(q,p) = partial(Bb(q))/partial(dLocal(p)).
    delta_Bb(:,:,p) = ...
        N2_p.' * D1 * N1_3d ...
        + x_2_3d.' * delta_D1(:,:,p) * N1_3d ...
        + N2_p.' * D2 * N2_3d ...
        + x_2_3d.' * delta_D2(:,:,p) * N2_3d ...
        + n_p.' * N2_3d;

end

% Assemble the geometric bending tangent column by column.
for p = 1:nLocalDof
    localKgeo(:,p) = bendingMoment ...
        * squeeze(delta_Bb(:,:,p)).' * weight;
end

end

function levi = leviCivita(i,j,k)
if i==1 && j ==2 && k == 3
    delta = 1;
elseif i == 2 && j == 3 && k == 1
    delta = 1;
elseif i == 3 && j == 1 && k == 2
    delta = 1;
elseif i == 1 && j == 3 && k == 2
    delta = -1;
elseif i == 3 && j == 2 && k == 1
    delta = -1;
elseif i == 2 && j == 1 && k == 3
    delta = -1;
else
    delta = 0;
end

levi = delta;

end
