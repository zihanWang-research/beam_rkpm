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

nNode = 17;             % Replace with a value such as 17 after TODO 1.1.
R = 1;
E = 1e6;
totalForce = 1;         % Point-load magnitude F on the complete circular ring.

sectionArea = 0.2;

slenderness = 100;    % R / t
tEffective = R / slenderness;
sectionHeight = 2*sqrt(3)*tEffective;
sectionWidth = sectionArea / sectionHeight;

sectionInertia = sectionWidth * sectionHeight^3 / 12;

Dm = zeros(nNode,1);
Db = zeros(nNode,1);
nodalArcWeight = zeros(nNode,1);
rkpmOrder = 2;
pointsPerCell = 1;
alphaPCA = 3.5;
alphaRKPM = 2.5;  % Dimensionless alpha in rho = alpha*h.


% Keep the unknown ordering fixed throughout the project:
% d = [u_x(1); u_y(1); u_x(2); u_y(2); ... ; u_x(nNode); u_y(nNode)]
thetaNode = linspace(0,pi/2,nNode)';
X = zeros(nNode,3);
X(:,1) = R * cos(thetaNode);
X(:,2) = R * sin(thetaNode);
plot3(X(:,1),X(:,2),X(:,3))

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

% rkpm part diff pca part
maxNeighbor = ceil(2*alphaRKPM) + 1;
N = zeros(nNode,maxNeighbor);
Neta = zeros(nNode,maxNeighbor);
Netaeta = zeros(nNode,maxNeighbor);
rkpmActiveLocalIds = zeros(nNode,maxNeighbor);
rkpmNeighborIds = zeros(nNode,maxNeighbor);
rkpmNeighborCount = zeros(nNode,1);

supportRKPM = alphaRKPM * h * ones(1,nNode);

for i = 1:nNode
    rkpmCount = 0;
    for j = 1:nNode
        distance = norm(X(i,:) - X(j,:));
        if distance <= supportRKPM(j)
            rkpmCount = rkpmCount + 1;
            rkpmNeighborIds(i,rkpmCount) = j;
        end
    end
    rkpmNeighborCount(i) = rkpmCount;
end

etaRKPM = zeros(size(rkpmNeighborIds));
for i = 1:nNode
    for j = 1:rkpmNeighborCount(i)
        etaRKPM(i,j) = (X(rkpmNeighborIds(i,j),:) - X(i,:)) * tPCA(i,:)';
    end
end

for i = 1:nNode
    rkpmNodeLocal = etaRKPM(i,1:rkpmNeighborCount(i));
    supportRKPM = alphaRKPM * h * ones(1,rkpmNeighborCount(i));
    [active,NLocal,NetaLocal,NetaetaLocal] = ...
        rkpm_shape_1d(0,rkpmNodeLocal,supportRKPM,rkpmOrder);
    N(i,1:rkpmNeighborCount(i)) = NLocal';
    Neta(i,1:rkpmNeighborCount(i)) = NetaLocal';
    Netaeta(i,1:rkpmNeighborCount(i)) = NetaetaLocal';
    rkpmActiveLocalIds(i,1:rkpmNeighborCount(i)) = active';
end

% Reconstructed reference geometry derivatives and their metric quantities.
Xeta = zeros(nNode,3);
Xetaeta = zeros(nNode,3);
metric = zeros(nNode,1);
XetaetaNorm = zeros(nNode,1);

for i = 1:nNode
    Xeta(i,:) = Neta(i,1:rkpmNeighborCount(i)) ...
        * X(rkpmNeighborIds(i,1:rkpmNeighborCount(i))',:);
    Xetaeta(i,:) = Netaeta(i,1:rkpmNeighborCount(i)) ...
        * X(rkpmNeighborIds(i,1:rkpmNeighborCount(i))',:);
    metric(i) = norm(Xeta(i,:));
    XetaetaNorm(i) = norm(Xetaeta(i,:));
end

% construction the 3 frame
tRef   = zeros(nNode,3);
bRef   = zeros(nNode,3);
nRef   = zeros(nNode,3);
rRef   = zeros(nNode,3);
hRef   = zeros(nNode,3);
for i = 1:nNode
    tRef(i,:) = Xeta(i,:) / metric(i);
    rRef(i,:) = cross(Xeta(i,:),Xetaeta(i,:));
    bRef(i,:) = rRef(i,:) / norm(rRef(i,:));
    hRef(i,:) = cross(bRef(i,:),Xeta(i,:));
    nRef(i,:) = hRef(i,:) / norm(hRef(i,:));
end

% delta b or n
Pb = zeros(3,3,nNode);
Pn = zeros(3,3,nNode);

for i = 1:nNode
    Pb(:,:,i) = eye(3) - bRef(i,:)' * bRef(i,:);
    Pb(:,:,i) = Pb(:,:,i) / norm(rRef(i,:));
    Pn(:,:,i) = eye(3) - nRef(i,:)' * nRef(i,:);
    Pn(:,:,i) = Pn(:,:,i) / norm(hRef(i,:));
end

A1frame = zeros(3,3,nNode);
A2frame = zeros(3,3,nNode);
B1frame = zeros(3,3,nNode);
B2frame = zeros(3,3,nNode);

for i = 1:nNode
    aCross = crossMatrix(Xeta(i,:));
    cCross = crossMatrix(Xetaeta(i,:));
    bCross = crossMatrix(bRef(i,:));
    A1frame(:,:,i) = -Pn(:,:,i) * aCross;
    A2frame(:,:,i) = Pn(:,:,i) * bCross;
    B1frame(:,:,i) = -Pb(:,:,i) * cCross;
    B2frame(:,:,i) = Pb(:,:,i) * aCross;
end

M1frame = zeros(3,3,nNode);
M2frame = zeros(3,3,nNode);
for i = 1:nNode
    M1frame(:,:,i) = A1frame(:,:,i) * B1frame(:,:,i) + A2frame(:,:,i);
    M2frame(:,:,i) = A1frame(:,:,i) * B2frame(:,:,i);
end

Kmem = zeros(nNode*2);
Dm = E * sectionArea * metric.^4;
Db = E * sectionInertia * metric.^4;

arcSegmentLength = (thetaNode(2:nNode) - thetaNode(1:nNode-1)) * R;
nodalArcWeight(2:nNode-1) = arcSegmentLength(1:nNode-2)/2 ...
    + arcSegmentLength(2:nNode-1)/2;
nodalArcWeight(1) = arcSegmentLength(1)/2;
nodalArcWeight(nNode) = arcSegmentLength(end)/2;

for i = 1:nNode
    local_ids = rkpmNeighborIds(i,1:rkpmNeighborCount(i))';
    xieta = Xeta(i,1:2);
    nieta = Neta(i,1:rkpmNeighborCount(i))';
    Bm_local = nieta * xieta;
    Bm_local = reshape(Bm_local',[],1);
    bm_ids = [local_ids*2-1,local_ids*2];
    bm_ids = reshape(bm_ids',[],1);
    Kmem(bm_ids,bm_ids) = Kmem(bm_ids,bm_ids) + Bm_local*Bm_local'*Dm(i)*nodalArcWeight(i);
end
Kbend = zeros(nNode*2);
for i = 1:nNode
    local_ids = rkpmNeighborIds(i,1:rkpmNeighborCount(i))';
    xietaeta = Xetaeta(i,1:2);
    nieta = Neta(i,1:rkpmNeighborCount(i))';
    nietaeta = Netaeta(i,1:rkpmNeighborCount(i))';
    Bb_local = nietaeta * (nRef(i,1:2) + xietaeta * M2frame(1:2,1:2,i)) + nieta * (xietaeta * M1frame(1:2,1:2,i));
    Bb_local = reshape(Bb_local',[],1);
    bb_ids = [local_ids*2-1,local_ids*2];
    bb_ids = reshape(bb_ids',[],1);
    Kbend(bb_ids,bb_ids) = Kbend(bb_ids,bb_ids) + Bb_local*Bb_local'*Db(i)*nodalArcWeight(i);
end
K = Kbend + Kmem;

% BC I = 1  u_1y = 0 u_1xeta = 0   I = nnode u_nnodex = 0 u_nnodeyeta = 0
C = zeros(4,nNode*2);
bcNodalBegin = 1;
NetaNodalBegin = Neta(bcNodalBegin,1:rkpmNeighborCount(bcNodalBegin));
NNodalBegin = N(bcNodalBegin,1:rkpmNeighborCount(bcNodalBegin));
C(1,rkpmNeighborIds(bcNodalBegin,1:rkpmNeighborCount(bcNodalBegin))*2) = NNodalBegin;
C(2,rkpmNeighborIds(bcNodalBegin,1:rkpmNeighborCount(bcNodalBegin))*2-1) = NetaNodalBegin;
bcNodalEnd = nNode;

NetaNodalEnd = Neta(bcNodalEnd,1:rkpmNeighborCount(bcNodalEnd));
NNodalEnd = N(bcNodalEnd,1:rkpmNeighborCount(bcNodalEnd));
C(3,rkpmNeighborIds(bcNodalEnd,1:rkpmNeighborCount(bcNodalEnd))*2-1) = NNodalEnd;
C(4,rkpmNeighborIds(bcNodalEnd,1:rkpmNeighborCount(bcNodalEnd))*2) = NetaNodalEnd;

g = zeros(4,1);

% force BC
F = zeros(nNode*2,1);
quarterRingLoad = totalForce / 2;
F(nNode*2) = -quarterRingLoad;

Kaug = [K,C';
        C,zeros(4,4)];
raug = [F;
        g];
daug = Kaug \ raug;
d = daug(1:nNode*2);
dx = d(1:2:end);
dy = d(2:2:end);
dz = zeros(nNode,1);
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

uA_RKPM = DX(1);

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
