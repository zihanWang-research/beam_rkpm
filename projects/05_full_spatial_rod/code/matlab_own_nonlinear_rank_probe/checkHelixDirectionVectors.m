function report = checkHelixDirectionVectors(varargin)
%CHECKHELIXDIRECTIONVECTORS Validate the PCA/RKPM frame of a 3-D helix.
%
% This diagnostic is independent of the nonlinear solver.  It checks the
% three directions used by the centerline geometry:
%   t : tangent, n : normal, b : binormal.
%
% The PCA calculation supplies only tPCA.  The reference frame tRef, nRef,
% bRef is reconstructed from X_,S and X_,SS, exactly as in the mechanics
% code.  The latter is compared with the analytical helix frame.
%
% Usage from this directory:
%   report = checkHelixDirectionVectors();
%   report = checkHelixDirectionVectors(99,1.0,0.20,3*pi,3.5,2.5,2,false);
%
% Inputs (all optional)
%   nNode, R, pitch, thetaMax, alphaPCA, alphaRKPM, rkpmOrder, showPlot
%
% The function does not modify the nonlinear main program.

if nargin < 1, nNode = 99; else, nNode = varargin{1}; end
if nargin < 2, R = 1.0; else, R = varargin{2}; end
if nargin < 3, pitch = 0.20; else, pitch = varargin{3}; end
if nargin < 4, thetaMax = 3*pi; else, thetaMax = varargin{4}; end
if nargin < 5, alphaPCA = 3.5; else, alphaPCA = varargin{5}; end
if nargin < 6, alphaRKPM = 2.5; else, alphaRKPM = varargin{6}; end
if nargin < 7, rkpmOrder = 2; else, rkpmOrder = varargin{7}; end
if nargin < 8, showPlot = false; else, showPlot = varargin{8}; end

if nNode < 5 || nNode ~= floor(nNode)
    error('nNode must be an integer greater than or equal to 5.');
end

%% Reference helix
thetaNode = linspace(0,thetaMax,nNode).';
X = [R*cos(thetaNode), R*sin(thetaNode), pitch*thetaNode];
dX = diff(X,1,1);
segmentLength = vecnorm(dX,2,2);
h = mean(segmentLength);

%% PCA tangent domains
supportPCA = alphaPCA*h;
pcaNeighborIds = cell(nNode,1);
pcaEigenvalues = zeros(nNode,3);
tPCA = zeros(nNode,3);

for i = 1:nNode
    distance = vecnorm(X-X(i,:),2,2);
    pcaIds = find(distance <= supportPCA).';
    if numel(pcaIds) < rkpmOrder + 1
        error('PCA domain at node %d has too few neighbors.',i);
    end
    pcaNeighborIds{i} = pcaIds;

    XbarLocal = mean(X(pcaIds,:),1);
    YLocal = X(pcaIds,:) - XbarLocal;
    covariancePCA = YLocal.'*YLocal;
    [pcaVectors,pcaValues] = eig(covariancePCA);
    [sortedValues,sortIds] = sort(diag(pcaValues),'descend');
    pcaEigenvalues(i,:) = sortedValues.';
    tPCA(i,:) = pcaVectors(:,sortIds(1)).';
end

% Eigenvectors have arbitrary signs.  Orient the PCA tangent consistently.
for i = 2:nNode
    if dot(tPCA(i,:),tPCA(i-1,:)) < 0
        tPCA(i,:) = -tPCA(i,:);
    end
end
startEdge = X(2,:) - X(1,:);
if dot(tPCA(1,:),startEdge) < 0
    tPCA = -tPCA;
end

%% RKPM reconstruction in the local PCA coordinate
NCell = cell(nNode,1);
NetaCell = cell(nNode,1);
NetaetaCell = cell(nNode,1);
rkpmNeighborIds = cell(nNode,1);
etaCell = cell(nNode,1);

Xeta = zeros(nNode,3);
Xetaeta = zeros(nNode,3);
metric = zeros(nNode,1);
partitionError = zeros(nNode,1);
firstDerivativeSum = zeros(nNode,1);
secondDerivativeSum = zeros(nNode,1);

supportRKPM = alphaRKPM*h;
for i = 1:nNode
    distance = vecnorm(X-X(i,:),2,2);
    ids = find(distance <= supportRKPM).';
    eta = (X(ids,:) - X(i,:))*tPCA(i,:).';

    [active,NLocal,NetaLocal,NetaetaLocal] = ...
        rkpm_shape_1d(0,eta,supportRKPM*ones(size(eta)),rkpmOrder);
    ids = ids(active);
    eta = eta(active);

    NCell{i} = NLocal(:);
    NetaCell{i} = NetaLocal(:);
    NetaetaCell{i} = NetaetaLocal(:);
    rkpmNeighborIds{i} = ids(:);
    etaCell{i} = eta(:);

    Xeta(i,:) = NetaLocal(:).'*X(ids,:);
    Xetaeta(i,:) = NetaetaLocal(:).'*X(ids,:);
    metric(i) = norm(Xeta(i,:));

    partitionError(i) = abs(sum(NLocal)-1.0);
    firstDerivativeSum(i) = abs(sum(NetaLocal));
    secondDerivativeSum(i) = abs(sum(NetaetaLocal));
end

if any(metric <= eps)
    error('The reconstructed reference tangent has zero metric.');
end

% Convert derivatives from the local PCA coordinate eta to arc length S.
metricEta = sum(Xeta.*Xetaeta,2)./metric;
X_S = Xeta./metric;
X_SS = Xetaeta./metric.^2 ...
    - Xeta.*metricEta./metric.^3;

%% Reference frame reconstructed from the numerical geometry
tRef = X_S./vecnorm(X_S,2,2);
bRef = cross(X_S,X_SS,2);
bRef = bRef./vecnorm(bRef,2,2);
nRef = cross(bRef,tRef,2);

%% Analytical helix frame
gamma = sqrt(R^2 + pitch^2);
tExact = [-R*sin(thetaNode), R*cos(thetaNode), pitch*ones(nNode,1)]/gamma;
nExact = [-cos(thetaNode), -sin(thetaNode), zeros(nNode,1)];
bExact = [pitch*sin(thetaNode), -pitch*cos(thetaNode), R*ones(nNode,1)]/gamma;
curvatureExact = R/gamma^2;
XSSExact = curvatureExact*nExact;

%% Errors and frame checks
interior = 2:nNode-1;  % endpoint support domains are one-sided

report.nNode = nNode;
report.meanChordLength = h;
report.pcaEigenvalues = pcaEigenvalues;
report.metric = metric;
report.tPCA = tPCA;
report.tRef = tRef;
report.nRef = nRef;
report.bRef = bRef;
report.tExact = tExact;
report.nExact = nExact;
report.bExact = bExact;
report.X_S = X_S;
report.X_SS = X_SS;
report.rkpmNeighborIds = rkpmNeighborIds;
report.eta = etaCell;

report.maxPcaTangentError = max(1-abs(sum(tPCA(interior,:).*tExact(interior,:),2)));
report.maxTangentError = max(vecnorm(tRef(interior,:)-tExact(interior,:),2,2));
report.maxNormalError = max(vecnorm(nRef(interior,:)-nExact(interior,:),2,2));
report.maxBinormalError = max(vecnorm(bRef(interior,:)-bExact(interior,:),2,2));
report.maxTangentDerivativeError = max(vecnorm(X_S(interior,:)-tExact(interior,:),2,2));
report.maxSecondDerivativeError = max(vecnorm(X_SS(interior,:)-XSSExact(interior,:),2,2));

frameDotTN = sum(tRef.*nRef,2);
frameDotTB = sum(tRef.*bRef,2);
frameDotNB = sum(nRef.*bRef,2);
frameRightHanded = cross(tRef,nRef,2)-bRef;
report.maxOrthonormalityError = max(abs([frameDotTN;frameDotTB;frameDotNB]));
report.maxRightHandedError = max(vecnorm(frameRightHanded,2,2));
report.minTangentContinuity = min(sum(tRef(1:end-1,:).*tRef(2:end,:),2));
report.minNormalContinuity = min(sum(nRef(1:end-1,:).*nRef(2:end,:),2));
report.minBinormalContinuity = min(sum(bRef(1:end-1,:).*bRef(2:end,:),2));

report.maxPartitionError = max(partitionError);
report.maxFirstDerivativeSum = max(firstDerivativeSum);
report.maxSecondDerivativeSum = max(secondDerivativeSum);
report.maxNormalPlaneEigenvalueRatio = ...
    max(pcaEigenvalues(interior,2)./max(pcaEigenvalues(interior,1),eps));
report.maxThirdEigenvalueRatio = ...
    max(pcaEigenvalues(interior,3)./max(pcaEigenvalues(interior,1),eps));

fprintf('\n3-D helix direction-vector check\n');
fprintf('nodes = %d, interior nodes = %d\n',nNode,numel(interior));
fprintf('max PCA tangent 1-|dot| error      = %.6e\n',report.maxPcaTangentError);
fprintf('max reconstructed t error           = %.6e\n',report.maxTangentError);
fprintf('max reconstructed n error           = %.6e\n',report.maxNormalError);
fprintf('max reconstructed b error           = %.6e\n',report.maxBinormalError);
fprintf('max X_S error                       = %.6e\n',report.maxTangentDerivativeError);
fprintf('max X_SS error                      = %.6e\n',report.maxSecondDerivativeError);
fprintf('max frame orthonormality error      = %.6e\n',report.maxOrthonormalityError);
fprintf('max frame right-handedness error    = %.6e\n',report.maxRightHandedError);
fprintf('minimum adjacent t/n/b dot products = [%.6e %.6e %.6e]\n', ...
    report.minTangentContinuity,report.minNormalContinuity,report.minBinormalContinuity);
fprintf('RKPM [sum(N)-1, sum(Neta), sum(Netaeta)] = [%.6e %.6e %.6e]\n', ...
    report.maxPartitionError,report.maxFirstDerivativeSum,report.maxSecondDerivativeSum);
fprintf('max PCA eigenvalue ratios [lambda2/lambda1, lambda3/lambda1] = [%.6e %.6e]\n', ...
    report.maxNormalPlaneEigenvalueRatio,report.maxThirdEigenvalueRatio);

if showPlot
    figure('Name','Helix direction-vector check');
    plot3(X(:,1),X(:,2),X(:,3),'k.-','LineWidth',1.0); hold on;
    stride = max(1,round(nNode/20));
    quiver3(X(1:stride:end,1),X(1:stride:end,2),X(1:stride:end,3), ...
        tRef(1:stride:end,1),tRef(1:stride:end,2),tRef(1:stride:end,3), ...
        0.25,'r','LineWidth',1.0);
    quiver3(X(1:stride:end,1),X(1:stride:end,2),X(1:stride:end,3), ...
        nRef(1:stride:end,1),nRef(1:stride:end,2),nRef(1:stride:end,3), ...
        0.25,'g','LineWidth',1.0);
    quiver3(X(1:stride:end,1),X(1:stride:end,2),X(1:stride:end,3), ...
        bRef(1:stride:end,1),bRef(1:stride:end,2),bRef(1:stride:end,3), ...
        0.25,'b','LineWidth',1.0);
    axis equal; grid on; view(35,25);
    xlabel('X'); ylabel('Y'); zlabel('Z');
    legend('helix','t','n','b','Location','best');
    title('Reconstructed reference triad');
end
end
