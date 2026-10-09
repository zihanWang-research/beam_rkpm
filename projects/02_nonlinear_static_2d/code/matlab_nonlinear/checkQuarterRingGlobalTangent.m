%% Finite-difference check for the assembled quarter-ring tangent
% The main program first computes a converged nonlinear state. This script
% then differentiates the assembled internal force with respect to every
% global displacement degree of freedom.

clear; clc;

mainFile = fullfile(fileparts(mfilename('fullpath')), ...
    'main_quarter_ring_nonlinear_learning.m');
run(mainFile);

dBase = d;
Kanalytic = KtTrial;

nDof = numel(dBase);
KfiniteDifference = zeros(nDof,nDof);
relativeStep = 1.0e-6;

for column = 1:nDof
    direction = zeros(nDof,1);
    direction(column) = 1.0;
    step = relativeStep * max(1.0,abs(dBase(column)));

    forcePlus = assembleQuarterRingInternalForce( ...
        dBase + step*direction, nNode, rkpmNeighborIds, ...
        rkpmNeighborCount, Neta, Netaeta, Xeta, Xetaeta, metric, ...
        E, sectionArea, sectionInertia, nodalArcWeight);
    forceMinus = assembleQuarterRingInternalForce( ...
        dBase - step*direction, nNode, rkpmNeighborIds, ...
        rkpmNeighborCount, Neta, Netaeta, Xeta, Xetaeta, metric, ...
        E, sectionArea, sectionInertia, nodalArcWeight);

    KfiniteDifference(:,column) = ...
        (forcePlus - forceMinus) / (2.0*step);
end

relativeError = norm(Kanalytic-KfiniteDifference,'fro') / ...
    max(1.0,norm(KfiniteDifference,'fro'));

fprintf('\nAssembled tangent finite-difference check\n');
fprintf('||Kanalytic||_F = %.6e\n', norm(Kanalytic,'fro'));
fprintf('||KFD||_F       = %.6e\n', norm(KfiniteDifference,'fro'));
fprintf('relative error   = %.6e\n', relativeError);

function internalForce = assembleQuarterRingInternalForce( ...
    displacement, nNode, rkpmNeighborIds, rkpmNeighborCount, ...
    Neta, Netaeta, Xeta, Xetaeta, metric, youngsModulus, sectionArea, ...
    sectionInertia, nodalArcWeight)

internalForce = zeros(2*nNode,1);

for node = 1:nNode
    localNodeIds = rkpmNeighborIds( ...
        node,1:rkpmNeighborCount(node))';
    nLocalNode = numel(localNodeIds);
    localDofIds = reshape([2*localNodeIds-1, 2*localNodeIds]',[],1);
    dLocal = displacement(localDofIds);

    NetaLocal = Neta(node,1:rkpmNeighborCount(node))';
    NetaetaLocal = Netaeta(node,1:rkpmNeighborCount(node))';
    XetaRef = Xeta(node,1:2)';
    XetaetaRef = Xetaeta(node,1:2)';

    Heta = zeros(2,2*nLocalNode);
    Hetaeta = zeros(2,2*nLocalNode);
    for localNode = 1:nLocalNode
        columns = 2*localNode-1:2*localNode;
        Heta(:,columns) = NetaLocal(localNode)*eye(2);
        Hetaeta(:,columns) = NetaetaLocal(localNode)*eye(2);
    end

    metricEta = dot(XetaRef,XetaetaRef)/metric(node);
    H1 = Heta/metric(node);
    H2 = Hetaeta/metric(node)^2 ...
        - metricEta*Heta/metric(node)^3;
    aReference = XetaRef/metric(node);
    aCurrent = aReference + H1*dLocal;

    membraneStrain = 0.5*(dot(aCurrent,aCurrent)-1.0);
    axialForce = youngsModulus*sectionArea*membraneStrain;
    Bm = aCurrent'*H1;
    localMembraneForce = Bm'*axialForce*nodalArcWeight(node);

    [localBendingForce,~,~] = evaluateQuarterRingBending2D( ...
        dLocal, XetaRef, XetaetaRef, metric(node), NetaLocal, ...
        NetaetaLocal, youngsModulus, sectionInertia, nodalArcWeight(node));

    internalForce(localDofIds) = internalForce(localDofIds) ...
        + localMembraneForce + localBendingForce;
end

end
