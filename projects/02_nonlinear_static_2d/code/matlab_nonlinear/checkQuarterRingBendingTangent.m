%% Finite-difference check for the local nonlinear bending tangent
% This check compares (Kmat + Kgeo)*p with the directional derivative of
% the bending internal force. It is intentionally independent of the load
% stepping loop in main_quarter_ring_nonlinear_learning.m.

clear; clc;

radius = 1.0;
theta = pi/4;
metric = radius;
X_eta_ref = [-radius*sin(theta); radius*cos(theta)];
X_etaeta_ref = [-radius*cos(theta); -radius*sin(theta)];

% A deterministic local RKPM derivative pattern for this tangent check.
NetaLocal = [-0.40; -0.20; 0.00; 0.20; 0.40];
NetaetaLocal = [0.30; -0.10; -0.40; -0.10; 0.30];
nLocalDof = 2*numel(NetaLocal);

youngsModulus = 1.0e6;
sectionInertia = 2.0e-4;
weight = 1.0;

dLocal = 1.0e-3 * [ ...
     0.40; -0.20; -0.10;  0.30;  0.20; -0.10; ...
    -0.30;  0.10;  0.15; -0.25];

[force0, Kmat, Kgeo] = evaluateQuarterRingBending2D( ...
    dLocal, X_eta_ref, X_etaeta_ref, metric, NetaLocal, NetaetaLocal, ...
    youngsModulus, sectionInertia, weight);

Kanalytic = Kmat + Kgeo;
direction = (1:nLocalDof)';
direction = direction / norm(direction);

fprintf('Bending tangent finite-difference check\n');
fprintf('||force|| = %.6e\n', norm(force0));
fprintf('||Kmat||  = %.6e\n', norm(Kmat,'fro'));
fprintf('||Kgeo||  = %.6e\n', norm(Kgeo,'fro'));

for step = [1e-3, 1e-4, 1e-5, 1e-6]
    [forcePlus,~,~] = evaluateQuarterRingBending2D( ...
        dLocal + step*direction, X_eta_ref, X_etaeta_ref, metric, ...
        NetaLocal, NetaetaLocal, youngsModulus, sectionInertia, weight);
    [forceMinus,~,~] = evaluateQuarterRingBending2D( ...
        dLocal - step*direction, X_eta_ref, X_etaeta_ref, metric, ...
        NetaLocal, NetaetaLocal, youngsModulus, sectionInertia, weight);

    finiteDifference = (forcePlus - forceMinus) / (2*step);
    analytic = Kanalytic * direction;
    relativeError = norm(analytic - finiteDifference) / ...
        max(1.0, norm(finiteDifference));

    fprintf('h = %.1e, relative error = %.6e\n', step, relativeError);
end
