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
