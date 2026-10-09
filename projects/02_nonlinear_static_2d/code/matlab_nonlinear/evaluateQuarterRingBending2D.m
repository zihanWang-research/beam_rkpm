%% Local nonlinear bending force and tangent for the RKPM quarter ring
% This file contains the same manuscript-based implementation used by
% main_quarter_ring_nonlinear_learning.m. It is separated so that the
% bending tangent can be checked independently.

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

function vMatrix = crossMatrix(v)
if length(v) == 3
    vMatrix = [0, -v(3), v(2);
               v(3), 0, -v(1);
               -v(2), v(1), 0];
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
