function [active, N, N_eta, N_etaeta] = rkpm_shape_1d(eta, node_eta, support, order)
%RKPM_SHAPE_1D One-dimensional RKPM shape functions and derivatives.
%
% Inputs
%   eta      evaluation coordinate in the local PCA parameter space
%   node_eta coordinates of candidate control points in that same space
%   support  RKPM support radius for every candidate control point
%   order    polynomial reproduction order: 1, 2, or 3
%
% Outputs are restricted to active candidate indices.  If this function is
% called with a locally filtered candidate list, active indexes that list,
% not the original global control-point list.

node_eta = node_eta(:);
support = support(:);

if numel(node_eta) ~= numel(support)
    error('node_eta and support must have the same length.');
end
if any(support <= 0)
    error('All support radii must be positive.');
end
if order < 1 || order > 3
    error('This function supports order = 1, 2, or 3.');
end

p_count = order + 1;
dx_all = eta - node_eta;
active = find(abs(dx_all) <= support);
if numel(active) < p_count
    error('Not enough active nodes at eta = %.6g. Need %d, got %d.', ...
        eta, p_count, numel(active));
end

dx = dx_all(active);
active_eta = node_eta(active);
active_support = support(active);
[w, w_eta, w_etaeta] = cubic_spline_weight_1d(dx, active_support);

P = polynomial_basis_1d(active_eta, order);
p = polynomial_basis_1d(eta, order);
p_eta = polynomial_basis_eta_1d(eta, order);
p_etaeta = polynomial_basis_etaeta_1d(eta, order);

W = diag(w);
W_eta = diag(w_eta);
W_etaeta = diag(w_etaeta);

A = P * W * P.';
A_eta = P * W_eta * P.';
A_etaeta = P * W_etaeta * P.';
B = P * W;
B_eta = P * W_eta;
B_etaeta = P * W_etaeta;

gamma = A \ p;
gamma_eta = A \ (p_eta - A_eta * gamma);
gamma_etaeta = A \ ( ...
    p_etaeta - 2.0 * A_eta * gamma_eta - A_etaeta * gamma);

N = (gamma.' * B).';
N_eta = (gamma_eta.' * B + gamma.' * B_eta).';
N_etaeta = ( ...
    gamma_etaeta.' * B ...
    + 2.0 * gamma_eta.' * B_eta ...
    + gamma.' * B_etaeta).';
end


function [w, w_eta, w_etaeta] = cubic_spline_weight_1d(dx, support)
q = abs(dx) ./ support;
direction = sign(dx);
w = zeros(size(dx));
w_eta = zeros(size(dx));
w_etaeta = zeros(size(dx));

first_half = q <= 0.5;
second_half = q > 0.5 & q <= 1.0;
q1 = q(first_half);
q2 = q(second_half);

w(first_half) = 2.0/3.0 - 4.0*q1.^2 + 4.0*q1.^3;
w(second_half) = 4.0/3.0 - 4.0*q2 + 4.0*q2.^2 - 4.0/3.0*q2.^3;
w_eta(first_half) = (-8.0*q1 + 12.0*q1.^2) ...
    .* direction(first_half) ./ support(first_half);
w_eta(second_half) = (-4.0 + 8.0*q2 - 4.0*q2.^2) ...
    .* direction(second_half) ./ support(second_half);
w_etaeta(first_half) = (-8.0 + 24.0*q1) ./ support(first_half).^2;
w_etaeta(second_half) = (8.0 - 8.0*q2) ./ support(second_half).^2;
end


function p = polynomial_basis_1d(x, order)
x = x(:).';
if order == 1
    p = [ones(size(x)); x];
elseif order == 2
    p = [ones(size(x)); x; x.^2];
elseif order == 3
    p = [ones(size(x)); x; x.^2; x.^3];
else
    error('order must be 1, 2, or 3.');
end
if numel(x) == 1
    p = p(:);
end
end


function p_eta = polynomial_basis_eta_1d(x, order)
if order == 1
    p_eta = [0; 1];
elseif order == 2
    p_eta = [0; 1; 2.0*x];
elseif order == 3
    p_eta = [0; 1; 2.0*x; 3.0*x^2];
else
    error('order must be 1, 2, or 3.');
end
end


function p_etaeta = polynomial_basis_etaeta_1d(x, order)
if order == 1
    p_etaeta = [0; 0];
elseif order == 2
    p_etaeta = [0; 0; 2.0];
elseif order == 3
    p_etaeta = [0; 0; 2.0; 6.0*x];
else
    error('order must be 1, 2, or 3.');
end
end
