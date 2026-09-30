%% W39 snapshot: 2D geometrically nonlinear membrane bar scaffold
% The active learning file is:
%   code/matlab/nonlinear/static_2d/main_nonlinear_membrane_bar_learning.m
%
% This weekly file records the 2026-09-20 starting point.  Continue editing
% the active file rather than this snapshot.

clear; clc;

L = 1.0;
E = 1.0e8;
radius = 0.05;
sectionArea = pi*radius^2;
EA = E*sectionArea;

nNode = 17;
order = 2;
supportScale = 3.2;
nGaussPerCell = 2;

SNode = linspace(0,L,nNode).';
X = [SNode,zeros(nNode,1)];
h = L/(nNode-1);
support = supportScale*h*ones(nNode,1);
nDof = 2*nNode;
d = zeros(nDof,1);

totalLoad = 0.01*EA;
nLoadStep = 5;

% Nonlinear learning order:
%   1. strain = 0.5*(a.'*a-1)
%   2. Bm = a.'*H
%   3. Fint = integral(Bm.'*EA*strain)
%   4. Kmat = integral(Bm.'*EA*Bm)
%   5. Kgeo = integral(EA*strain*H.'*H)
%   6. R = Fint-loadFactor*Fext
%   7. Newton update and analytical stretch comparison

fprintf('W39 nonlinear membrane scaffold initialized: %d nodes, %d DOFs.\n', ...
    nNode,nDof);
