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
