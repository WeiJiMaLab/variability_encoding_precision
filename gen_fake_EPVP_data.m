% function [error_vec N_vec] = gen_fake_EPVP_data(modelpars,modelnr,nTrials,N)
%
% Produces synthetic data created with the EP or VP model (modelnr 3 = EP,
% modelnr 4=VP).
function [error_vec N_vec] = gen_fake_EPVP_data(modelpars,modelnr,nTrials,N,kmap,Jmap)

% modelnr: 3=EP, 4=VP

% parse model parameters
J1bar = modelpars(1);
power = modelpars(2);
tau   = modelpars(3);
km    = modelpars(4);

if ~exist('kmap','var')
    % precompute mapping between J and kappa
    kmap = linspace(0,700,1e5);
    Jmap = kmap.*besseli(1,kmap)./besseli(0,kmap);
end    

% simulate trials
Jbar = J1bar*N^power;
if tau==0 || modelnr==3
    J = ones(1,nTrials)*Jbar;  % EP
else
    J = gamrnd(Jbar/tau,tau,1,nTrials); % VP
end
J = min(J,max(Jmap));
kappa = interp1(Jmap,kmap,J);
x = circ_vmrnd(0,kappa);  % estimation errors
error_vec = circ_vmrnd(x,km); % add response noise
N_vec = ones(1,nTrials)*N;