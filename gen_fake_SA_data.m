% function [error_vec N_vec] = gen_fake_SA_data(modelpars,nTrials,N)
%
% Produces synthetic data created using the SA model

function [error_vec N_vec] = gen_fake_SA_data(modelpars,nTrials,N,kmap,Jmap)
% parse model pars
J1 = modelpars(1);
km = modelpars(2);
K  = modelpars(3);

% precompute mapping kappa <--> J
if ~exist('kmap','var')
    kmap = linspace(0,700,1e5);
    Jmap = kmap.*besseli(1,kmap)./besseli(0,kmap);
end

for ii=1:nTrials
    
    % decide which stimuli to encode
    if K>=N
        p_high = mod(K,N)/N;
        if rand<p_high
            J_high = min(max(Jmap),J1*(floor(K/N)+1));
            kappa_high = interp1(Jmap,kmap,J_high);
            x(ii) = circ_vmrnd(0,kappa_high); % encode with kappa_high
        else
            J_low  = min(max(Jmap),J1*(floor(K/N)));
            kappa_low = interp1(Jmap,kmap,J_low);
            x(ii) = circ_vmrnd(0,kappa_low); % encode with kappa_low
        end
    else
        Pmem = K/N;
        if rand<Pmem
            kappa1 = interp1(Jmap,kmap,J1);
            x(ii) = circ_vmrnd(0,kappa1); % encode with kappa corresponding to 1 chunk of resource
        else
            x(ii) = rand*2*pi-pi;                              % random guess
        end
    end    
end

error_vec = circ_vmrnd(x,km);  % add response noise
N_vec = ones(1,nTrials)*N;


