% function [error_vec N_vec] = gen_fake_IL_data(modelpars,nTrials,N)
%
% Produces synthetic data created with the IL model.
function [error_vec N_vec] = gen_fake_IL_data(modelpars,nTrials,N)

% parse model pars
km = modelpars(1);
K  = modelpars(2);

% simulate
for ii=1:nTrials
    % decide which stimuli to encode
    if K>=N
        x(ii) = 0; % item was remembered; estimation is perfect
    else
        p_mem = K/N;
        if rand<p_mem
            x(ii) = 0; % item was remembered; estimation is perfect
        else
            x(ii) = rand*2*pi-pi; % item was not remembered; random guess
        end
    end
end

error_vec = circ_vmrnd(x,km); % add response noise
N_vec = ones(1,nTrials)*N;
