% function create_EPVP_table(N,nSteps)
%
% Creates a table with predictions for continuous resource model with
% J(bar) a power law function of N
%
% INPUT
%   N       : set size for which to create a prediction table
%   nSteps  : the number of bins to discretize each parameter dimension

function create_EPVP_table(N,nSteps)

% precompute mapping kappa <--> J
kmap = linspace(0,700,1e5);
Jmap = kmap.*besseli(1,kmap)./besseli(0,kmap);

% initalize ranges
J1bar_range = logspace(log10(1),log10(500),nSteps);
power_range = linspace(-4,0,nSteps);
tau_range   = [0 logspace(log10(1),log10(250),nSteps-1)];
km_range    = logspace(log10(2),log10(500),nSteps);
error_range = linspace(0,pi,91);                       % discretization of the error space
error_range = error_range(1:end-1)+diff(error_range(1:2))/2;
nSamples    = 32; % number of samples to use to approximate VP predictions (was 256 for the results in paper, but to speed things up it was set to 32 here)

error_table = zeros(length(J1bar_range),length(power_range),length(tau_range),length(km_range),length(error_range));

% precompute I_0(km) for all values of km
bessel_km = besseli(0,km_range);

% compute and store predictions
fname = ['precomputed_tables/table_EPVP_' num2str(N) '_' num2str(nSteps) '.mat'];
if ~exist(fname,'file')
    fprintf('Creating table for continuous resource model, N=%d\n',N);
    itCnt=0;
    tic;
    for ii=1:length(J1bar_range)    
        for jj=1:length(power_range)
            Jbar = J1bar_range(ii)*N.^power_range(jj);
            for kk=1:length(tau_range)
                % draw precision values and compute corresponding values of kappa
                if tau_range(kk)<1e-10
                    J_sample = Jbar;
                else
                    J_sample = gamrnd(Jbar/tau_range(kk),tau_range(kk),1,nSamples);
                    J_sample = min(J_sample,max(Jmap));
                end
                kappa_sample = interp1(Jmap,kmap,J_sample);
                kappa_sample_sq = kappa_sample.^2;
                denom_tmp = 2*pi*besseli(0,kappa_sample);
                
                % loop over motor noise (recycle kappa's across motor noise values)
                for ll=1:length(km_range)
                    denom = denom_tmp*bessel_km(ll);
                    k_c = zeros(length(error_range),nSamples);
                    for mm=1:length(error_range)
                        k_c(mm,:) = sqrt(km_range(ll)^2 + kappa_sample_sq + 2*km_range(ll)*kappa_sample*cos(error_range(mm)));
                    end
                    k_c = min(k_c,700);
                    error_table(ii,jj,kk,ll,:) = mean(bsxfun(@rdivide,besseli(0,k_c),denom),2);
                    itCnt=itCnt+1;
                end
            end
        end
        fprintf('ETL=%2.1f minutes\n',(toc/itCnt)*(numel(error_table)/length(error_range)-itCnt)/60);
    end
    save(fname,'error_table','J1bar_range','power_range','tau_range','km_range','error_range','nSamples');
end

