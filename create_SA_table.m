% function create_SA_table(N,nSteps)
%
% Creates a table with SA model predictions for set size N.
%
% nSteps indicates the number of bins to discretize each parameter dimension

function create_SA_table(N,nSteps)

% precompute mapping kappa <--> J
kmap = linspace(0,700,1e5);
Jmap = kmap.*besseli(1,kmap)./besseli(0,kmap);

% initalize ranges 
J1_range    = logspace(log10(1),log10(40),nSteps);           % (mean) precision provided by 1 slot
K_range     = 1:8;                                           % number of available slots to encode items 
km_range    = logspace(log10(5),log10(500),nSteps);    % kappa of the motor noise distribution
error_range = linspace(0,pi,91);                             % discretization of the error space
error_range = error_range(1:end-1)+diff(error_range(1:2))/2;

error_table = zeros(length(K_range),length(J1_range),length(km_range),length(error_range));

% precompute I_0(km) for all values of km
bessel_km = besseli(0,km_range);

% compute and store predictions
fname = ['precomputed_tables/table_SA_' num2str(N) '_' num2str(nSteps) '.mat'];
if ~exist(fname,'file')
    fprintf('Creating table for SA model, N=%d\n',N);
    itCnt=0;    
    tic;    
    for ii=1:length(K_range)
        for jj=1:length(J1_range)
            k1 = interp1(Jmap,kmap,J1_range(jj));
            bessel_k1 = besseli(0,k1);            
            for kk=1:length(km_range)                
                if N<=K_range(ii)
                    % Case N<=K
                    J_high = min(max(Jmap),J1_range(jj)*(floor(K_range(ii)/N)+1));
                    J_low = min(max(Jmap),J1_range(jj)*floor(K_range(ii)/N));
                    kappa_high = interp1(Jmap,kmap,J_high);
                    kappa_low = interp1(Jmap,kmap,J_low);
                    p_high = mod(K_range(ii),N)/N;
                    kc_high = min(699,sqrt(kappa_high^2 + km_range(kk)^2 + 2*kappa_high*km_range(kk)*cos(error_range)));
                    kc_low = min(699,sqrt(kappa_low^2 + km_range(kk)^2 + 2*kappa_low*km_range(kk)*cos(error_range)));
                    bessel_ratio_high = besseli(0,kc_high)./bessel_km(kk);
                    bessel_ratio_low = besseli(0,kc_low)./bessel_km(kk);
                    error_table(ii,jj,kk,:) = p_high * bessel_ratio_high .* 1./(2*pi*besseli(0,kappa_high)) + (1-p_high) * bessel_ratio_low .* 1./(2*pi*besseli(0,kappa_low));
                else
                    % Case N>K
                    p_guess = 1-K_range(ii)/N;
                    kc = min(max(Jmap),sqrt(k1^2 + km_range(kk)^2 + 2*k1*km_range(kk)*cos(error_range)));
                    error_table(ii,jj,kk,:) = p_guess*(1/2/pi) + (1-p_guess)*besseli(0,kc)./(2*pi*bessel_k1*bessel_km(kk));
                end                               
                itCnt=itCnt+1;
            end            
        end
        fprintf('ETL=%2.1f minutes\n',(toc/itCnt)*(numel(error_table)/length(error_range)-itCnt)/60);
    end
    save(fname,'error_table','J1_range','K_range','km_range','error_range');
end

