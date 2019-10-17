% function fit_SA_model(expnr,subjidx,ktype,nSteps)
%
% Fits the SA model to the data of the specified subject.
%
% nSteps indicates the number of bins used to discretize each parameter
% dimension

function L = fit_SA_model(expnr,subjidx,nSteps,makeplot)

if ~exist('nSteps','var')
    nSteps=15;
end
if ~exist('makeplot','var')
    makeplot = 1;
end

% check if not done yet
results_fname = ['saved_results/exp' num2str(expnr) '/results_exp' num2str(expnr) '_' num2str(subjidx) '_' num2str(nSteps) '_2.mat'];
if ~exist(results_fname,'file')    
    % precompute mapping between J and kappa
    kmap = linspace(0,700,1e5);
    Jmap = kmap.*besseli(1,kmap)./besseli(0,kmap);
    
    Kpar_range = 1:8;
    
    % load subject data
    [error_vec N_vec] = readdata(expnr,subjidx);
    error_vec = abs(error_vec);
    uN = unique(N_vec);
    
    % create tables if not exist yet
    for ii=1:length(uN)
        create_SA_table(uN(ii),nSteps);
    end
    
    % load tables
    fprintf('Loading tables...\n');
    for ii=1:max(N_vec)
        load(['precomputed_tables/table_SA_' num2str(ii) '_' num2str(nSteps) '.mat'],'error_table','J1_range','K_range','km_range','error_range');
        all_error_table{ii} = shiftdim(error_table,3);   % shiftdim, to avoid having to use squeeze thousands of times later on
    end
    
    % find estimation error indices for each set size
    for ii=1:length(uN)
        error_idx{uN(ii)} = interp1(error_range,1:length(error_range),error_vec(N_vec==uN(ii)),'nearest','extrap'); % error w.r.t. target
    end
    
    % initialize all_LLH matrix (for models without tau and/or Kpar, these dimensions will be singleton)
    all_LLH = zeros(length(J1_range),length(km_range),length(Kpar_range));
    itTotal = numel(all_LLH)*length(uN);
    itCnt=0;
    tic;
    fprintf('Computing p(data | SA model, parameters) for all parameter combinations...\n');
    for ii=1:length(J1_range)
        for jj=1:length(km_range)
            for ll=1:length(uN)
                N = uN(ll);
                
                for kk=1:length(Kpar_range)
                    K = Kpar_range(kk);
                    if K<=N
                        pvec = (1-K/N)*1/2/pi + K/N*all_error_table{N}(:,N,ii,jj); % predictions for encoding with 1 chunk of resource
                    else
                        pvec = all_error_table{N}(:,K,ii,jj);
                    end
                    p_resp = pvec(error_idx{N});
                    p_resp = max(p_resp,eps);
                    all_LLH(ii,jj,kk) = all_LLH(ii,jj,kk) + sum(log(p_resp));
                    itCnt=itCnt+1;
                end
            end
        end
    end
    
    %-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
    % find maximum and compute model likelihood %
    %-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
    [maxLLH I] = max(all_LLH(:));
    [idx1 idx2 idx3] = ind2sub(size(all_LLH),I);
    J1_hat = J1_range(idx1);
    km_hat = km_range(idx2);
    Kpar_hat = Kpar_range(idx3);
    fitpars = [J1_hat km_hat Kpar_hat];
    parnames = {'J1','km','K'};
    
    % compute model likelihood (integrated of all parameters)
    M = exp(all_LLH-maxLLH);
    L = log(trapz(J1_range,trapz(km_range,trapz(Kpar_range,M,3),2),1)) + maxLLH - log(range(J1_range)*range(km_range)*range(Kpar_range));
    
    %-%-%-%-%-%-%-%-%-%-%-%-%-%-%-
    % compute summary statistics %
    %-%-%-%-%-%-%-%-%-%-%-%-%-%-%-
    % estimate descriptive mixture fit parameters for data and for model
    fprintf('Computing summary statistics for subject and SA model.')
    for ii=1:length(uN)
        fprintf('.');
        % dm fit to subject data
        dm_fitpars_emp(ii,:) = dmfit(error_vec(N_vec==uN(ii)));
        % dm fit to model data
        nTrials = sum(N_vec==uN(ii));
        for jj=1:10
            error_vec_sim = gen_fake_SA_data(fitpars,nTrials,uN(ii),kmap,Jmap);
            fitpars_tmp(jj,:) = dmfit(error_vec_sim);
        end
        csd_fit(ii) = mean(sqrt(1 - besseli(1,fitpars_tmp(:,2))./besseli(0,fitpars_tmp(:,2))));
        w_fit(ii) = mean(1-fitpars_tmp(:,1));
    end
    fprintf('\n');
    w_emp = 1-dm_fitpars_emp(:,1);
    csd_emp = sqrt(1 - besseli(1,dm_fitpars_emp(:,2))./besseli(0,dm_fitpars_emp(:,2)));
    
    %-%-%-%-%-%-%-%-
    % save results %
    %-%-%-%-%-%-%-%-
    save(results_fname,'L','fitpars','w_emp','w_fit','csd_emp','csd_fit','uN','parnames');    
else
    load(results_fname);
end

%-%-%-%-%-%-%-%-
% plot results %
%-%-%-%-%-%-%-%-
if makeplot
    close
    figure
    set(gcf,'Position',get(gcf,'Position').*[1 1 1 0.6]);
    set(gcf,'PaperPosition',get(gcf,'PaperPosition').*[1 1 1 0.6]);
    subplot(1,2,1);
    plot(uN,w_emp,'ko');
    hold on;
    plot(uN,w_fit,'r-');
    ylim([0 1]);
    xlim([0.5 8.5]);
    title(['SA, subjidx=' num2str(subjidx) ', L=' num2str(L)]);
    subplot(1,2,2);
    plot(uN,csd_emp,'ko');
    hold on;
    plot(uN,csd_fit,'r-');
    ylim([0 1]);
    xlim([0.5 8.5]);        
    drawnow()
end
