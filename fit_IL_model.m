% function fit_IL_model(expnr,subjidx,nSteps)
%
% Fits the IL model to the data of the specified subject.
%
% nSteps indicates the number of bins used to discretize each parameter
% dimension

function L = fit_IL_model(expnr,subjidx,nSteps,makeplot)

if ~exist('nSteps','var')
    nSteps=15;
end
if ~exist('makeplot','var')
    makeplot = 1;
end

% check if not done yet
results_fname = ['saved_results/exp' num2str(expnr) '/results_exp' num2str(expnr) '_' num2str(subjidx) '_' num2str(nSteps) '_1.mat'];
if ~exist(results_fname,'file');
    Kpar_range = 1:20;
    
    % load subject data
    [error_vec N_vec] = readdata(expnr,subjidx);
    error_vec = abs(error_vec);
    uN = unique(N_vec);

    % create tables if not exist yet
    for ii=1:length(uN)
        create_IL_table(uN(ii),nSteps);
    end    
    
    % load tables
    fprintf('Loading tables...\n');
    for ii=1:max(N_vec)
        load(['precomputed_tables/table_IL_' num2str(ii) '_' num2str(nSteps) '.mat'],'error_table','km_range','error_range');
        all_error_table{ii} = shiftdim(error_table,1);   % shiftdim, to avoid having to use squeeze thousands of times later on
    end
    
    % find estimation error indices for each set size
    for ii=1:length(uN)
        error_idx{uN(ii)} = interp1(error_range,1:length(error_range),error_vec(N_vec==uN(ii)),'nearest','extrap'); % error w.r.t. target
    end
    
    % initialize all_LLH matrix (for models without tau and/or Kpar, these dimensions will be singleton)
    all_LLH = zeros(length(km_range),length(Kpar_range));
    
    fprintf('Computing p(data | IL model, parameters) for all parameter combinations...\n');
    for jj=1:length(km_range)
        for kk=1:length(Kpar_range)
            for ll=1:length(uN)
                N = uN(ll);
                K = Kpar_range(kk);
                if K<N
                    pvec = (1-K/N)*1/2/pi + K/N*all_error_table{N}(:,jj);
                else
                    pvec = all_error_table{N}(:,jj);
                end
                p_resp = pvec(error_idx{N});
                p_resp = max(p_resp,eps);  % avoid log(0)
                all_LLH(jj,kk) = all_LLH(jj,kk) + sum(log(p_resp));
                
            end
        end
    end
    
    %-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
    % find maximum and compute model likelihood %
    %-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
    [maxLLH I] = max(all_LLH(:));
    [idx1 idx2] = ind2sub(size(all_LLH),I);
    km_hat = km_range(idx1);
    Kpar_hat = Kpar_range(idx2);
    fitpars = [km_hat Kpar_hat];
    parnames = {'km','K'};
    
    % compute model likelihood (integrated of all parameters)
    M = exp(all_LLH-maxLLH);
    L = log(trapz(km_range,trapz(Kpar_range,M,2),1)) + maxLLH - log(range(km_range)*range(Kpar_range));
    
    %-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
    % compute predicted w and CSD %
    %-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
    fprintf('Computing summary statistics for subject and IL model\n')
    % estimate descriptive mixture fit parameters for data and for model
    for ii=1:length(uN)
        % dm fit to subject data
        dm_fitpars_emp(ii,:) = dmfit(error_vec(N_vec==uN(ii)));
        % dm fit to model data
        nTrials = sum(N_vec==uN(ii));
        for jj=1:10
            error_vec_sim = gen_fake_IL_data(fitpars,nTrials,uN(ii));
            fitpars_tmp(jj,:) = dmfit(error_vec_sim);
        end
        csd_fit(ii) = mean(sqrt(1 - besseli(1,fitpars_tmp(:,2))./besseli(0,fitpars_tmp(:,2))));
        w_fit(ii) = mean(1-fitpars_tmp(:,1));
    end
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
    xlabel('Set size');
    ylabel('w');
    title(['IL, subjidx=' num2str(subjidx) ', L=' num2str(L)]);
    subplot(1,2,2);
    plot(uN,csd_emp,'ko');
    hold on;
    plot(uN,csd_fit,'r-');
    ylim([0 1]);
    xlim([0.5 8.5]);        
    xlabel('Set size');
    ylabel('CSD');
    drawnow()
end

