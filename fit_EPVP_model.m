% function fit_EPVP_model(expnr,subjidx,nSteps)
%
% Fits the EP and VP models to the data of the specified subject.
%
% nSteps indicates the number of bins used to discretize each parameter
% dimension

function [L_EP L_VP] = fit_EPVP_model(expnr,subjidx,nSteps,makeplot)

if ~exist('nSteps','var')
    nSteps=15;
end
if ~exist('makeplot','var')
    makeplot = 1;
end

% check if results not exist yet
results_fname_3 = ['saved_results/exp' num2str(expnr) '/results_exp' num2str(expnr) '_' num2str(subjidx) '_' num2str(nSteps) '_3.mat'];
results_fname_4 = ['saved_results/exp' num2str(expnr) '/results_exp' num2str(expnr) '_' num2str(subjidx) '_' num2str(nSteps) '_4.mat'];
if ~(exist(results_fname_3,'file') && exist(results_fname_4,'file'))
    
    % precompute mapping between J and kappa
    kmap = linspace(0,700,1e5);
    Jmap = kmap.*besseli(1,kmap)./besseli(0,kmap);
    
    % load subject data
    [error_vec N_vec] = readdata(expnr,subjidx);
    error_vec = abs(error_vec); % model predictions only depend on error magnitude, not sign
    uN = unique(N_vec);

    % create tables if not exist yet
    for ii=1:length(uN)
        create_EPVP_table(uN(ii),nSteps);
    end
    
    if max(uN)>8
        error('Code assumes that 8 is the maximum set size. Need to change Poisson weights for data sets with N>8');
    end
    
    % load tables
    fprintf('Loading tables...\n');
    for ii=1:max(N_vec)
        load(['precomputed_tables/table_EPVP_' num2str(ii) '_' num2str(nSteps) '.mat'],'error_table','J1bar_range','tau_range','km_range','error_range','power_range');
        all_error_table{ii} = shiftdim(error_table,4);   % shiftdim, to avoid having to use squeeze thousands of times later on
    end
    
    % find estimation error indices for each set size
    for ii=1:length(uN)
        error_idx{uN(ii)} = interp1(error_range,1:length(error_range),error_vec(N_vec==uN(ii)),'nearest','extrap'); % error w.r.t. target
    end
    
    % initialize all_LLH matrix (for models without tau and/or Kpar, these dimensions will be singleton)
    all_LLH = zeros(length(J1bar_range),length(power_range),length(tau_range),length(km_range));
    itTotal = numel(all_LLH)*length(uN);
    itCnt=0;
    
    tic;
    fprintf('Computing p(data | EP/VP model, parameters) for all parameter combinations...\n');
    for ii=1:length(J1bar_range)
        for jj=1:length(power_range)
            for kk=1:length(tau_range)
                for ll=1:length(km_range)
                    for mm=1:length(uN)
                        N = uN(mm);
                        pvec = all_error_table{N}(:,ii,jj,kk,ll); % get model predictions for current (J1bar,power,tau,km,N)-combination
                        p_resp = pvec(error_idx{N}); % get p(response | data) for all trials
                        p_resp = max(p_resp,eps); % avoid log(0)
                        all_LLH(ii,jj,kk,ll) = all_LLH(ii,jj,kk,ll) + sum(log(p_resp));
                        itCnt=itCnt+1;
                    end
                end
            end
        end
        fprintf('ETL=%2.1f minutes\n',(toc/itCnt)*(itTotal-itCnt)/60);
    end
    
    %---- VP ---%
    
    %-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
    % find maximum and compute model likelihood %
    %-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
    [maxLLH I] = max(all_LLH(:));
    [idx1 idx2 idx3 idx4] = ind2sub(size(all_LLH),I);
    J1bar_hat = J1bar_range(idx1);
    power_hat = power_range(idx2);
    tau_hat = tau_range(idx3);
    km_hat = km_range(idx4);
    fitpars = [J1bar_hat power_hat tau_hat km_hat];
    parnames = {'J1bar','power','tau','km'};
    
    % compute model likelihood (integrated of all parameters)
    M = exp(all_LLH-maxLLH);
    L = log(trapz(J1bar_range,trapz(power_range,trapz(tau_range,trapz(km_range,M,4),3),2),1)) + maxLLH - log(range(J1bar_range)*range(power_range)*range(tau_range)*range(km_range));
    L_VP = L;
    
    %-%-%-%-%-%-%-%-%-%-%-%-%-%-%-
    % compute summary statistics %
    %-%-%-%-%-%-%-%-%-%-%-%-%-%-%-
    % estimate descriptive mixture fit parameters for data and for model
    fprintf('Computing summary statistics for subject and VP model.')
    for ii=1:length(uN)
        fprintf('.');
        % dm fit to subject data
        dm_fitpars_emp(ii,:) = dmfit(error_vec(N_vec==uN(ii)));
        % dm fit to model data
        nTrials = sum(N_vec==uN(ii));
        for jj=1:10
            error_vec_sim = gen_fake_EPVP_data(fitpars,4,nTrials,uN(ii),kmap,Jmap);
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
    save(results_fname_4,'L','fitpars','w_emp','w_fit','csd_emp','csd_fit','uN','csd_fit','parnames');
    
    %---- EP ---%
    if tau_range(1) ~= 0
        error('First value in tau_range should be 0');
    end
    all_LLH = squeeze(all_LLH(:,:,1,:)); % EP is special case of VP (namely tau==0)
    
    %-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
    % find maximum and compute model likelihood %
    %-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
    [maxLLH I] = max(all_LLH(:));
    [idx1 idx2 idx3] = ind2sub(size(all_LLH),I);
    J1bar_hat = J1bar_range(idx1);
    power_hat = power_range(idx2);
    km_hat = km_range(idx3);
    fitpars = [J1bar_hat power_hat 0 km_hat];
    
    % compute model likelihood (integrated of all parameters)
    M = exp(all_LLH-maxLLH);
    L = log(trapz(J1bar_range,trapz(power_range,trapz(km_range,M,3),2),1)) + maxLLH - log(range(J1bar_range)*range(power_range)*range(km_range));
    L_EP = L;
    
    %-%-%-%-%-%-%-%-%-%-%-%-%-%-%-
    % compute summary statistics %
    %-%-%-%-%-%-%-%-%-%-%-%-%-%-%-
    % estimate descriptive mixture fit parameters for data and for model
    fprintf('Computing summary statistics for subject and EP model.')
    for ii=1:length(uN)
        fprintf('.');
        % dm fit to subject data
        dm_fitpars_emp(ii,:) = dmfit(error_vec(N_vec==uN(ii)));
        % dm fit to model data
        nTrials = sum(N_vec==uN(ii));
        for jj=1:10
            error_vec_sim = gen_fake_EPVP_data(fitpars,3,nTrials,uN(ii),kmap,Jmap);
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
    save(results_fname_3,'L','fitpars','w_emp','w_fit','csd_emp','csd_fit','uN','csd_fit','parnames');
else
    load(results_fname_3,'L');
    L_EP = L;
    load(results_fname_4,'L');
    L_VP = L;    
end

%-%-%-%-%-%-%-%-
% plot results %
%-%-%-%-%-%-%-%-
if makeplot
    close
    figure
    set(gcf,'Position',get(gcf,'Position').*[1 1 1 0.6]);
    set(gcf,'PaperPosition',get(gcf,'PaperPosition').*[1 1 1 0.6]);
    load(results_fname_3)
    subplot(2,2,1);
    plot(uN,w_emp,'ko');
    hold on;
    plot(uN,w_fit,'r-');
    ylim([0 1]);
    xlim([0.5 8.5]);
    title(['EP, subjidx=' num2str(subjidx) ', L=' num2str(L)]);
    subplot(2,2,2);
    plot(uN,csd_emp,'ko');
    hold on;
    plot(uN,csd_fit,'r-');
    ylim([0 1]);
    xlim([0.5 8.5]);
    
    load(results_fname_4)
    subplot(2,2,3);
    plot(uN,w_emp,'ko');
    hold on;
    plot(uN,w_fit,'r-');
    ylim([0 1]);
    xlim([0.5 8.5]);
    title(['VP, subjidx=' num2str(subjidx) ', L=' num2str(L)]);
    subplot(2,2,4);
    plot(uN,csd_emp,'ko');
    hold on;
    plot(uN,csd_fit,'r-');
    ylim([0 1]);
    xlim([0.5 8.5]);
    
    drawnow()
end

