function plot_BMC_results(expnr,nSteps) 

if ~exist('nSteps','var')
    nSteps = 15;
end

ei = getExperimentInfo(expnr);
mnames = {'IL','SA','EP','VP'};

% get model likelihoods
for ii=1:length(ei.subjids)
    L_mat(ii,1) = fit_IL_model(expnr,ii,nSteps,0);
    L_mat(ii,2) = fit_SA_model(expnr,ii,nSteps,0);
    [A B] = fit_EPVP_model(expnr,ii,nSteps,0);
    L_mat(ii,3:4) = [A B];
end

% express model likelihoods relative to best model
L_mean = mean(L_mat);
winidx = find(L_mean==max(L_mean));
L_mat = bsxfun(@minus,L_mat,L_mat(:,winidx));

% compute mean and std.err.
L_mean = mean(L_mat);
L_stderr = std(L_mat)/sqrt(size(L_mat,1));

% plot
figure
bar(1:4,L_mean,'k');
hold on
errorbar(1:4,L_mean,L_stderr,'k','LineStyle','none');
set(gca,'XTick',1:4,'XTickLabel',mnames);
ylabel('Relative model log likelihood');
