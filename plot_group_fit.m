% function plot_group_fit(expnr,mtype,nSteps)
% 
% Fit specified model to data of specified experiment and plots the results.
%
% mtype: 1=IL, 2=SA, 3=EP, 4=VP

function plot_group_fit(expnr,mtype,nSteps)

if ~exist('nSteps','var')
    nSteps = 15;
end

expinfo = getExperimentInfo(expnr);

for ii=1:length(expinfo.subjids)
    
    fname = ['saved_results/exp' num2str(expnr) '/results_exp' num2str(expnr) '_' num2str(ii) '_' num2str(nSteps) '_' num2str(mtype) '.mat'];
    if mtype==1 && ~exist(fname,'file')
        fit_IL_model(expnr,ii,nSteps);
    elseif mtype==2 && ~exist(fname,'file')
        fit_SA_model(expnr,ii,nSteps);
    elseif (mtype==3 || mtype==4)  && ~exist(fname,'file')
        fit_EPVP_model(expnr,ii,nSteps);        
    end
    load(fname,'w_emp','w_fit','csd_emp','csd_fit','L','fitpars','uN');
   
    all_w_emp(ii,:) = w_emp;
    all_w_fit(ii,:) = w_fit;
    all_csd_emp(ii,:) = csd_emp;
    all_csd_fit(ii,:) = csd_fit;
    all_L(ii) = L;
    
    all_fitpars(ii,:) = fitpars(:);
end

save(['saved_results/exp' num2str(expnr) '/results_' num2str(expnr) '_group_' num2str(nSteps) '_' num2str(mtype) '.mat'],'all_fitpars','all_L','mtype');

% compute means and std errors
w_emp_mean = mean(all_w_emp);
w_emp_stderr = std(all_w_emp)/sqrt(length(expinfo.subjids));
w_fit_mean = mean(all_w_fit);
w_fit_stderr = std(all_w_fit)/sqrt(length(expinfo.subjids));

csd_emp_mean = mean(all_csd_emp);
csd_emp_stderr = std(all_csd_emp)/sqrt(length(expinfo.subjids));
csd_fit_mean = mean(all_csd_fit);
csd_fit_stderr = std(all_csd_fit)/sqrt(length(expinfo.subjids));

% plot summary stats
mnames = {'IL','SA','EP','VP'};

% close
figure
set(gcf,'Position',get(gcf,'Position').*[1 1 1 .6]);
set(gcf,'PaperPosition',get(gcf,'PaperPosition').*[1 1 1 .6]);

subplot(1,2,1);
hold on;
uN=uN(:)';
XX=[uN uN(end:-1:1)];
YY=[w_fit_mean+w_fit_stderr w_fit_mean(end:-1:1)-w_fit_stderr(end:-1:1)];
if length(XX)>2
    patch(XX,YY,[.7 .7 .7],'LineStyle','none');
else
    errorbar(uN,w_fit_mean,w_fit_stderr,'rx');
end
errorbar(uN,w_emp_mean,w_emp_stderr,'ko');
ylim([0 1]);
xlim([min(uN)-0.5 max(uN)+0.5]);
R2(1) = 1 - sum((w_emp_mean(:)-w_fit_mean(:)).^2)/sum((w_emp_mean(:)-mean(w_emp_mean)).^2);
R2(2) = 1 - sum((csd_emp_mean(:)-csd_fit_mean(:)).^2)/sum((csd_emp_mean(:)-mean(csd_emp_mean)).^2);
text(min(uN),0.1,['RMSE=' num2str(sqrt(mean((all_w_emp(:)-all_w_fit(:)).^2)),4)]);
text(min(uN),0.2,['R^2=' num2str(R2(1),4)]);
title(['Experiment ' num2str(expnr) ', ' mnames{mtype}]);
set(gca,'Xtick',uN);
xlabel('Set size');
ylabel('w')

subplot(1,2,2);
hold on;
XX=[uN uN(end:-1:1)];
YY=[csd_fit_mean+csd_fit_stderr csd_fit_mean(end:-1:1)-csd_fit_stderr(end:-1:1)];
if length(XX)>2
    patch(XX,YY,[.7 .7 .7],'LineStyle','none');
else
    errorbar(uN,csd_fit_mean,csd_fit_stderr,'rx');
end
errorbar(uN,csd_emp_mean,csd_emp_stderr,'ko');
ylim([0.15 .6]);
xlim([min(uN)-0.5 max(uN)+0.5]);
text(min(uN),0.50,['RMSE=' num2str(sqrt(mean((all_csd_emp(:)-all_csd_fit(:)).^2)),4)]);
text(min(uN),0.55,['R^2=' num2str(R2(2),4)]);
set(gca,'Xtick',uN);
xlabel('Set size');
ylabel('CSD')
