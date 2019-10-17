% function create_IL_table(N,nSteps)
%
% Creates a table with IL model predictions for set size N.
%
% nSteps indicates the number of bins to discretize each parameter dimension

function create_IL_table(N,nSteps)

% initalize ranges
km_range    = logspace(log10(1),log10(50),nSteps);     % kappa of the motor noise distribution
error_range = linspace(0,pi,91);                       % discretization of the error space
error_range = error_range(1:end-1)+diff(error_range(1:2))/2;

error_table = zeros(length(km_range),length(error_range));

% precompute I_0(km) for all values of km
bessel_km = besseli(0,km_range);

% compute and store predictions
fname = ['precomputed_tables/table_IL_' num2str(N) '_' num2str(nSteps) '.mat'];
if ~exist(fname,'file')
    fprintf('Creating table for IL model, N=%d\n',N);
    for jj=1:length(km_range)
        error_table(jj,:) = 1/(2*pi*bessel_km(jj)) * exp(km_range(jj)*cos(error_range));
    end
    save(fname,'error_table','km_range','error_range');
end

