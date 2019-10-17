% function [error_vec N_vec] = readdata(expnr,subjidx)
%
% Returns the data of a particular subject in two vectors:
%  error_vec : vector with estimation errors
%  N_vec     : vector with corresponding set sizes
%
% If you added your own data to getExperimentInfo.m, you should add code
% here to read those data from file.

function [error_vec N_vec] = readdata(expnr,subjidx)

% get information about file location and subject identifiers
expinfo = getExperimentInfo(expnr);

if expnr == 1
    %-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
    % Color experiment in Van den Berg et al. 2012 (PNAS) - response by scrolling %
    %-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
    files = dir([expinfo.datadir expinfo.subjids{subjidx} '*.mat']);
    fprintf('Subject %s, reading %d files\n',expinfo.subjids{subjidx},length(files));
    target_vec = [];
    N_vec = [];
    nontarget_vec = {};
    resp_vec = [];
    for ii=1:length(files)
        load([expinfo.datadir files(ii).name]);
        if ABBA_array(1)==1
            resp_vec = [resp_vec recording1.sub_picked];
            for jj=1:length(recording1.set_size)
                N_vec(end+1) = recording1.set_size(jj);
                target_vec(end+1) = recording1.c_array(recording1.test(jj),jj);
            end
        end        
        if ABBA_array(2)==1
            resp_vec = [resp_vec recording2.sub_picked];
            for jj=1:length(recording2.set_size)
                N_vec(end+1) = recording2.set_size(jj);
                target_vec(end+1) = recording2.c_array(recording2.test(jj),jj);
            end
        end        
        if ABBA_array(3)==1
            resp_vec = [resp_vec recording3.sub_picked];
            for jj=1:length(recording3.set_size)
                N_vec(end+1) = recording3.set_size(jj);
                target_vec(end+1) = recording3.c_array(recording3.test(jj),jj);
            end
        end        
        if ABBA_array(4)==1
            resp_vec = [resp_vec recording4.sub_picked];
            for jj=1:length(recording4.set_size)
                N_vec(end+1) = recording4.set_size(jj);
                target_vec(end+1) = recording4.c_array(recording4.test(jj),jj);
            end
        end
    end    
    % convert [1,180] to [-pi,pi] range 
    target_vec = target_vec/180*2*pi-pi; 
    resp_vec = resp_vec/180*2*pi-pi;
    error_vec = circ_dist(resp_vec,target_vec);
    
elseif expnr == 2
    %-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
    % Orientation experiment in Van den Berg et al. 2012 (PNAS) %
    %-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%-%
    files = dir([expinfo.datadir expinfo.subjids{subjidx} '*.mat']);    
    fprintf('Subject %s, reading %d files\n',expinfo.subjids{subjidx},length(files));
    N_vec = [];
    resp_vec = [];
    nontarget_vec = {};    
    target_vec = [];    
    for ii=1:length(files)
        load([expinfo.datadir files(ii).name]);
        for jj=1:5
            if jj==1
                data = block_one;
            elseif jj==2
                data = block_two;
            elseif jj==3
                data = block_three;
            elseif jj==4
                data = block_four;
            elseif jj==5
                data = block_five;
            end
            for kk=1:length(data.set_size)
                N_vec(end+1) = data.set_size(kk);
                resp_vec(end+1) = data.subject_responses(kk);
                target_vec(end+1) = data.true_values(kk);
            end
        end
    end
    
    % convert [1,180] to [-pi,pi] range 
    target_vec = target_vec/180*2*pi-pi;
    resp_vec = resp_vec/180*2*pi-pi;  
    error_vec = circ_dist(resp_vec,target_vec);
elseif expnr == 3
    % Add code here to return the data for the specified subject
    N_vec = [];
    error_vec = [];    
end


