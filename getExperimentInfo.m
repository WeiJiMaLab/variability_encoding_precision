% function info = getExperimentInfo(expnr)
% 
% Returns a structure with information about a data set

function info = getExperimentInfo(expnr)

if expnr==1
    info.datadir = 'data/vandenberg2012_color/';  
    info.subjids = {'cc','clc','ela','hml','jv','kw','mbc','mt','rjj','ss','stp','wc','wjm'};    
    info.expname = {'Van den berg et al. 2012, color (scrolling)'};
elseif expnr==2
    info.datadir = 'data/vandenberg2012_orientation/';
    info.subjids = {'AA','ACO','ELA','RGG','TCS','WJM'};
    info.expname = {'Van den berg et al. 2012, orientation'};
elseif expnr==3
    % To add a data set, modify the following:
    info.datadir = 'data/add_your_directory_here/';
    info.subjids = {'1','2','etc'};
    info.expname = {'Just a string to describe the data set'};
    % Note that you should also adjust readdata.m 
end
