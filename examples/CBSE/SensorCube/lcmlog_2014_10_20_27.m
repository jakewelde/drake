function [d imFnames]=lcmlog_2014_10_20_27()
full_fname = '/home/gizatt/Downloads/lcmlog_2014_10_20_27.mat';
fname = '/home/gizatt/Downloads/lcmlog_2014_10_20_27.mat';
if (exist(full_fname,'file'))
    filename = full_fname;
else
    filename = fname;
end
d = load(filename);
