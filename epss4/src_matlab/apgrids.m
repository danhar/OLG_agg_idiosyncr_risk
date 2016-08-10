% Copyright (C) 2016 Daniel Harenberg - All Rights Reserved
global gplot_all
%clear;
close all;


if (isempty(gplot_all))    
    visibility = 'on ';
    dir = 'FD1_KK/msge';  % 'ge'
else
    visibility = 'off';
    dir = getenv('EPSSDIR');
end

cd('../model_output/')
% jr = dlmread('params.txt','', [40 11 40 13]);

cd(dir)
apgrid=dlmread('apgrid_mean.txt');
scrsz = get(0,'ScreenSize');

for j=[1,35]
    gen_str=num2str(j);
    figure('OuterPosition',[1 1 scrsz(3)/2 scrsz(4)],'visible',visibility)
    plot(apgrid(j,:),0,'x');
    set(gca,'FontSize',16);
    title(['apgrid for generation ',gen_str,', ',dir]);
    xlabel('ap','FontSize',16);
end

figure('OuterPosition',[1 1 scrsz(3)/2 scrsz(4)],'visible',visibility)
plot(apgrid);
ylabel('ap','FontSize',16);
set(gca,'FontSize',16);
title(['Visualising apgrids for all generations, ',dir]);
print('-depsc', ['graphs/apgrids']);
system(['epstopdf graphs/apgrids.eps']);

if (~isempty(strfind(pwd,'insurance'))) 
    cd('..')
end
cd('../../../src_matlab/')
