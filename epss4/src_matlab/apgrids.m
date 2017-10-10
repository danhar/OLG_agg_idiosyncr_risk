%*******************************************************************************
% Copyright (c) 2016 Daniel Harenberg - All rights reserved.
%*******************************************************************************
global gplot_all
%clear;
close all;


if (isempty(gplot_all))    
    visibility = 'on ';
    dir = 'FD1_KK/msge';  % 'ge'
    scrsz = get(0,'ScreenSize');
    fig_pos = 'OuterPosition';
else
    visibility = 'off';
    dir = getenv('EPSSDIR');
    scrsz = [1 1 1920 1080];
    fig_pos = 'PaperPosition';   
end

cd('../model_output/')
% jr = dlmread('params.txt','', [40 11 40 13]);

cd(dir)
apgrid=dlmread('apgrid_mean.txt');

for j=[1,35]
    gen_str=num2str(j);
    figure(fig_pos,[1 1 scrsz(3)/2 scrsz(4)],'visible',visibility,'PaperPositionMode','auto')
    plot(apgrid(j,:),0,'x');
    set(gca,'FontSize',16);
    title(['apgrid for generation ',gen_str,', ',dir]);
    xlabel('ap','FontSize',16);
end

figure(fig_pos,[1 1 scrsz(3)/2 scrsz(4)],'visible',visibility,'PaperPositionMode','auto')
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
