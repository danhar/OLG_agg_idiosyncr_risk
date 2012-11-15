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
jr =43;
cd(dir)
apgrid=dlmread('apgrid_mean.txt');

scrsz = get(0,'ScreenSize');

for j=[1,35]
    gen_str=num2str(j);
    i= floor((j+1)/jr)+1;
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

cd('..')
cd('../../src_matlab/')
