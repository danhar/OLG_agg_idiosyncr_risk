% Copyright (C) 2016 Daniel Harenberg - All Rights Reserved

global gplot_all
%clear;
close all;

if (isempty(gplot_all))    
    visibility = 'on ';
    dir = 'FD1_KK/ge';  % 'ge'
else
    visibility = 'off';
    dir = getenv('EPSSDIR');
end

cd(['../model_output/',dir])
Phi = dlmread('Phi_tilde.txt');
xgrid_mean = dlmread('xgrid_mean.txt');

grid_error = (xgrid_mean(:,1)==xgrid_mean(:,2));
% grid_error can happen because a) numerical Euler equation can become
% inaccurate at the very small values, and b) the Fortran results are saved
% with small precision (6 digits). b) is much more likely
if (any(grid_error))
    xgrid_mean(grid_error,2) = (xgrid_mean(grid_error,1) + xgrid_mean(grid_error,3))/2;
end

[nj nx] = size(Phi);
nx_interp = min([nx*4, 100]); % min so that it doesn't become too large, which is a problem when printing graph on cluster
Phi_cutoff=0.00001;

scrsz = get(0,'ScreenSize');
if (visibility=='off')
    fontsize=12;
    linewidth =1.5;
    markersize = 7;
else
    fontsize=16;
    linewidth =2.5;
    markersize = 9;
end

nx2 = size(xgrid_mean,2);
if (nx2 ~= nx) % usually they differ, bc nx_factor makes finer grids for simulations, and because usually a fixed grid is created for the simulations
    xgrid_new = zeros(nj,nx);
    c =2; % 2 is the value used in interpolate_xgrid.f90 (which creates the xgrid for the simulations)
    for j=1:size(xgrid_mean,1)
        xgrid_new(j,1) = xgrid_mean(j,1);
        xgrid_new(j,nx) = xgrid_mean(j,nx2);
  		scalefact = xgrid_new(j,nx)-xgrid_new(j,1);
        for i=2:nx-1
		    xgrid_new(j,i) = xgrid_new(j,1) + scalefact*((i-1.0)/(nx-1.0))^c;
        end
    end
    % This should be very close to the 'true' (average) grid over the simulations
    xgrid_mean = xgrid_new;
end

figure('OuterPosition',[1 1 scrsz(3)/2 scrsz(4)],'visible',visibility)
count=0;
for i=[1,2,3,4] %[1,15,30,45] %[1, 20, 32, 40, 60]
    count=count+1;
    subplot(2,2,count)
    %bar(xgrid_mean(i,:),Phi(i,:));
    plot(xgrid_mean(i,:),Phi(i,:),'-x','MarkerSize',markersize,'LineWidth',linewidth);
    title(['Phi_{',dir,'}^{',num2str(i),'}(x)']);
end

print('-depsc', ['graphs/Phi_j']);
system(['epstopdf graphs/Phi_j.eps']);

% Need to create a grid and interpolate, bc
% (a) first and foremost, want to throw out very small values to get a nice graph
% (b) as mentioned above, we do not have the grid that exactly corresponds to Phi, bc Phi is average over simulations
%     However, xgrid_mean should be a very good approximation, so (b) is not the reason for interpolating here
Phi_interp=zeros(nj-1,nx_interp); % nj-1, because all dlmread variables are read in with one line too much...
xmax=max(max(xgrid_mean(Phi>Phi_cutoff)));
xmin=min(min(xgrid_mean(Phi>Phi_cutoff)));
xgrid=(xmin:((xmax-xmin)/(nx_interp-1)):xmax)';
for j=1:size(Phi_interp,1)
    Phi_interp(j,:)=interp1(xgrid_mean(j,:),Phi(j,:), xgrid,'linear',0);
end
figure('OuterPosition',[scrsz(3)/2 1 scrsz(3)/2 scrsz(4)],'visible',visibility)
surf((1:size(Phi_interp,1)),xgrid,Phi_interp');
set(gca,'FontSize',fontsize);
title(['\Phi_{',dir,'}(x,j), for \Phi>',num2str(Phi_cutoff)]);
view([-25 60]);

print('-depsc', ['graphs/Phi_tilde']);
system(['epstopdf graphs/Phi_tilde.eps']);

if (~isempty(strfind(pwd,'insurance'))) 
    cd('..')
end
cd('../../../src_matlab/')
