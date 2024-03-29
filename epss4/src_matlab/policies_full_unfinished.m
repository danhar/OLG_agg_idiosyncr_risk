%*******************************************************************************
% Copyright (c) 2016 Daniel Harenberg - All rights reserved.
%*******************************************************************************

clear;
%close all;

setenv('EPSSDIR','ge')
dir = getenv('EPSSDIR');
plot_bonds = 0;

nmu=10;
nk=10;
nj=70;
nz=4;
neta=5;
nx=20;

j0= 0; %21;
mean_mu = floor(nmu/2);
mean_k  = floor(nk/2);

cd(['../model_output/',dir])
apgrid_in=dlmread('apgrid.txt')';
kappa_in = dlmread('kappa.txt')';
xgrid_in=dlmread('xgrid.txt')';
agg_grid_k=dlmread('agg_grids.txt','\t',['A2..A',num2str(nk+1)])';
%cons = dlmread('cons_mean.txt');
cd('../../src_matlab/')

xgrid=zeros(nx,neta,nz,nj,nk,nmu);
apgrid=xgrid;
kappa=xgrid;

line=0;
for muc=1:nmu
for kc=1:nk
for jc=1:nj
    for zc=1:nz
        for ec=1:neta
            line = line+1;
            xgrid(:,ec,zc,jc,kc,muc) = xgrid_in(:,line);
            apgrid(:,ec,zc,jc,kc,muc) = apgrid_in(:,line);
            kappa(:,ec,zc,jc,kc,muc) = apgrid_in(:,line);
        end
    end  
end
end
end
cons = xgrid-apgrid;

mut = mean_mu;
kt  = mean_k;
zt  = 1;

scrsz = get(0,'ScreenSize');
for jc=[1,10,20,30,40,50,60]%[39, 40, 41,42,43,44] % %[40,45,50,57] 
    gen_str=num2str(jc+j0);

    xmin = xgrid(1,1,1,jc,:,mean_mu);
    xgrid_new

    
    figure('OuterPosition',[1 1 scrsz(3)/2 scrsz(4)])
        
    subplot(2,2,1)
    surf(xgrid(:,1,1,jc,:,mean_mu),agg_grid_k,cons(:,1,1,jc,:,mean_mu)); hold;
    surf(xgrid(:,2,1,jc,:,mean_mu),agg_grid_k,cons(:,2,1,jc,:,mean_mu)); hold
    
    h=legend(['cons^{',gen_str,'}(x)'],'Location','SouthEast');
    xlabel('x','FontSize',16); ylabel(['cons^{',gen_str,'}(x)'],'FontSize',16);
    set(gca,'FontSize',16);
    title(['cons^{',gen_str,'}(x,z=1) ',dir]);

    subplot(2,2,3)
    surf(xgrid(:,1,4,jc,:,mean_mu),agg_grid_k,cons(:,1,4,jc,:,mean_mu)); hold;
    surf(xgrid(:,2,4,jc,:,mean_mu),agg_grid_k,cons(:,2,4,jc,:,mean_mu)); hold
    
    h=legend(['cons^{',gen_str,'}(x)'],'Location','SouthEast');
    xlabel('x','FontSize',16); ylabel(['cons^{',gen_str,'}(x)'],'FontSize',16);
    set(gca,'FontSize',16);
    title(['cons^{',gen_str,'}(x,z=4) ',dir]);
    
%    subplot(2,2,2)
end
