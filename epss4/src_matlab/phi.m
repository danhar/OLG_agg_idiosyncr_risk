% Copyright (C) 2016 Daniel Harenberg - All Rights Reserved

global gplot_all
%clear;
close all;

if (isempty(gplot_all))    
    visibility = 'on ';
    dir = 'FD1/ge';  % 'ge'
else
    visibility = 'off';
    dir = getenv('EPSSDIR');
end

disp('attention: not well implemented for GE!');


nj=70;%55;%70; %63;
nx=20;
neta=2;
nz=4;
nx_interp=40;
Phi_cutoff=0.0001/2.0;

cd(['../model_output/',dir])
Phi_in = dlmread('Phi.txt')';
xgrid_in= dlmread('xgrid.txt')';
xgrid_mean = dlmread('xgrid_mean.txt')';
Phi_tilde = dlmread('Phi_tilde.txt')';

line1=0;
line2=0;
Phi = zeros(nx,neta,nj);
xgrid_full = zeros(nx,neta,nz,nj);
for jc=1:nj
    for ec=1:neta
        line1=line1+1;
        Phi(:,ec,jc) = Phi_in(:,line1);        
    end
    for zc=1:nz
        for ec=1:neta
            line2 = line2+1;
            xgrid_full(:,ec,zc,jc) = xgrid_in(:,line2);
        end
    end  
end

xgrid=zeros(nx,neta,nj);
for zc=1:nz
    xgrid= xgrid + (1.0/nz)*squeeze(xgrid_full(:,:,zc,:));
end

Phi_tilde2 = reshape(sum(Phi,2),nx,nj);
diff = Phi_tilde - Phi_tilde2;
max_err_Phi=max(max(diff))

xgrid_mean2 = squeeze(1.0/neta*sum(xgrid,2));
diff = xgrid_mean - xgrid_mean2;
max_err_xgrid=max(max(diff))

scrsz = get(0,'ScreenSize');
figure('OuterPosition',[1 1 scrsz(3)/2 scrsz(4)])
count=0;
for i=[1,2,3,4] %[1,15,30,45] %[1, 20, 32, 40, 60]
    count=count+1;
    subplot(2,2,count)
    %bar(xgrid_mean(i,:),Phi(i,:));
    plot(xgrid(:,1,i),Phi(:,1,i),'-kx',xgrid(:,neta,i),Phi(:,neta,i),'-rx','MarkerSize',9,'LineWidth',2);
    title(['Phi_{',dir,'}^{',num2str(i),'}(x,eta)']);
    legend('eta=1', ['eta=',num2str(neta)]);
end

Phi_interp=zeros(nx_interp,nj);
xmax=max(max(max(xgrid(Phi>Phi_cutoff))));
xmin=min(min(min(xgrid(Phi>Phi_cutoff))));
xgrid_new=(xmin:((xmax-xmin)/(nx_interp-1)):xmax)';
for etac=1:neta
    for jc=1:nj
        Phi_interp(:,jc)=interp1(xgrid(:,etac,jc),Phi(:,etac,jc), xgrid_new,'linear',0);
    end
    figure('OuterPosition',[scrsz(3)/2 1 scrsz(3)/2 scrsz(4)])
    surf((1:nj),xgrid_new,Phi_interp);
    set(gca,'FontSize',16);
    title(['\Phi_{',dir,'}(x,',num2str(etac),',j), for \Phi>',num2str(Phi_cutoff)]);
end

for jc=1:nj
    Phi_interp(:,jc)=interp1(xgrid_mean2(:,jc),Phi_tilde2(:,jc), xgrid_new,'linear',0);
end
figure('OuterPosition',[scrsz(3)/2 1 scrsz(3)/2 scrsz(4)])
surf((1:nj),xgrid_new,Phi_interp);
set(gca,'FontSize',16);
title(['\Phi_{',dir,'}(x,mean,j), for \Phi>',num2str(Phi_cutoff)]);
