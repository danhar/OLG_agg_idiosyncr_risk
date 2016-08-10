% Copyright (C) 2016 Daniel Harenberg - All Rights Reserved

clear;
%close all;

setenv('EPSSDIR','pe')
dir = getenv('EPSSDIR');

cd(['../model_output/',dir])
simvars  = dlmread('simvars.txt','',1, 0);
z=simvars(1,:);
K=simvars(2,:);
mu=simvars(3,:);
B=simvars(4,:);
C=simvars(5,:);
r=simvars(6,:);
rf=simvars(7,:);
wage=simvars(8,:);
pens=simvars(9,:);
cd('../../src_matlab')

scrsz = get(0,'ScreenSize');
figure('OuterPosition',[1 1 scrsz(3)/2 scrsz(4)])

subplot(2,2,1)
plot(K); %,'LineWidth',2.5
%axis([agegrid(1) agegrid(end) min(kappa_lc) max(kappa_lc)])
xlabel('time','FontSize',16); ylabel('K','FontSize',16);
set(gca,'FontSize',16);
title(['Simulated aggregate capital, ',dir]);

subplot(2,2,3)
plot(mu);
%axis([agegrid(1) agegrid(end) -1 2])
xlabel('time','FontSize',16); ylabel('mu','FontSize',16);
set(gca,'FontSize',16);
title(['Simulated (realized) equity premium, ',dir]);

subplot(2,2,2)
plot(C);
%axis([agegrid(1) agegrid(end) min(stock_lc) max(stock_lc)])
xlabel('time','FontSize',16); ylabel('C','FontSize',16);
set(gca,'FontSize',16);
title(['Simulated aggregate consumption, ',dir]);

subplot(2,2,4)
plot(z,'x');
%axis([agegrid(1) agegrid(end) min(bond_lc) max(bond_lc)])
xlabel('time','FontSize',16); ylabel('z','FontSize',16);
set(gca,'FontSize',16);
title(['Simulated shock process, ',dir]);


figure('OuterPosition',[scrsz(3)/2 1 scrsz(3)/4 scrsz(4)])

subplot(2,1,1)
plot(r);
%axis([agegrid(1) agegrid(end) min(cons_lc) max(cons_lc)])
xlabel('time','FontSize',16); ylabel('r','FontSize',16);
set(gca,'FontSize',16);
title(['Simulated realized risky return, ',dir]);

subplot(2,1,2)
plot(rf);
%axis([agegrid(1) agegrid(end) min(ap_lc) max(ap_lc)])
xlabel('time','FontSize',16); ylabel('rf','FontSize',16);
set(gca,'FontSize',16);
title(['Simulated realized riskfree return, ',dir]);