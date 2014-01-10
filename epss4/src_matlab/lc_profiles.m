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

%setenv('EPSSINPUTFILE','calib_ep/FD1.txt') % 'ge'
%dir_input = getenv('EPSSINPUTFILE');
%fid = fopen(['../model_input/',dir_input]','r');
%param_data = textscan(fid,'%s%n%*[^\n]','HeaderLines',2,'CommentStyle','!');
%fclose(fid);
%j0 = param_data{econ_life_start}{1}

cd(['../model_output/',dir])

all_generations_unit_mass = 1;
ap_lc  = dlmread('ap_lc.txt');
nj = length(ap_lc);

if (all_generations_unit_mass)   
    pop_frac = dlmread('pop_frac.txt');
    pop_frac_rescale=pop_frac*double(nj);
    disp('Plotting lifecycles conditional on surviving to max age!');
else
   pop_frac_rescale = 1.0;
   disp('Plotting expected lifecycle at birth (i.e. including survival rates, i.e. population shrinks)');
end
ap_lc           = dlmread('ap_lc.txt')          ./pop_frac_rescale;
kappa_lc        = dlmread('kappa_lc.txt');
stock_lc        = dlmread('stock_lc.txt')       ./pop_frac_rescale;
cons_lc         = dlmread('cons_lc.txt')        ./pop_frac_rescale;
consvar_lc      = dlmread('consvar_lc.txt')     ./pop_frac_rescale;
return_lc       = dlmread('return_lc.txt')      ./pop_frac_rescale;
returnvar_lc    = dlmread('return_var_lc.txt')  ./pop_frac_rescale;
log_cons_lc     = dlmread('log_cons_lc.txt')    ./pop_frac_rescale;
var_log_cons_lc = dlmread('var_log_cons_lc.txt')./pop_frac_rescale;
bond_lc    = ap_lc - stock_lc;
% cd('../../../src_matlab/')

j0 = 21;

scrsz = get(0,'ScreenSize');
offset = -35;
agegrid=[j0+1:j0+nj];

if (visibility=='off')
    fontsize=12;
    linewidth =1.5;
else
    fontsize=16;
    linewidth =2.5;
end

figure('OuterPosition',[1 1 scrsz(3)/2 scrsz(4)],'visible',visibility)

subplot(3,2,1)
plot(agegrid,kappa_lc,'-o','LineWidth',linewidth);
if (min(kappa_lc) < max(kappa_lc))
    axis([agegrid(1) agegrid(end) min(kappa_lc) max(kappa_lc)])
end
%xlabel('age','FontSize',fontsize); ylabel('kappa','FontSize',fontsize);
set(gca,'FontSize',fontsize);
title('Portfolio choice');

subplot(3,2,2)
plot(agegrid,kappa_lc,'LineWidth',linewidth);
axis([agegrid(1) agegrid(end) -1 2])
%xlabel('age','FontSize',fontsize); ylabel('kappa, zoomed','FontSize',fontsize);
set(gca,'FontSize',fontsize);
%title(['Portfolio choice life cycle profile (zoomed), ',dir]);
title('Portfolio choice (zoomed)');

subplot(3,2,3)
plot(agegrid,bond_lc,'LineWidth',linewidth);
if (min(bond_lc) < max(bond_lc))
    axis([agegrid(1) agegrid(end) min(bond_lc) max(bond_lc)])
end
%xlabel('age','FontSize',fontsize); ylabel('bonds','FontSize',fontsize);
set(gca,'FontSize',fontsize);
%title(['Bonds life cycle profile (=ap-stock), ',dir]);
title('Bonds');

subplot(3,2,4)
plot(agegrid,stock_lc,'LineWidth',linewidth);
if (min(stock_lc) < max(stock_lc))
    axis([agegrid(1) agegrid(end) min(stock_lc) max(stock_lc)])
end
%xlabel('age','FontSize',fontsize); ylabel('stock','FontSize',fontsize);
set(gca,'FontSize',fontsize);
title('Stocks');

subplot(3,2,5)
plot(agegrid,ap_lc,'LineWidth',linewidth);
if (min(ap_lc) < max(ap_lc))
    axis([agegrid(1) agegrid(end) min(ap_lc) max(ap_lc)])
end
xlabel('age','FontSize',fontsize); %ylabel('ap','FontSize',fontsize);
set(gca,'FontSize',fontsize);
title('Total savings');


if (all_generations_unit_mass == 1)
    note='Profiles are conditional on surviving to max age';
else
    note='Expected profiles at birth, incl survival risk';    
end
stringmatrix=char(...
'{\bf            Lifecycle profiles over age}',...
'  ',...
note, ...
['Calibration is:  ',dir]);

text(108,0.15,stringmatrix);

print('-depsc', ['graphs/lifecycles']);
system(['epstopdf graphs/lifecycles.eps']);


figure('OuterPosition',[scrsz(3)/2 1 scrsz(3)/2 scrsz(4)],'visible',visibility)

subplot(2,2,1)
plot(agegrid,cons_lc,'LineWidth',linewidth);
if (min(cons_lc) < max(cons_lc))
    axis([agegrid(1) agegrid(end) min(cons_lc) max(cons_lc)])
end
%xlabel('age','FontSize',fontsize); ylabel('cons','FontSize',fontsize);
set(gca,'FontSize',fontsize);
title('Consumption');

subplot(2,2,2)
plot(agegrid,consvar_lc,'LineWidth',linewidth);
if (min(consvar_lc) < max(consvar_lc))
    axis([agegrid(1) agegrid(end) min(consvar_lc) max(consvar_lc)])
end
xlabel('age','FontSize',fontsize); %ylabel('var(cons)','FontSize',fontsize);
set(gca,'FontSize',fontsize);
title('Cross-sectional variance of consumption');

subplot(2,2,3)
plot(agegrid,log_cons_lc,'LineWidth',linewidth);
if (min(log_cons_lc) < max(log_cons_lc))
    axis([agegrid(1) agegrid(end) min(log_cons_lc) max(log_cons_lc)])
end
%xlabel('age','FontSize',fontsize); ylabel('cons','FontSize',fontsize);
set(gca,'FontSize',fontsize);
title('Log Consumption');

subplot(2,2,4)
plot(agegrid,var_log_cons_lc,'LineWidth',linewidth);
if (min(var_log_cons_lc) < max(var_log_cons_lc))
    axis([agegrid(1) agegrid(end) min(var_log_cons_lc) max(var_log_cons_lc)])
end
xlabel('age','FontSize',fontsize); %ylabel('var(cons)','FontSize',fontsize);
set(gca,'FontSize',fontsize);
title('Cross-sectional variance of log consumption');

print('-depsc', ['graphs/lifecycles2']);
system(['epstopdf graphs/lifecycles2.eps']);

figure('OuterPosition',[scrsz(3)/2 1 scrsz(3)/2 scrsz(4)],'visible',visibility)

subplot(2,2,1)
plot(agegrid,return_lc,'LineWidth',linewidth);
if (min(return_lc) < max(return_lc))
    axis([agegrid(1) agegrid(end) min(return_lc) max(return_lc)])
end
%xlabel('age','FontSize',fontsize); ylabel('avg portf return','FontSize',fontsize);
set(gca,'FontSize',fontsize);
title('Portfolio return');

subplot(2,2,2)
plot(agegrid,returnvar_lc,'LineWidth',linewidth);
if (min(returnvar_lc) < max(returnvar_lc))
    axis([agegrid(1) agegrid(end) min(returnvar_lc) max(returnvar_lc)])
end
xlabel('age','FontSize',fontsize); %ylabel('var(portf return)','FontSize',fontsize);
set(gca,'FontSize',fontsize);
title('Variance of portfolio return');

print('-depsc', ['graphs/lifecycles3']);
system(['epstopdf graphs/lifecycles3.eps']);

if (~isempty(strfind(pwd,'insurance'))) 
    cd('..')
end
cd('../../../src_matlab/')
