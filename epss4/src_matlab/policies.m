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

plot_bonds = 0;

cd(['../model_output/',dir])
apgrid=dlmread('apgrid_mean.txt');
kappa = dlmread('kappa_mean.txt');
stocks=dlmread('stocks_mean.txt');
xgridM=dlmread('xgrid_mean.txt');
%cons = dlmread('cons_mean.txt');
cons = xgridM - apgrid;

[nj,nx] = size(apgrid);
j0= 0; %21;

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

jplot=[1,10:10:nj];
if (jplot(end) ~= nj) 
    jplot= [jplot,nj];
end
for j=jplot%[39, 40, 41,42,43,44] % %[40,45,50,57] 
    gen_str=num2str(j+j0);
    xgrid=xgridM(j,:);

    figure('OuterPosition',[1 1 scrsz(3)/2 scrsz(4)],'visible',visibility)
        
    subplot(2,2,1)
    plot(xgrid,cons(j,:),'LineWidth',linewidth);
    h=legend(['cons^{',gen_str,'}(x)'],'Location','SouthEast');
    %xlabel('x','FontSize',fontsize); ylabel(['cons^{',gen_str,'}(x)'],'FontSize',fontsize);
    set(gca,'FontSize',fontsize);
    axis([xgrid(1) xgrid(end) min(cons(j,:)) max(cons(j,:))])
    title(['\bf Calibration ',dir]);

    subplot(2,2,2)
    plot(xgrid,kappa(j,:),'bx','MarkerSize',markersize,'LineWidth',linewidth); % 
    h=legend(['\kappa^{',gen_str,'}(x)'],'Location','SouthEast');
    %xlabel('x','FontSize',fontsize); ylabel(['\kappa^{',gen_str,'}(x)'],'FontSize',fontsize);
    set(gca,'FontSize',fontsize);
    if (min(kappa(j,:)) < max(kappa(j,:)))
        axis([xgrid(1) xgrid(end) min(kappa(j,:)) max(kappa(j,:))])
    else
        axis([xgrid(1) xgrid(end) -1 1])
    end
        
    %title(['\kappa^{',gen_str,'}(x), ',dir]);
   
    subplot(2,2,3)
    plot(xgrid,apgrid(j,:),'LineWidth',linewidth);
    h=legend(['ap^{',gen_str,'}(x)'],'Location','SouthEast');
    %xlabel('x','FontSize',fontsize); ylabel(['ap^{',gen_str,'}(x)'],'FontSize',fontsize);
    set(gca,'FontSize',fontsize);
    if (min(apgrid(j,:)) < max(apgrid(j,:)))
        axis([xgrid(1) xgrid(end) min(apgrid(j,:)) max(apgrid(j,:))])
    else
        axis([xgrid(1) xgrid(end) -1 1])
    end    
    %title(['ap^{',gen_str,'}(x), ',dir]);
    
    subplot(2,2,4)
    plot(xgrid,stocks(j,:),'LineWidth',linewidth);
    h=legend(['stocks^{',gen_str,'}(x)'],'Location','SouthEast');
    %xlabel('x','FontSize',fontsize); ylabel(['stocks^{',gen_str,'}(x)'],'FontSize',fontsize);
    set(gca,'FontSize',fontsize);    
    if (min(stocks(j,:)) < max(stocks(j,:)))
        axis([xgrid(1) xgrid(end) min(stocks(j,:)) max(stocks(j,:))])        
    else
        axis([xgrid(1) xgrid(end) -1 1])
    end
    %title(['stocks^{',gen_str,'}(x), ',dir]);
    
    if (plot_bonds)
        bonds=apgrid(j,:)-stocks(j,:);
        figure('OuterPosition',[scrsz(3)/2 1 scrsz(3)/4 scrsz(4)/2])
        plot(xgrid,bonds,'LineWidth',linewidth);
        h=legend(['bonds^{',gen_str,'}(x)'],'Location','SouthEast');
     %   xlabel('x','FontSize',fontsize); ylabel(['stocks^{',gen_str,'}(x)'],'FontSize',fontsize);
        set(gca,'FontSize',fontsize);
        if (min(bonds) < max(bonds))
            axis([xgrid(1) xgrid(end) min(bonds) max(bonds)])            
        else
            axis([xgrid(1) xgrid(end) -1 1])
        end
%        title(['bonds^{',gen_str,'}(x), ',dir]);
    end
    
    print('-depsc', ['graphs/policies_',gen_str]);
    system(['epstopdf graphs/policies_',gen_str,'.eps']);
end

cd('../../../src_matlab/')
