% Copyright (C) 2016 Daniel Harenberg - All Rights Reserved
%clear;
close all;

calibname = getenv('CALIBNAME');
if (strcmp(calibname,''))    
    visibility = 'on ';
    calibname = 'FD1_KK';  % 'ge'
else
    visibility = 'off';
end

cd('../model_output/')

data  = dlmread(['welfare_',calibname,'.txt'],'', 1,0);

% cd('../../../src_matlab/')

risk_scale = data(:,1);
welfare_0 = data(:,2);
if (length(welfare_0) <2)
    disp('Only one data point, cannot regress; stopping execution');
    cd('../src_matlab/')
    return
end
if (size(data,2)>2)
    cev_present = true;    
    welfare_1 = data(:,3);
    cev = data(:,4);    
end

%%%%%%%%%%%%%%%%%%% Calculations  %%%%%%%%%%%%%%%%%%%%%
index_1_0 = find(risk_scale == 1.0); % index of baseline, i.e. risk_scale = 1.0
sample_size = 20; %length(welfare_0); %21; %floor(length(welfare_0)/2); % length(welfare_0); %
degree = 3; % largest degree of monomial in interpolation scheme

if (index_1_0 == 1) 
    sample_index_1 = 1;
    sample_index_2 = 1 + sample_size;
else
    sample_index_1 = index_1_0 - floor(sample_size/2);
    sample_index_2 = index_1_0 + floor(sample_size/2);
end

if (sample_index_1 < 1)
    sample_index_1 = 1;
    disp('Warning: sample_size too large, setting index_1 = 1');
end

if (sample_index_2 > length(welfare_0))
    sample_index_2 = length(welfare_0);
    disp('Warning: sample_size too large, setting index_2 = length(welfare_0)');
end

pr_x=linspace(0.0,risk_scale(1),100)';


if (degree == 1)
    Xmat=[ones(length(risk_scale(sample_index_1:sample_index_2)),1),risk_scale(sample_index_1:sample_index_2)];
    pr_Xmat=[ones(length(pr_x),1),pr_x];
elseif (degree == 2)
    Xmat=[ones(length(risk_scale(sample_index_1:sample_index_2)),1),risk_scale(sample_index_1:sample_index_2),risk_scale(sample_index_1:sample_index_2).^2];
    pr_Xmat=[ones(length(pr_x),1),pr_x,pr_x.^2];
elseif (degree == 3)
    Xmat=[ones(length(risk_scale(sample_index_1:sample_index_2)),1),risk_scale(sample_index_1:sample_index_2),risk_scale(sample_index_1:sample_index_2).^2,risk_scale(sample_index_1:sample_index_2).^3];
    pr_Xmat=[ones(length(pr_x),1),pr_x,pr_x.^2,pr_x.^3];
elseif (degree == 4)
    Xmat=[ones(length(risk_scale(sample_index_1:sample_index_2)),1),risk_scale(sample_index_1:sample_index_2),risk_scale(sample_index_1:sample_index_2).^2,risk_scale(sample_index_1:sample_index_2).^3,risk_scale(sample_index_1:sample_index_2).^4];
    pr_Xmat=[ones(length(pr_x),1),pr_x,pr_x.^2,pr_x.^3,pr_x.^4];

else
    disp('Regression of degree =',num2str(degree),' not implemented!')
    error('Stopping execution')
end    

betta_0=Xmat\welfare_0(sample_index_1:sample_index_2);
pr_welfare_0=pr_Xmat*betta_0;

if (cev_present)
    betta_1=Xmat\welfare_1(sample_index_1:sample_index_2);
    pr_welfare_1=pr_Xmat*betta_1;
    cev0=betta_1(1)/betta_0(1)-1;
    
    betta_cev=Xmat\cev(sample_index_1:sample_index_2);
    pr_cev = pr_Xmat*betta_cev;    
end

disp(' ')
disp(['predicted value for risk_scale=0,tau=0 ', num2str(betta_0(1))]);

if (cev_present) % this if goes to end of file

disp(['predicted value for risk_scale=0,tau=0.02 ', num2str(betta_1(1))]);
disp(['cev(risk_scale=0) from welfare predictions ', num2str(cev0)]);
disp(['cev(risk_scale=0) from cev prediction ', num2str(betta_cev(1))]);
disp(' ')

copyfile(['welfare_',calibname,'.txt'],['welfare_',calibname,'_pr.txt']); % only need this for testing 
file1 = fopen(['welfare_',calibname,'_pr.txt'],'a');
fprintf(file1,'\n\n');
fprintf(file1,'Predicted values from welfare regressions of degree=%d:\n',degree);
fprintf(file1,' 0.00* %15.6E  %15.6E  %15.6E\n', betta_0(1), betta_0(1), cev0);
fprintf(file1,'and from cev regression:                 %15.6E', betta_cev(1));
fclose(file1);
%%%%%%%%%%%%%%%%%%% Plotting  %%%%%%%%%%%%%%%%%%%%%

scrsz = get(0,'ScreenSize');
offset = -35;

if (visibility=='off')
    fontsize=10;
    linewidth =1.0;
    markersize=6;
else
    fontsize=16;
    linewidth =2.0;
    markersize=8;
end

figure('OuterPosition',[1 1 scrsz(3)/2 scrsz(4)],'visible',visibility)

subplot(2,2,1)
plot(risk_scale,welfare_0,'bx',risk_scale,welfare_1,'rx','LineWidth',linewidth,'MarkerSize',markersize);
%xlabel('age','FontSize',fontsize); ylabel('kappa','FontSize',fontsize);
set(gca,'FontSize',fontsize);
%xlim([min(risk_scale), max(risk_scale)])
title('Welfare indeces');
legend('tau = 0.00', 'tau = 0.02','Location','Best')
yL = get(gca,'YLim');
line([risk_scale(sample_index_1) risk_scale(sample_index_1)], yL)
line([risk_scale(sample_index_2) risk_scale(sample_index_2)], yL)

subplot(2,2,2)
plot(pr_x,pr_welfare_0,'b-',risk_scale,welfare_0,'bx',pr_x,pr_welfare_1,'r--',risk_scale,welfare_1,'rx','LineWidth',linewidth,'MarkerSize',markersize);
%xlabel('age','FontSize',fontsize); ylabel('kappa, zoomed','FontSize',fontsize);
set(gca,'FontSize',fontsize);
%title(['Portfolio choice life cycle profile (zoomed), ',dir]);
title('Welfare indeces, predicted');
legend('tau = 0.00','tau = 0.00', 'tau = 0.02','tau = 0.02','Location','Best')
yL = get(gca,'YLim');
line([risk_scale(sample_index_1) risk_scale(sample_index_1)], yL)
line([risk_scale(sample_index_2) risk_scale(sample_index_2)], yL)


subplot(2,2,3)
plot(pr_x,(pr_welfare_1./pr_welfare_0)-1,'b-',pr_x,pr_cev,'r--',risk_scale,cev,'rx','LineWidth',linewidth,'MarkerSize',markersize);
%xlabel('age','FontSize',fontsize); ylabel('kappa, zoomed','FontSize',fontsize);
set(gca,'FontSize',fontsize);
%title(['Portfolio choice life cycle profile (zoomed), ',dir]);
title('CEV, predicted');
legend('from welfare', 'from cev', 'datapoints','Location','Best')
yL = get(gca,'YLim');
line([risk_scale(sample_index_1) risk_scale(sample_index_1)], yL)
line([risk_scale(sample_index_2) risk_scale(sample_index_2)], yL)

print('-depsc', ['welfare_',calibname]);
system(['epstopdf welfare_',calibname,'.eps']);

end % ends if (cev_present)

cd('../src_matlab/')
