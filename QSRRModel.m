%% The data and codes for 
% Title: Application of Bayesian Multilevel Modeling in Quantitative Structure-Retention Relationship Studies of Heterogeneous Compounds 
% Authors: Pawe³ Wiczling, Agnieszka Kamedulska, £ukasz Kubik
% Adress: Department of Biopharmaceutics and Pharmacodynamics, Medical University of Gdañsk, Gen. J. Hallera 107, 80-416 Gdañsk, Poland
% Data: 30/10/2020
% Version 1.0
%% Load data
clc;
clear all;
data = readtable('Data\database_logk_1026.csv');
analyte_names = readtable('Data\database_logk_1026_analyte_names.csv');
functional_groups = readtable('Data\checkmol_functional_groups.csv');
functional_groups_names = readtable('Data\checkmol_functional_group_names.csv');

% combine nr of caroboxylic acid and carboxyalic acid salt functional groups
functional_groups{:,76}=functional_groups{:,76}+functional_groups{:,77};       
functional_groups{functional_groups{:,202}>8.1,202} = 8; % heterocyclic compounds with more than 8 heterocycles are treated as if they have 8 (strychnine)

% exclude functional groups that repeat itself (some groups are nested)
idx_excluded = [1 2 3 6 27 28 37 47 48 51 55 61 62 67 73 74 75 77 80 91 99 109 116 117 121 125 129 142 153 154 160 161 168 173 178 181 182 186 187 191 196];
writetable(functional_groups_names(idx_excluded,:),'Tables/functional_groups_excluded.csv','Delimiter',',','QuoteStrings',false)
functional_groups_names(idx_excluded,:) = []; functional_groups(:,idx_excluded) = []; clear idx_excluded

% exclude functional groups not present on any analyte from the dataset
idx_not_present = find(sum(functional_groups{:,:})'==0);
writetable(functional_groups_names(idx_not_present,:),'Tables/functional_groups_not_present.csv','Delimiter',',','QuoteStrings',false)
functional_groups_names(idx_not_present,:) = []; functional_groups(:,idx_not_present) = []; clear idx_not_present

%% Raw data
%Figure S1. Relationship between the logarithm of retention factor (log k)
%and acetonitrile content in the mobile phase. Lines connect measurements
%obtained for a particular analyte.

figure('Color',[1 1 1])
h1 = gscatter(data.concentration,data.logk,data.ID);
set(h1,'linestyle', '-')
xlabel('$$\varphi$$ (ACN)','Interpreter','latex')
ylabel('logk')
legend off
box off
clear h1 

savefig('Figures/FigureS1_RawData.fig')
set(gcf,'paperunits','centimeters','paperposition',[0 0 16.5 18])
print -dtiff -r300 Figures/FigureS1_RawData.tif
%% Functional groups
% Figure 1. Functional groups identified by Checkmol. Figures show the
% number of analytes having at least one functional group of a given type.

[SortedSum,I] = sort(sum(functional_groups{:,:}>0.5));
figure('Color',[1 1 1])
subplot(1,2,1)
plot(1:1:50,SortedSum([1:1:50]),'-o')
xlabel('Functional group')
ylabel('                                                                      Number of analytes having at least one functional group of a given type')
view(90,90)
set(gca,'Xtick',[1:1:50],'XTickLabelRotation',0,'XTickLabel',functional_groups_names{I([1:1:50]),2})
set(gca,'Yscale','lin','FontSize',8)
subplot(1,2,2)
plot(51:1:100,SortedSum([51:1:100]),'-o')
view(90,90)
set(gca,'Xtick',[51:1:100],'XTickLabelRotation',0,'XTickLabel',functional_groups_names{I([51:1:100]),2})
set(gca,'Yscale','log','FontSize',8)
clear I SortedSum 

savefig('Figures/Figure1_FunctionalGroups.fig')
set(gcf,'paperunits','centimeters','paperposition',[0 0 16.5 18])
print -dtiff -r300 Figures/Figure1_FunctionalGroups.tif
%% Aproximate individual parameters by least-sqaure method (two-stage aproach)
nAnalytes = length(unique(data.ID));
initial_paramA = zeros(nAnalytes,2);
logkPred = 0.*data.logk;

for i = 1:nAnalytes
    if length(data.concentration(data.ID==i))==1
    initial_paramA(i, 2) = polyfit(-16*data.concentration(data.ID==i) ./ (1 + 2 .* data.concentration(data.ID==i)), data.logk(data.ID==i), 0);
    initial_paramA(i, 1) = -16;
    else
	initial_paramA(i, :) = polyfit(data.concentration(data.ID==i) ./ (1 + 2 .* data.concentration(data.ID==i)), data.logk(data.ID==i), 1);
    end
    logkPred(data.ID==i) = initial_paramA(i, 2) + initial_paramA(i, 1) .* data.concentration(data.ID==i) ./ (1 + 2 .* data.concentration(data.ID==i));  
end
initial_paramA(:,3) = 2*ones(nAnalytes,1);

res_std = sqrt(mean((data.logk - logkPred).^2));

[~,i1,~]=unique(data.ID,'first');

ID = data.ID(i1);
MW_ACD = (data.MW_ACD(i1)-300)/100;
logkw = initial_paramA(:,2);
logka = initial_paramA(:,2) + initial_paramA(:,1)./3;
logS2 = log10(initial_paramA(:,3));

Parameters_basic_fit = table(ID,MW_ACD,logkw,logka,logS2);

% Visualise classical fits for 10 randomly selected analytes
figure('Color', [1 1 1]);
rng(3333)
k = datasample(1:nAnalytes,10,'Replace',false);
fi= 0:0.01:1;
for i = 1:numel(k)
    subplot(5, 2, i)
	hold on
	plot(data.concentration(data.ID==k(i)), data.logk(data.ID==k(i)),  'k.')
	plot(fi, initial_paramA(k(i), 2) + initial_paramA(k(i), 1) .* fi ./ (1 + 2 .* fi), 'k-')
	xlim([0 1])
	ylim([-2.5 4])
    title(analyte_names.Analyte(k(i)),'FontSize',8)
    if i==5
        ylabel('logk')
    end
     if (i==9) || (i==10)
       xlabel('$$\varphi$$ (ACN)','Interpreter','latex')
     end   
end

clear k clear logS2 logka logkw MW_ACD ID i i1 z1 fi logkPred
%% Determine the center and scale. It is used to construct priors.
figure('Color', [1 1 1]);
subplot(2,1,1)
hold on
plot(Parameters_basic_fit.MW_ACD, Parameters_basic_fit.logkw,'.')
z1 = fitglm(Parameters_basic_fit.MW_ACD(~isnan(Parameters_basic_fit.MW_ACD)), Parameters_basic_fit.logkw(~isnan(Parameters_basic_fit.MW_ACD)));
plot([-3:0.1:4]',predict(z1,[-3:0.1:4]'))
Priors.p_neutral_logkw = z1.Coefficients{:,1}';
Priors.p_std_logkw = std(z1.Residuals{:,1});
xlabel('(Molecular Weight - 300)/100')
ylabel('logk_w')
subplot(2,1,2)
hold on
plot(Parameters_basic_fit.MW_ACD, Parameters_basic_fit.logka,'.')
z1 = fitglm(Parameters_basic_fit.MW_ACD(~isnan(Parameters_basic_fit.MW_ACD)), Parameters_basic_fit.logka(~isnan(Parameters_basic_fit.MW_ACD)));
plot([-3:0.1:4]',predict(z1,[-3:0.1:4]'))
Priors.p_neutral_logka = z1.Coefficients{:,1}';
Priors.p_std_logka = std(z1.Residuals{:,1});
xlabel('(Molecular Weight - 300)/100')
ylabel('logk_a')
clear z1
Priors.resstd = res_std;
Priors
%% Initialize variables and parameters
clc
nObs = length(data.ID);
nAnalytes = length(unique(data.ID));
[~,i1,j]=unique(data.ID,'first');

datastruct = struct(...
    'nObs',nObs, ...
	'nAnalytes', nAnalytes, ...
    'K', size(functional_groups,2),...
	'start',i1', ...
	'analyte',j',...
	'Mmolx',(data.MW_ACD'-300)/100,...
    'nrfungroups',functional_groups{:,:},...
	'fi',data.concentration',...
    'run_estimation', 0, ...
	'logkObs', data.logk);

clear init0
% Initialize the values for each variable in each chain
param = [initial_paramA(:,2) initial_paramA(:,2) + initial_paramA(:,1)./3 log10(initial_paramA(:,3))];
for i=1:4
    S.logkwHat =  normrnd(6.6,1,1);
	S.logkaHat  = normrnd(1.3,1,1) ;
	S.logS2Hat = normrnd(log(2),0.1,1);
    S.beta = [1.4 0.2 0] .* exp(normrnd(0, 1, 1, 3));
    S.rho = diag([1 1 1]);
    S.sigma  = lognrnd(log(0.05),0.2,1,1);
    S.nu    = gamrnd(2,1./0.1);
    S.nuobs  = gamrnd(2,1./0.1);
    S.nupi   = gamrnd(2,1./0.1);
    S.pilogkw = lognrnd(log(0.2),0.5, 1, size(functional_groups,2));
    S.pidlogk = normrnd(0.2, 0.5, 1, size(functional_groups,2));
    S.pilogS2 = normrnd(0, 0.1, 1, size(functional_groups,2));
    S.omega = [1, 1, 0.1] .* exp(normrnd(0, 0.2, 1, 3));
    S.spilogkw= abs(normrnd(0,0.5));
    S.spidlogk= abs(normrnd(0,0.5));
    S.spilogS2= abs(normrnd(0,0.5));
    S.mpilogkw = abs(normrnd(0,0.5));
    S.mpidlogk = normrnd(0,0.5);
    S.param  = param;
    init0(i) = S;
end
clear S i i1 j kaHat kwHat nAnalytes nObs fi param
%% Use Stan. 
% For prior predictive checks
clc
fprintf( 'Running Stan...\n' );
datastruct.run_estimation=0;
fit101p = stan('file','ACNQSRR101.stan','data', datastruct, ...
              'working_dir','Tmpstan','iter',1000,'warmup',1000,'chains',4,'init',init0, ...
              'stan_home', 'D:\cmdstan-2.18.1\cmdstan-2.18.1');
clc

% Include likelihood
fprintf( 'Running Stan...\n' );
datastruct.run_estimation=1;
fit101= stan('file','ACNQSRR101.stan','data', datastruct, ...
              'working_dir','Tmpstan','iter',1000,'warmup',1000,'chains',4,'init',init0, ...
              'stan_home', 'D:\cmdstan-2.18.1\cmdstan-2.18.1');

%% 10-fold cross-validation. Include single observation per analyte with logk about 1
rng(3333)
k = datasample(1:length(unique(data.ID)),100,'Replace',false);
idxcv = 1:1:datastruct.nObs;

for i =1:length(k)
idxcv(datastruct.analyte==k(i))=NaN;    
[~,I]=min(abs(datastruct.logkObs(datastruct.analyte==k(i))-1));
idxcv(datastruct.start(k(i)) + I-1) = datastruct.start(k(i)) + I-1;
end

idxcv(isnan(idxcv))=[];

datastructcv1 = datastruct;

datastructcv1.nEst= length(idxcv);
datastructcv1.cvidx= idxcv;
datastructcv1.logkObsEst=datastruct.logkObs(idxcv)';
datastructcv1.logkObs=datastructcv1.logkObs';

fprintf( 'Running Stan...\n' );
datastruct.run_estimation=1;
fit101cv1= stan('file','ACNQSRR101CV.stan','data', datastructcv1, ...
              'working_dir','Tmpstan','iter',1000,'warmup',1000,'chains',4,'init',init0, ...
              'stan_home', 'D:\cmdstan-2.18.1\cmdstan-2.18.1');
%% 10-fold crossvalidation (1-fold shown for simplicity). Include signle observation per analyte with min value of logk
rng(3333)
k = datasample(1:length(unique(data.ID)),100,'Replace',false);

idxcv = 1:1:datastruct.nObs;

for i =1:length(k)
idxcv(datastruct.analyte==k(i))=NaN;    
[~,I]=min(abs(datastruct.logkObs(datastruct.analyte==k(i))-(-3)));
idxcv(datastruct.start(k(i)) + I-1) = datastruct.start(k(i)) + I-1;
end

idxcv(isnan(idxcv))=[];

datastructcv2 = datastruct;
datastructcv2.nEst= length(idxcv);
datastructcv2.cvidx= idxcv;
datastructcv2.logkObsEst=datastruct.logkObs(idxcv)';
datastructcv2.logkObs=datastructcv2.logkObs';

fprintf( 'Running Stan...\n' );
datastruct.run_estimation=1;
fit101cv2= stan('file','ACNQSRR101CV.stan','data', datastructcv2, ...
              'working_dir','Tmpstan','iter',1000,'warmup',1000,'chains',4,'init',init0, ...
              'stan_home', 'D:\cmdstan-2.18.1\cmdstan-2.18.1');
%% 10-fold crossvalidation (1-fold shown for simplicity). Include signle observation per analyte with min value of logk
rng(3333)
k = datasample(1:length(unique(data.ID)),100,'Replace',false);

idxcv = 1:1:datastruct.nObs;
idxcvmax =[];
idxcvmin =[];
for i =1:length(k)
idxcv(datastruct.analyte==k(i))=NaN; 
I = find(datastruct.fi(datastruct.analyte==k(i))==0.3);
if ~isempty(I)
idxcv(datastruct.start(k(i)) + I-1) = datastruct.start(k(i)) + I-1;
end

if isempty(I) && max(datastruct.logkObs(datastruct.analyte==k(i)))>0.75
idxcvmax = [idxcvmax k(i)];
end
if isempty(I) && max(datastruct.logkObs(datastruct.analyte==k(i)))<=0.75
idxcvmin = [idxcvmin k(i)];
end
end

idxcv(isnan(idxcv))=[];

datastructcv3 = datastruct;
datastructcv3.nEst= length(idxcv);
datastructcv3.cvidx=idxcv;
datastructcv3.nEstmin= length(idxcvmin);
datastructcv3.cvidxmin=idxcvmin;
datastructcv3.nEstmax= length(idxcvmax);
datastructcv3.cvidxmax=idxcvmax;
datastructcv3.logkObsEst=datastruct.logkObs(idxcv)';
datastructcv3.logkObs=datastructcv3.logkObs';
datastructcv3.fix = 0.3;

fprintf( 'Running Stan...\n' );
datastruct.run_estimation=1;
fit101cv3= stan('file','ACNQSRR101CVULLIM.stan','data', datastructcv3, ...
              'working_dir','Tmpstan','iter',1000,'warmup',1000,'chains',4,'init',init0, ...
              'stan_home', 'D:\cmdstan-2.18.1\cmdstan-2.18.1');
fit101cv3.block();
%% Summary of model parameters. Save to file
% Table 1. Summary of the MCMC Simulations of the Marginal Posterior
% Distributions of Population-Level Model Parameters. Mean Denotes Sample
% Mean, MCSE Denotes Monte Carlo Standard Error, StdDev Denotes Sample
% Standard Deviation, 5%, 50%, 95% Denote Corresponding Quantiles, N_Eff
% Denotes Effective Sample Size, R_Hat Denotes a Measure of Chain
% Equilibrium, must be within 0.05 of 1.0.

diary ACNFits101.txt
fit101.print();
diary off
save ACNFits101 '-v7.3'
%% load saved data
load ACNFits101
%% get samples
samplesp = fit101p.extract;
samples  = fit101.extract;
samplescv1 = fit101cv1.extract;
samplescv2 = fit101cv2.extract;
samplescv3 = fit101cv3.extract;
%% Posterior summary (used later as priors)
pmpilogkw = mean(samples.pilogkw)';
pspilogkw = std(samples.pilogkw)';
pmpidlogk = mean(samples.pidlogk)';
pspidlogk = std(samples.pidlogk)';
pmpilogS2 = mean(samples.pilogS2)';
pspilogS2 = std(samples.pilogS2)';

for i=1:10
    figure
    subplot(3,1,1)
    hold on
    histogram(samples.pilogkw(:,i), 'Normalization','pdf'); 
    plot(0:0.01:3,normpdf(-0:0.01:3,pmpilogkw(i),pspilogkw(i)))
    subplot(3,1,2)
    hold on
    histogram(samples.pidlogk(:,i), 'Normalization','pdf'); 
    plot(-1:0.01:3,normpdf(-1:0.01:3,pmpidlogk(i),pspidlogk(i)))
    subplot(3,1,3)
    hold on
    histogram(samples.pilogS2(:,i), 'Normalization','pdf'); 
    plot(-1:0.01:1,normpdf(-1:0.01:1,pmpilogS2(i),pspilogS2(i)))
end

Posterior_summary = table(pmpilogkw,pspilogkw,pmpidlogk,pspidlogk,pmpilogS2,pspilogS2);

writetable(Posterior_summary,'Tables/Posterior_summary.csv','Delimiter',',','QuoteStrings',false)

clear pmpilogkw pspilogkw pmpilogka pspilogka pmpilogS2 pspilogS2
%% Prior predictive check (not used)
% Visual predictive check:
figure('Color', [1 1 1]);
VPC(samplesp.logkCond', datastruct.logkObs, datastruct.fi, 1)
title('Prior predicitve check')

% Individaul and population predictions:
logkCond_p = prctile((samplesp.logkCond),[5 50 95],1);
logkPred_p = prctile((samplesp.logkPred),[5 50 95],1);
figure('Color', [1 1 1]);
rng(3333)
k = datasample(1:1026,10,'Replace',false);

for i = 1:10
    subplot(5, 2, i)
	hold on
	plot(datastruct.fi(datastruct.analyte==k(i)), datastruct.logkObs(datastruct.analyte==k(i)),  'k.')
	plot(datastruct.fi(datastruct.analyte==k(i)), logkCond_p(1,datastruct.analyte==k(i)),  'k:')
    plot(datastruct.fi(datastruct.analyte==k(i)), logkCond_p(2,datastruct.analyte==k(i)),  'k-')
    plot(datastruct.fi(datastruct.analyte==k(i)), logkCond_p(3,datastruct.analyte==k(i)),  'k:')
    plot(datastruct.fi(datastruct.analyte==k(i)), logkPred_p(1,datastruct.analyte==k(i)),  'r:')
    plot(datastruct.fi(datastruct.analyte==k(i)), logkPred_p(2,datastruct.analyte==k(i)),  'r-')
    plot(datastruct.fi(datastruct.analyte==k(i)), logkPred_p(3,datastruct.analyte==k(i)),  'r:')
    set(gca,'XTick',0:0.2:1)
    xlim([0 1])
	ylim([-6 6])	
    if i==5
        ylabel('logk')
    end
    if (i==9) || (i==10)
         xlabel('$$\varphi$$ (ACN)','Interpreter','latex')
    end 
   title(analyte_names.Analyte(k(i)),'FontSize',8)
end

savefig('Figures/PriorPredictionsModel.fig')
set(gcf,'paperunits','centimeters','paperposition',[0 0 16.5 18])
print -dtiff -r300 Figures/PriorPredictionsModel.tif

clear i logkPred_p logkCond_p
%% Goodness of Fit Plots, GOF
% Figure S7. Goodness-of-fit plots. The observed vs. the mean
% population-predicted retention factors (i.e., the a posteriori means of
% predictive distributions corresponding to the future observations of a
% new analyte) and the observed vs the mean individual-predicted retention
% times (i.e., the a posteriori mean of a predictive distribution
% conditioned on the observed data from the same analyte).
logkPred_mean  = mean(samples.logkPred);
logkCond_mean  = mean(samples.logkCond);
figure('Color', [1 1 1]);
subplot(2,1,1)
hold on
plot(logkPred_mean,datastruct.logkObs','.')
xlabel('Population predicted logk')
ylabel('Observed logk')
plot(xlim,xlim,':')
subplot(2,1,2)
hold on
plot(logkCond_mean,datastruct.logkObs','.')
plot(xlim,xlim,':')
xlabel('Individual Predicted logk')
ylabel('Observed logk')
clear logkPred_mean logkCond_mean

savefig('Figures/FigureS8_GOF.fig')
set(gcf,'paperunits','centimeters','paperposition',[0 0 16.5 18])
print -dtiff -r300 Figures/FigureS8_GOF.tif
%% Trace plots (not used)
samples_np = fit101.extract('permuted',false);
Param = 'beta';

[~,m]=size(samples_np(1).(Param));

for i=1:min(m,10)
figure('Color', [1 1 1]);
for z=1:4
hold on
plot(samples_np(z).(Param)(:,i),'-');
xlabel('Iteration')
end
ylabel([Param '(:,' num2str(i) ')'],'fontsize',12);
end

clear samples_np Param z i n m
%% Prior Posterior comparisons (not used):
Param = 'pilogkw';
[~,m]=size(samplesp.(Param));

for i=1:min(m,10)
figure('Color', [1 1 1]);
hold on
[f,xi] = ksdensity((samplesp.(Param)(:,i))); 
plot(xi,f);
[f,xi] = ksdensity((samples.(Param)(:,i)));  
plot(xi,f);
legend('Prior','Posterior')
xlabel([Param '(:,' num2str(i) ')'],'fontsize',12);
ylabel('Probability density estimate')
end

clear Param z i n m f xi 
%% Individaul and population predictions:
% Figure 3.  Individual and population predictions represented as posterior
% medians (lines) and 5th-95th percentiles (dotted lines) for a random set
% of 10 analytes. Observed retention factors are shown as dots. Black
% corresponds to future observations on the same analyte, and red
% corresponds to future observations of a new analyte.

logkCond_p = prctile((samples.logkCond),[5 50 95],1);
logkPred_p = prctile((samples.logkPred),[5 50 95],1);
figure('Color', [1 1 1]);
rng(3333)
k = datasample(1:1026,10,'Replace',false);

for i = 1:10
    subplot(5, 2, i)
	hold on
	plot(datastruct.fi(datastruct.analyte==k(i)), datastruct.logkObs(datastruct.analyte==k(i)),  'k.')
	plot(datastruct.fi(datastruct.analyte==k(i)), logkCond_p(1,datastruct.analyte==k(i)),  'k:')
    plot(datastruct.fi(datastruct.analyte==k(i)), logkCond_p(2,datastruct.analyte==k(i)),  'k-')
    plot(datastruct.fi(datastruct.analyte==k(i)), logkCond_p(3,datastruct.analyte==k(i)),  'k:')
    plot(datastruct.fi(datastruct.analyte==k(i)), logkPred_p(1,datastruct.analyte==k(i)),  'r:')
    plot(datastruct.fi(datastruct.analyte==k(i)), logkPred_p(2,datastruct.analyte==k(i)),  'r-')
    plot(datastruct.fi(datastruct.analyte==k(i)), logkPred_p(3,datastruct.analyte==k(i)),  'r:')
    set(gca,'XTick',0:0.2:1)
    xlim([0 1])
	ylim([-2.5 4])	
    if i==5
        ylabel('logk')
    end
    if (i==9) || (i==10)
         xlabel('$$\varphi$$ (ACN)','Interpreter','latex')
    end 
   title(analyte_names.Analyte(k(i)),'FontSize',8)
end

savefig('Figures/FigureS3_PredictionsModel.fig')
set(gcf,'paperunits','centimeters','paperposition',[0 0 16.5 18])
print -dtiff -r300 Figures/FigureS3_PredictionsModel.tif

clear i logkPred_p logkCond_p
%% 10-fold crossvalidation. Include signle observation per analyte (logk about 1)
% % Figure S8A. Predictions represented as posterior medians (lines) and
% 5th-95th percentiles (dotted lines) for a random set of 10 analytes.
% Observed retention factors are shown as dots. Predictions corresponding
% to future observations given single retention time measurements are shown
% as an open symbol.

logkCond_p = prctile((samplescv1.logkCond),[5 50 95],1);

rng(3333)
k = datasample(1:1026,10,'Replace',false);
datastructcv1.fiEst = datastructcv1.fi(datastructcv1.cvidx);
datastructcv1.analyteEst = datastructcv1.analyte(datastructcv1.cvidx);

figure('Color', [1 1 1]);
for i = 1:length(k)
    subplot(5, 2, i)
	hold on
	plot(datastructcv1.fi(datastructcv1.analyte==k(i)), datastructcv1.logkObs(datastructcv1.analyte==k(i)),  'k.')
    plot(datastructcv1.fiEst(datastructcv1.analyteEst==k(i)), datastructcv1.logkObsEst(datastructcv1.analyteEst==k(i)),  'ko')
	plot(datastructcv1.fi(datastructcv1.analyte==k(i)), logkCond_p(1,datastructcv1.analyte==k(i)),  'k:')
    plot(datastructcv1.fi(datastructcv1.analyte==k(i)), logkCond_p(2,datastructcv1.analyte==k(i)),  'k-')
    plot(datastructcv1.fi(datastructcv1.analyte==k(i)), logkCond_p(3,datastructcv1.analyte==k(i)),  'k:')

    set(gca,'XTick',0:0.2:1)
    xlim([0 1])
	ylim([-2.5 4])	
    if i==5
        ylabel('logk')
    end
    if (i==9) || (i==10)
        xlabel('$$\varphi$$ (ACN)','Interpreter','latex')
    end 
   title(analyte_names.Analyte(k(i)),'FontSize',8)
end

savefig('Figures/FigureS8A_PredictionsCV1.fig')
set(gcf,'paperunits','centimeters','paperposition',[0 0 16.5 18])
print -dtiff -r300 Figures/FigureS8A_PredictionsCV1.tif

clear i logkPred_p logkCond_p
%% 10-fold crossvalidation. Include signle observation per analyte with min(logk)
% Figure S8B. Predictions represented as posterior medians (lines) and
% 5th-95th percentiles (dotted lines) for a random set of 10 analytes.
% Observed retention factors are shown as dots. Predictions corresponding
% to future isocratic observations given single retention time measurements
% are shown as an open symbol.
logkCond_p = prctile((samplescv2.logkCond),[5 50 95],1);

rng(3333)
k = datasample(1:1026,10,'Replace',false);
datastructcv2.fiEst = datastructcv2.fi(datastructcv2.cvidx);
datastructcv2.analyteEst = datastructcv2.analyte(datastructcv2.cvidx);

figure('Color', [1 1 1]);
for i = 1:length(k)
    subplot(5, 2, i)
	hold on
	plot(datastructcv2.fi(datastructcv2.analyte==k(i)), datastructcv2.logkObs(datastructcv2.analyte==k(i)),  'k.')
    plot(datastructcv2.fiEst(datastructcv2.analyteEst==k(i)), datastructcv2.logkObsEst(datastructcv2.analyteEst==k(i)),  'ko')
	plot(datastructcv2.fi(datastructcv2.analyte==k(i)), logkCond_p(1,datastructcv2.analyte==k(i)),  'k:')
    plot(datastructcv2.fi(datastructcv2.analyte==k(i)), logkCond_p(2,datastructcv2.analyte==k(i)),  'k-')
    plot(datastructcv2.fi(datastructcv2.analyte==k(i)), logkCond_p(3,datastructcv2.analyte==k(i)),  'k:')
    
    set(gca,'XTick',[0:0.2:1])
    xlim([0 1])
	ylim([-2.5 4])	
    if i==5
        ylabel('logk')
    end
    if (i==9) || (i==10)
        xlabel('$$\varphi$$ (ACN)','Interpreter','latex')
    end 
   title(analyte_names.Analyte(k(i)),'FontSize',8)
end

savefig('Figures/FigureS8B_PredictionsCVmin.fig')
set(gcf,'paperunits','centimeters','paperposition',[0 0 16.5 18])
print -dtiff -r300 Figures/FigureS8B_PredictionsCVmin.tif

clear fi i

%% 10-fold crossvalidation. Include signle observation logk at fi=0.3.
% In not observed the observation is assummed to be censored 

% % % Figure S8C. Predictions represented as posterior median (line) and
% 5th-95th percentiles (dotted lines) for a random set of 10 analytes.
% Observed retention factors are shown as dots. Predictions corresponding
% to future isocratic observations given single retention time measurements
% are shown as an open symbol.

logkCond_p = prctile((samplescv3.logkCond),[5 50 95],1);

rng(3333)
k = datasample(1:1026,10,'Replace',false);
datastructcv3.fiEst = datastructcv3.fi(datastructcv3.cvidx);
datastructcv3.analyteEst = datastructcv3.analyte(datastructcv3.cvidx);

figure('Color', [1 1 1]);
for i = 1:length(k)
    subplot(5, 2, i)
	hold on
	plot(datastructcv3.fi(datastructcv3.analyte==k(i)), datastructcv3.logkObs(datastructcv3.analyte==k(i)),  'k.')
    plot(datastructcv3.fiEst(datastructcv3.analyteEst==k(i)), datastructcv3.logkObsEst(datastructcv3.analyteEst==k(i)),  'ko')
	plot(datastructcv3.fi(datastructcv3.analyte==k(i)), logkCond_p(1,datastructcv3.analyte==k(i)),  'k:')
    plot(datastructcv3.fi(datastructcv3.analyte==k(i)), logkCond_p(2,datastructcv3.analyte==k(i)),  'k-')
    plot(datastructcv3.fi(datastructcv3.analyte==k(i)), logkCond_p(3,datastructcv3.analyte==k(i)),  'k:')
    
    set(gca,'XTick',[0:0.2:1])
    xlim([0 1])
	ylim([-2.5 4])	
    if i==5
        ylabel('logk')
    end
    if (i==9) || (i==10)
        xlabel('$$\varphi$$ (ACN)','Interpreter','latex')
    end 
   title(analyte_names.Analyte(k(i)),'FontSize',8)
end

savefig('Figures/FigureS8C_PredictionsCVfix.fig')
set(gcf,'paperunits','centimeters','paperposition',[0 0 16.5 18])
print -dtiff -r300 Figures/FigureS8C_PredictionsCVfix.tif
clear fi i
%% Posterior predictive checks (not used)
figure('Color', [1 1 1]);
subset = datastruct.Mmolx>-inf;
VPC(samples.logkPred(:,subset)', datastruct.logkObs(subset), datastruct.fi(subset),0)
clear subset
%% Influence of functional groups
% Figure 2. Graphical display of the marginal posterior distributions for
% the effects of each functional group on logkw, logka, and logS2.
figure('Color', [1 1 1]);
subplot(1,4,2)
hold on
boxplot_pwhisker(samples.pilogkw(:,:),{'Labels',functional_groups_names{:,2}},5,95);
plot(xlim,[0 0],':')
ylim([0 3])
view(90,90)
set(gca,'FontSize',5)
set(gca,'Position', [0.2139    0.1100    0.2138    0.8150])
ylabel('\pi_{logk_w}','FontSize',8)
subplot(1,4,3)
hold on
boxplot_pwhisker(samples.pilogka(:,:),{'Labels',functional_groups_names{:,1}},5,95);
plot(xlim,[0 0],':')
ylim([-1 3])
view(90,90)
set(gca,'FontSize',5)
set(gca,'Position', [0.4854    0.1100    0.2178    0.8150])
ylabel('\pi_{logk_a}','FontSize',8)
subplot(1,4,4)
hold on
boxplot_pwhisker(samples.pilogS2(:,:),{'Labels',functional_groups_names{:,1}},5,95);
plot(xlim,[0 0],':')
ylim([-0.5 0.5])
view(90,90)
set(gca,'FontSize',5)
set(gca,'Position', [0.7334    0.1100    0.1708    0.8150])
ylabel('\pi_{logS_2}','FontSize',8)

savefig('Figures/Figure2_FunctionalGroupEffects.fig')
set(gcf,'paperunits','centimeters','paperposition',[0 0 16.5 18])
print -dtiff -r300 Figures/Figure2_FunctionalGroupEffects.tif

%% Influence of functional groups (pilogkw)
% Figure S2A. Graphical display of the marginal posterior distributions for
% the effects of each functional groups on pilogkw.
figure('Color', [1 1 1]);

[~,I] = sort(mean(samples.pilogkw(:,:)),'descend');

hold on
boxplot_pwhisker(samples.pilogkw(:,I),{'Labels',functional_groups_names{I,2}},5,95);
plot(xlim,[0 0],':')
ylim([0 3])
view(90,90)
set(gca,'FontSize',5)
set(gca,'Position', [0.2343    0.1100    0.6707    0.8150])
ylabel('\pi_{logk_w}','FontSize',8)

savefig('Figures/Figure2A_FunctionalGroupEffects_pilogkw.fig')
set(gcf,'paperunits','centimeters','paperposition',[0 0 16.5 18])
print -dtiff -r300 Figures/Figure2A_FunctionalGroupEffects_pilogkw.tif
%% Influence of functional groups (pilogka)
% Figure S2B. Graphical display of the marginal posterior distributions for
% the effects of each functional groups on pilogka.
figure('Color', [1 1 1]);

[~,I] = sort(mean(samples.pilogka(:,:)),'descend');

hold on
boxplot_pwhisker(samples.pilogka(:,I),{'Labels',functional_groups_names{I,2}},5,95);
plot(xlim,[0 0],':')
ylim([-1 3])
view(90,90)
set(gca,'FontSize',5)
set(gca,'Position', [0.2343    0.1100    0.6707    0.8150])
ylabel('\pi_{logk_a}','FontSize',8)

savefig('Figures/Figure2B_FunctionalGroupEffects_pilogka.fig')
set(gcf,'paperunits','centimeters','paperposition',[0 0 16.5 18])
print -dtiff -r300 Figures/Figure2B_FunctionalGroupEffects_pilogka.tif

%% Influence of functional groups (pilogS2)
% Figure S2C. Graphical display of the marginal posterior distributions for
% the effects of each functional groups on pilogS2.

figure('Color', [1 1 1]);
[~,I] = sort(mean(samples.pilogS2(:,:)),'descend');

hold on
boxplot_pwhisker(samples.pilogS2(:,I),{'Labels',functional_groups_names{I,2}},5,95);
plot(xlim,[0 0],':')
ylim([-0.2 0.2])
view(90,90)
set(gca,'FontSize',5)
set(gca,'Position', [0.2343    0.1100    0.6707    0.8150])
ylabel('\pi_{logS2}','FontSize',8)

savefig('Figures/Figure2C_FunctionalGroupEffects_pilogS2.fig')
set(gcf,'paperunits','centimeters','paperposition',[0 0 16.5 18])
print -dtiff -r300 Figures/Figure2C_FunctionalGroupEffects_pilogS2.tif
%% Influence of functional groups (pidlogk)
% Figure S2D. Graphical display of the marginal posterior distributions for the effects
% of each functional groups on the difference between pilogkw and pilogka.
figure('Color', [1 1 1]);

[~,I] = sort(mean(samples.pidlogk(:,:)),'descend');

hold on
boxplot_pwhisker(samples.pidlogk(:,I),{'Labels',functional_groups_names{I,2}},5,95);
plot(xlim,[0 0],':')
ylim([-1 2])
view(90,90)
set(gca,'FontSize',5)
set(gca,'Position', [0.2343    0.1100    0.6707    0.8150])
ylabel('\pi_{dlogk}=\pi_{logk_w}=\pi_{logk_a}','FontSize',8)

savefig('Figures/Figure2D_FunctionalGroupEffects_pidlogk.fig')
set(gcf,'paperunits','centimeters','paperposition',[0 0 16.5 18])
print -dtiff -r300 Figures/Figure2D_FunctionalGroupEffects_pidlogk.tif

%% Individual parameters and molecular mass
% Figure S3A. Scatter plots between individual chromatographic parameters
% and molecular mass. Diagonal subplots present histograms.
idata.param = squeeze(mean(samples.param));           % logkw, logka, logS2
idata.eta = squeeze(mean(samples.param-samples.miu)); % logkw, logka, logS2
idata.Mmolx = datastruct.Mmolx(datastruct.start)';
figure('Color', [1 1 1]);
xynames = {'log k_w','log k_a','log S_{2}','(Mmol-300)/100'};
gplotmatrix([idata.param(:,[1 2 3]) idata.Mmolx],[],0*idata.Mmolx,'kk',[],[],'on','stairs',xynames,xynames)

h=get(gcf,'children');
set(h(1),'Visible','off')
savefig('Figures/FigureS3A_IndividualParametersPlots.fig')
set(gcf,'paperunits','centimeters','paperposition',[0 0 16.5 18])
print -dtiff -r300 Figures/FigureS3A_IndividualParametersPlots.tif

clear xynames h

%% Etas and molecular mass
% Figure S3B. The scatter plots between eta values (difference between the
% analyte-specific chromatographic parameter and expected value). Diagonal
% subplots present histograms.
idata.param = squeeze(mean(samples.param));           % logkw, logka, logS2
idata.eta = squeeze(mean(samples.param-samples.miu)); % logkw, logka, logS2
idata.Mmolx = datastruct.Mmolx(datastruct.start)';
figure('Color', [1 1 1]);
xynames = {'\eta_{log k_w}','\eta_{log k_a}','\eta_{log S_2}','(Mmol-300)/100'};
gplotmatrix([idata.eta idata.Mmolx],[],0*idata.Mmolx,'kk',[],[],'on','stairs',xynames,xynames)
clear xynames

h=get(gcf,'children');
set(h(1),'Visible','off')

savefig('Figures/FigureS3B_EtaPlots.fig')
set(gcf,'paperunits','centimeters','paperposition',[0 0 16.5 18])
print -dtiff -r300 Figures/FigureS3B_EtaPlots.tif

clear h
%% Parameters vs. Mmol
% Figure S4. Effect of molecular mass on retention of compounds without
% functional groups (10 draws of the expected line). The analyte-specific
% chromatographic parameters are lower due to the presence of functional
% groups. For logS2, the effects are in both directions. A congeneric
% series of alkyl-substituted benzamides (benzamide and N-methyl-, N-ethyl,
% …, N-pentadecyl-, N-hexadecyl benzamides) is shown to highlight the
% effect of molecular mass.
idata.param = squeeze(mean(samples.param));  % logkw, logka, logS2
idata.Mmolx = datastruct.Mmolx(datastruct.start)';

figure('Color', [1 1 1]);
subplot(3,1,1)
hold on
plot(idata.Mmolx, idata.param(:,1),'.')
plot(idata.Mmolx(99:115), idata.param(99:115,1),'ro')
for i=1:100
plot([-2.3:0.1:3.6],samples.logkwHat(i) + samples.beta(i,1)*[-2.3:0.1:3.6],'-')
end
box off
ylabel('log k_w')

subplot(3,1,2)
hold on
plot(idata.Mmolx, idata.param(:,2),'.')
plot(idata.Mmolx(99:115), idata.param(99:115,2),'ro')
for i=1:100
plot([-2.3:0.1:3.6],samples.logkaHat(i) + samples.beta(i,2)*[-2.3:0.1:3.6],'-')
end
ylabel('log k_a')
box off
subplot(3,1,3)
hold on
plot(idata.Mmolx, idata.param(:,3),'.')
plot(idata.Mmolx(99:115), idata.param(99:115,3),'ro')
for i=1:100
plot([-2.3:0.1:3.6],samples.logS2Hat(i) + samples.beta(i,3)*[-2.3:0.1:3.6],'-')
end
xlabel('(Mmol-300)/100')
ylabel('log S_{2}')
box off


savefig('Figures/FigureS4_Param_MolecularMass.fig')
set(gcf,'paperunits','centimeters','paperposition',[0 0 16.5 18])
print -dtiff -r300 Figures/FigureS4_Param_MolecularMass.tif
%% Shrinkage
% Figure S5. Comparison of model parameters obtained using the two-stage
% approach and multilevel model. Shrinkage is shown as lines connecting the
% parameters estimated by both methods. The multilevel model shows more
% stable estimates, particularly for parameters that are difficult to
% estimate precisely due to lack of information.
idata.param = squeeze(mean(samples.param));  % logkw, logka, logS2
idata.Mmolx = datastruct.Mmolx(datastruct.start)';

param = [initial_paramA(:,2) initial_paramA(:,2) + initial_paramA(:,1)./3 log10(initial_paramA(:,3))];

figure('Color', [1 1 1]);
subplot(2,1,1)
hold on
quiver(idata.Mmolx, param(:,1),idata.Mmolx-idata.Mmolx, idata.param(:,1)-param(:,1),0,'MaxHeadSize',0.0)
plot(idata.Mmolx, idata.param(:,1),'k.')
plot(idata.Mmolx, param(:,1),'r.')
% legend({'multilevel model','two stage approach'}, 'Location','Northwest')
% hold on
% box off
ylabel('log k_w')

subplot(2,1,2)
hold on
quiver(idata.Mmolx, param(:,2),idata.Mmolx-idata.Mmolx, idata.param(:,2)-param(:,2),0,'MaxHeadSize',0.)
plot(idata.Mmolx, idata.param(:,2),'k.')
plot(idata.Mmolx, param(:,2),'r.')
ylabel('log k_a')
box off
xlabel('(Mmol-300)/100')

savefig('Figures/FigureS5_Shrinkage.fig')
set(gcf,'paperunits','centimeters','paperposition',[0 0 16.5 18])
print -dtiff -r300 Figures/FigureS5_Shrinkage.tif

%% Histogram of mean values of functional groups effects
% Figure S6. The histogram of mean posterior values of the effects of each
% functional group on chromatographic parameters.

figure('Color', [1 1 1]);
subplot(3,1,1)
hold on
for i=1:10
plot(0:0.01:3,lognpdf(0:0.01:3,log(samples.mpilogkw(i)),samples.spilogkw(i)))
end
histogram(mean(samples.pilogkw(:,:)), 20, 'Normalization', 'pdf')
xlabel('\pi_{logk_w}')
subplot(3,1,2)
hold on
for i=1:10
pd = makedist('tLocationScale','mu',samples.mpidlogk(i),'sigma',samples.spidlogk(i),'nu',samples.nupi(i));
plot(-2:0.01:2,pdf(pd,-2:0.01:2))
end
histogram(mean(samples.pidlogk(:,:)),20, 'Normalization', 'pdf')
xlabel('\pi_{logk_w}-\pi_{logk_a}')
ylabel('Frequency')
subplot(3,1,3)
hold on
for i=1:10
plot(-0.5:0.01:0.5,normpdf(-0.5:0.01:0.5,0,samples.spilogS2(i)))
end
histogram(mean(samples.pilogS2(:,:)),20, 'Normalization', 'pdf')
xlabel('\pi_{logS_2}')

savefig('Figures/FigureS6_FunctionalGroupEffectsHistogram.fig')
set(gcf,'paperunits','centimeters','paperposition',[0 0 16.5 18])
print -dtiff -r300 Figures/FigureS6_FunctionalGroupEffectsHistogram.tif
%% logkw - logka comparisons (not used in the manuscript)
figure('Color', [1 1 1]);
plot([mean(samples.pilogkw(:,:))' mean(samples.pilogka(:,:))']')
ylabel('\pi_{log k_w}')
yyaxis right
plot([2 2], [-1 2.5])
ylabel('\pi_{log k_a}');
box off
set(gca,'XTick',[1 2],'XTickLabel',{'\pi_{log k_w}','\pi_{log k_a}'});
savefig('Figures/logkw_logka_comparisons.fig')
set(gcf,'paperunits','centimeters','paperposition',[0 0 16.5 18])
print -dtiff -r300 Figures/logkw_logka_comparisons.tif

%% Plot uncertainity chromatograms expected isocratically at fi=0.5 for 10 randomly selected analytes
% Figure 4. Uncertainty chromatograms for isocratic conditions. It
% represents the range of analyte retention time along with uncertainty as
% predicted by the proposed model conditional on different experimental
% data. Colors correspond to different analytes that are identified at the
% bottom figure: 112: N-tridecylbenzamide, 122: Tetrabutylammonium,  241:
% Metaflumizone, 379: Apigenin, 498: CGS-21680 hydrochloride, 512:
% 6,7-dinitro-1,4-dihydroquinoxaline-2,3-dione, 626: Lidocaine N-ethyl
% bromide quaternary salt, 672:  Oxybutynin Chloride, 726: Ro 04-6790
% dihydrochloride, 772: Tolbutamide
figure('Color', [1 1 1]);
rng(3333)
k = datasample(1:length(unique(data.ID)),100,'Replace',false);

logkfi1 = hplcmodelizo(0.5, samples.paramPred(:,k,:));
logkfi2 = hplcmodelizo(0.5, samplescv1.param(:,k,:)); 
logkfi3 = hplcmodelizo(0.5, samplescv2.param(:,k,:));  
logkfi4 = hplcmodelizo(0.5, samplescv3.param(:,k,:)); 
logkfi5 = hplcmodelizo(0.5, samples.param(:,k,:)); 

for i=1:1:10
subplot(5,1,1)
hold on
[f,xi] = ksdensity(logkfi1(:,i));
plot(xi,f);
xlim([-3 3])
set(gca,'XTick',[-3:1:3], 'XTickLabel', [])
title('no preliminary experiments','FontWeight','normal')

subplot(5,1,2)
hold on
[f,xi] = ksdensity(logkfi2(:,i));
plot(xi,f);
xlim([-3 3])
set(gca,'XTick',[-3:1:3], 'XTickLabel', [])
title('one preliminary experiment (logk about 1)','FontWeight','normal')

subplot(5,1,3)
hold on
[f,xi] = ksdensity(logkfi3(:,i));
plot(xi,f);
xlim([-3 3])
ylabel('Uncertainity chromatogram, probability density estimate')
set(gca,'XTick',[-3:1:3], 'XTickLabel', [])
title('one preliminary experiment (min observed logk)','FontWeight','normal')

subplot(5,1,4)
hold on
[f,xi] = ksdensity(logkfi4(:,i));
plot(xi,f);
xlim([-3 3])
set(gca,'XTick',[-3:1:3], 'XTickLabel', [])
title('one preliminary experiment (fi=0.3)','FontWeight','normal')

subplot(5,1,5)
hold on
[f,xi] = ksdensity(logkfi5(:,i));
plot(xi,f);
xlim([-3 3])
xlabel('Retention factor, k')
[a,b] = max(f);
set(gca,'XTick',[-3:1:3], 'XTickLabel', {0.001 0.01 0.1 1 10 100 1000}, 'XTickLabelRotation', 45)
text(xi(b),5+a.*1.1,num2str(k(i)),'HorizontalAlignment', 'right','FontSize', 8)
title('all experimental data','FontWeight','normal')

end
% legend(analyte_names.Analyte(k(1:10)))

analyte_names(k(1:10),:)

savefig('Figures/Figure4_UncertainityChromatogram_fi_0p5.fig')
set(gcf,'paperunits','centimeters','paperposition',[0 0 16.5 18])
print -dtiff -r300 Figures/Figure4_UncertainityChromatogram_fi_0p5.tif
clear a b logkfi1 logkfi2 logkfi3 logkfi4 xi

%% Plot uncertainity chromatograms expected for a typical gradient for 10 randomly selected analytes
% Figure 5. Uncertainty chromatograms for selected gradient conditions. It
% represents the range of analyte retention time along with uncertainty as
% predicted by the proposed model conditional on different experimental
% data. Colors correspond to different analytes that are identified at the
% bottom figure: 112: N-tridecylbenzamide, 122: Tetrabutylammonium,  241:
% Metaflumizone, 379: Apigenin, 498: CGS-21680 hydrochloride, 512:
% 6,7-dinitro-1,4-dihydroquinoxaline-2,3-dione, 626: Lidocaine N-ethyl
% bromide quaternary salt, 672:  Oxybutynin Chloride, 726: Ro 04-6790
% dihydrochloride, 772: Tolbutamide.
figure('Color', [1 1 1]);
rng(3333)
k = datasample(1:length(unique(data.ID)),100,'Replace',false);

 tr1 = hplcmodelgra(samples.paramPred(:,k(1:10),:),20,0,1,1,0.2,0,linspace(0,20+10,2000));
 tr1(tr1>30)=30;
 tr2 = hplcmodelgra(samplescv1.param(:,k(1:10),:),20,0,1,1,0.2,0,linspace(0,20+10,2000));
 tr3 = hplcmodelgra(samplescv2.param(:,k(1:10),:),20,0,1,1,0.2,0,linspace(0,20+10,2000));
 tr4 = hplcmodelgra(samplescv3.param(:,k(1:10),:),20,0,1,1,0.2,0,linspace(0,20+10,2000));
 tr4(tr4>30)=30;
 tr5 = hplcmodelgra(samples.param(:,k(1:10),:),20,0,1,1,0.2,0,linspace(0,20+10,2000)); 

for i=1:1:10 
subplot(5,1,1)
hold on
[f,xi] = ksdensity(tr1(:,i));
plot(xi,f);
xlim([0 30])
set(gca,'XTick',[0:5:30], 'XTickLabel', [])
title('no preliminary experiments','FontWeight','normal')

subplot(5,1,2)
hold on
[f,xi] = ksdensity(tr2(:,i));
plot(xi,f);
xlim([0 30])
set(gca,'XTick',[0:5:30], 'XTickLabel', [])
title('one preliminary experiment (logk about 1)','FontWeight','normal')

subplot(5,1,3)
hold on
[f,xi] = ksdensity(tr3(:,i));
plot(xi,f);
xlim([0 30])
set(gca,'XTick',[0:5:30], 'XTickLabel', [])
ylabel('Uncertainity chromatogram, probability density estimate')
title('one preliminary experiment (min observed logk)','FontWeight','normal')

subplot(5,1,4)
hold on
[f,xi] = ksdensity(tr4(:,i));
plot(xi,f);
xlim([0 30])
set(gca,'XTick',[0:5:30], 'XTickLabel', [])
title('one preliminary experiment (fi=0.3)','FontWeight','normal')

subplot(5,1,5)
hold on
[f,xi] = ksdensity(tr5(:,i));
plot(xi,f);
xlim([0 30])
ylim([0 20])
set(gca,'XTick',[0:5:30])
xlabel('Retention time, t_R')
[a,b] = max(f);
text(xi(b),0+a.*1.2,num2str(k(i)),'HorizontalAlignment', 'left','FontSize', 8)
title('all experimental data','FontWeight','normal')

end
% legend(analyte_names.Analyte(k(1:10)))

analyte_names(k(1:10),:)

savefig('Figures/Figure5_UncertainityChromatogram_gra.fig')
set(gcf,'paperunits','centimeters','paperposition',[0 0 16.5 18])
print -dtiff -r300 Figures/Figure5_UncertainityChromatogram_gra.tif
% clear a b tr1 tr2 tr3 tr4 xi f
%% TOC


%% Version
ver
% ----------------------------------------------------------------------------------------------------
% MATLAB Version: 9.2.0.556344 (R2017a)
% MATLAB License Number: 261217
% Operating System: Microsoft Windows 10 Pro Version 10.0 (Build 19041)
% Java Version: Java 1.7.0_60-b19 with Oracle Corporation Java HotSpot(TM) 64-Bit Server VM mixed mode
% ----------------------------------------------------------------------------------------------------
% MATLAB                                                Version 9.2         (R2017a)
% Bioinformatics Toolbox                                Version 4.8         (R2017a)
% Curve Fitting Toolbox                                 Version 3.5.5       (R2017a)
% Global Optimization Toolbox                           Version 3.4.2       (R2017a)
% MATLAB Compiler                                       Version 6.4         (R2017a)
% MATLAB Compiler SDK                                   Version 6.3.1       (R2017a)
% Optimization Toolbox                                  Version 7.6         (R2017a)
% Parallel Computing Toolbox                            Version 6.10        (R2017a)
% SimBiology                                            Version 5.6         (R2017a)
% Statistics and Machine Learning Toolbox               Version 11.1        (R2017a)
% Symbolic Math Toolbox                                 Version 7.2         (R2017a)

%% Licenses 
%1) Code &copy; 2020, Pawe³ Wiczling, licensed under BSD-3.
%2) Text &copy; 2020, Pawe³ Wiczling, licensed under CC-BY-NC 4.0.