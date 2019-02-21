%% Analyses for Experiment 4 (continuous confidence data)

clear
clc
close all

% Subjects
subject_id = [1:13,15,16,18:22];
%subject_id = 10;
subjects = 1:length(subject_id);
numSessions = 3;
nratings = 2;

%flags
compute = 1;
save_measures = 1;

% Add helper functions
root = pwd;
parts = strsplit(root, '/');
cd(root(1:end-length(parts{end})));
addpath('helperFunctions','Mfunctions');
cd(root)

%Add data folder to path
addpath('Data')

if compute == 1
    
    for subject = subjects
        clear stim resp correct conf contrast sessionID
        
        totalTrials = 0;
        for session = 1:numSessions
            load(['results_sub',num2str(subject_id(subject)),'_s',num2str(session),'.mat'])
            [m,n] = size(params.data);
            numBlocks = m*n;
            main_blocks = [1:numBlocks];
            block_length = length(params.data{1,1}.stim);
            
            for block=1:length(main_blocks)
                trials = totalTrials + [1:block_length];
                totalTrials = totalTrials + block_length;
                stim(trials) = params.data{main_blocks(block)}.stim-1; %1: left, 2: right
                resp(trials) = params.data{main_blocks(block)}.answer-1;%1: left, 2: right
                correct(trials) = params.data{main_blocks(block)}.correct; %0: error, 1: correct
                conf_cont(trials) = params.data{main_blocks(block)}.conf; %0.5 to 1
                contrast(trials) = params.data{main_blocks(block)}.contrast; %1-3
                sessionID(trials) = session;
            end
        end
        
        %Vary criteria from .51 to .99
        criterion_vect = [.51:.01:.99];
        nCriteria = length(criterion_vect);
        
        for condition = 1:3
            for criterion = 1:length(criterion_vect)
                
                %Transform confidence into binary variable (1- low, 2-high)
                conf = (conf_cont>criterion_vect(criterion))+1;
                
                %Compute metacognition measures
                %Phi
                R = corrcoef(correct(contrast==condition),conf(contrast==condition));
                phi(subject,condition,criterion) = R(2);
                %Type-2 AUC
                [nR_S1, nR_S2] = trials2counts(stim(contrast==condition),resp(contrast==condition),conf(contrast==condition),nratings,[]);
                [t2AUC(subject,condition,criterion)] = type2ag(nR_S1, nR_S2, 1);
                %Mratio and meta-d'
                meta_MLE = type2_SDT_MLE(stim(contrast==condition),resp(contrast==condition),conf(contrast==condition),nratings, [], 1);
                mratio(subject,condition,criterion) = meta_MLE.M_ratio;
                meta_d(subject,condition,criterion) = meta_MLE.meta_da;
                %Type-2 d'
                t2_dprime(subject,condition,criterion) = data_analysis_resp(correct(contrast==condition),conf(contrast==condition)-1);
                
                %Compute HR and FAR for each criterion location
                [FAR_,HR_] = zROC(stim(contrast==condition),resp(contrast==condition),conf(contrast==condition),2,1);
                FAR(subject,nCriteria+1-criterion,condition) = FAR_(1); FAR(subject,nCriteria+1+criterion,condition) = FAR_(3);
                HR(subject,nCriteria+1-criterion,condition) = HR_(1); HR(subject,nCriteria+1+criterion,condition) = HR_(3);
            end
            %Decision criterion
            HR(subject,nCriteria+1,condition) = length(find(correct(contrast==condition)==1 & stim(contrast==condition)==1))./length(find(stim(contrast==condition)==1));
            FAR(subject,nCriteria+1,condition) = length(find(correct(contrast==condition)==0 & stim(contrast==condition)==0))./length(find(stim(contrast==condition)==0));
            
        end
        
    end
end

if save_measures
    save('metaMeasures','phi','t2AUC','mratio','meta_d','t2_dprime');
    save('ROC','HR','FAR')
else
    load metaMeasures
    load ROC
end

%% Figure - 3: Plot measures of metacognition against confidence criterion location

axes('box','off','tickdir','out','LineWidth',1.25,'FontSize',16); hold on;
criterion_vect = [51:1:99];
colors = {[.8 .1 .1],[.4 .6 .4],[.8 .8 1]};

measures = {'mratio','meta_d','t2AUC','phi'};
names = {'meta-d''/d''','meta-d''','Type-2 AUC','phi'};
numMeasures = length(names);

%y-limits
ylims = {[.4,1.15],[.4,2.6],[.5,.725],[.05,.325]};

for measure = 1:4
    y_ = eval(measures{measure});
    
    %Plot data
    subplot(2,2, measure); hold on;
    for condition = 1:3
        for criterion = 1:length(criterion_vect)
            y(criterion) = nanmean(y_(:,condition,criterion));
            y_sem(criterion) = nanstd(y_(:,condition,criterion))./sqrt(20);
        end
        h{condition} = plot(criterion_vect,y,'linewidth',2.5,'color',colors{condition}); hold on;
        shadedplot(criterion_vect,y-y_sem,y+y_sem,colors{condition})
    end
    
    if measure == 4
        l = legend([h{1},h{2},h{3}],'Contrast 1','Contrast 2','Contrast 3');
        set(l,'box','off','orientation','horizontal','fontsize',18)
    end
    
    ylabel(names{measure},'fontsize',18)
    ylim(ylims{measure})
    xlabel('Confidence criterion location','fontsize',18)
end

%% Curve fitting
% Generates linear fits for meta-d'/d' and meta-d' and quadratic fits for type-2 AUC
% and phi

fig_count = 0;

for subject = 1:20
    for condition = 1:3
        figure;
        fig_count = fig_count+1;
        for measure = 1:numMeasures
            y_ = eval(measures{measure});
            
            y = y_(subject,condition,:);
            if measure < 3 % For meta-d'/d' and meta-d'
                clear p p1 p2
                subplot(2,2,measure); hold on;
                % Linear fit
                [fitresult, gof] = PolynomialFit(criterion_vect,y,1);
                %Save the slope and intercept values
                a = fitresult.p1; %Slope
                b = fitresult.p2; %Intercept
                params = [a,b];
            else            %type-2 AUC and phi
                clear p p1 p2
                subplot(2,2,measure)
                %Quadratic fit
                [fitresult, gof] = PolynomialFit(criterion_vect,y,2);
                %Save the polynomial coefficients
                a = fitresult.p1; %degree of curvature - p1 < 0 --> downward curvature
                b = fitresult.p2; %location
                c = fitresult.p3; %y-intercept
                params = [a,b,c];
            end
            
            ylabel(names{measure},'fontsize',18)
            xlabel('Confidence criterion location','fontsize',18)
            curvefits{subject,measure,condition} = params;
        end
        % saveFigure(num2str(fig_count));
        delete gcf
    end
    
end

% Save fitting results
save('measures_curvefits.mat','curvefits')

%% Analysis on curve fits
% T-test for a < 0

load measures_curvefits.mat
for measure = 1:4
    for subject = subjects
        for condition = 1:3
            a(subject,condition) = curvefits{subject,measure,condition}(1)
        end
    end
    
    %One sample t-tests on a to check if a < 0
    [h,params,~,tstats] = ttest(a)
    tstat(measure,:) = tstats.tstat;
    pval(measure,:) = params;
end

%% Effect of task difficulty on metacognition measures

% Effect of contrast(task difficulty) on dprime
load DataForModeling.mat

y = exptData.dprime;
[params,t,stats] = onewayRepmeasuresANOVA(y,'Contrast');
%Save p and F values
F = t{2,6}; P = params(1);
%Pairwise comparisons
[pval,tstat] = pairwiseComparisons(y);

% Effect of task difficulty (contrast) on measures
for measure = 1:4
    clear y
    y_ = eval(measures{measure});
    for condition = 1:3
        for subject = 1:20
            y(subject,condition) = nanmean(y_(subject,condition,:));
        end
    end
    [params,t,stats] = onewayRepmeasuresANOVA(y,'Contrast');
    
    %Save p and F values
    F(measure,1) = t{2,6}; P(measure,1) = params(1);
    %Pairwise comparisons
    [pval(measure,:),tstat(measure,:)] = pairwiseComparisons(y);
    
end

%% Figure 4

% Confidence criteria used for model fitting
conf_binSize = 50/6;
for rating = 1:5
    model_confCriteria(rating) = 50+(rating)*conf_binSize;
    modelCriteria_index(rating) = find(criterion_vect<model_confCriteria(rating), 1, 'last' );
end

axes('box','off','tickdir','out','LineWidth',1.25,'FontSize',16); hold on;
for condition = 1:3
    
    %zHR and zFAR
    zHR = norminv(HR);
    zFAR = norminv(FAR);
    
    %Average across subjects
    zHR_avg = nanmean(zHR(:,:,condition));
    zFAR_avg = nanmean(zFAR(:,:,condition));
    nPoints = length(zHR_avg);
    
    %Plot confidence criteria
    h1(condition) = plot(zFAR_avg, zHR_avg, '.','markersize',25,'color',colors{condition},'linewidth',7);hold on;
    %Plot decision criterion in a different color
    h2 = plot(zFAR_avg(nCriteria+1), zHR_avg(nCriteria+1),'*k','markersize',10,'linewidth',2);
    %Plot the criteria used for model fitting in a separate color
    h3 = plot(zFAR_avg([modelCriteria_index, nCriteria+1 + modelCriteria_index]), zHR_avg([modelCriteria_index, nCriteria+1 + modelCriteria_index]),'x','color',[0 0 0],'markersize',10,'linewidth',3);
    %Plot the predicted linear zROC plot from d' estimated from decision criterion
    dprime = zHR_avg(nCriteria+1)-zFAR_avg(nCriteria+1);
    xlimits = [-3, 2];
    plot(xlimits, dprime+[xlimits], '--','color',[.5 .5 .5],'linewidth',1.5)
    
    xlim(xlimits);
    ylim([-2,3])
    xlabel({'zFAR';''},'fontsize',18);
    ylabel('zHR','fontsize',18);
    
end


l = legend([h2,h3,h1],'Decision criterion','Confidence criteria used for model fitting','Confidence criteria - Contrast 1','Confidence criteria - Contrast 2','Confidence criteria - Contrast 3','location','best');
set(l,'fontsize',18)

%% Generate curve fits for zROCs (after rotating clockwise by 45 deg)
for subject = subjects
    figure;
    for condition = 1:3
        clear p1 p2 p3 p
        zFAR = norminv(FAR(subject,:,condition));
        zHR = norminv(HR(subject,:,condition));
        
        %Rotate by 45 degrees
        [th, r] = cart2pol(zFAR, zHR);
        [nzFAR, nzHR] = pol2cart(th-pi/4, r);
        
        %Check rotation
        %plot(zFAR,zHR,'.b','markersize',30); hold on;
        plot(nzFAR,nzHR,'.k','markersize',30); hold on;
        
        %Polynomial fit
        [fitresult, gof] = PolynomialFit(nzFAR, nzHR, 2);
        %Save the polynomial coefficients
        a = fitresult.p1; %degree of curvature - a < 0 --> curves downwards
        b = fitresult.p2; %location
        c = fitresult.p3; %y-intercept
        params = [a,b,c];
        
        ylabel('zHR','fontsize',18)
        xlabel('zFAR','fontsize',18)
        curvefits{subject,condition} = params;        
    end    
end
save('zROC_curvefits.mat','curvefits')

%% One-sample t-test on the curve parameter (a < 0) for each contrast
load zROC_curvefits
    for subject = subjects
        for condition = 1:3
            a(subject,condition) = curvefits{subject,condition}(1);
        end
    end
    
    %1-way rep measures ANOVA on p1
    [p,t,stats] = onewayRepmeasuresANOVA(a,'Contrast');
    %Save p and F values
    F = t{2,6}; P = p(1);
    
    %Means
    mean_p = mean(a)';    
    %One sample t-tests
    [h,p,~,tstats] = ttest(a)
    
    %Pairwise comparisons 
    [pval,tstat] = pairwiseComparisons(a)




