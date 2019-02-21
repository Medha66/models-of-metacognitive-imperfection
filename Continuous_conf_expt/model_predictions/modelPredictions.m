% Generate model predictions - metacognition measures and zROC functions

clear
clc

global modelToFit
global criterion
global condition

% Select model to fit
modelToFit = 'lognormal'; %lognormal, normal, standard

%Add results folder and helper functions
addpath(genpath('helperFunctions'),'fittingResults','analyses')

% Load the data
load DataForModeling
N = 10000;
nratings = 2;
nCriteria = 49;

%Load metaNoise estimates
if ~strcmp(modelToFit,'standard')
    load(['fittingResults/fittingResults_' modelToFit])
end

%Flags
save_predictions = 0;
predict = 1;

subjects = [1:20];


if predict
    for subject = subjects
        
        %Compute the fixed parameters
        sigma = 1;
        c = exptData.criteria(subject,nCriteria);
        mu = exptData.dprime(subject,:);
        
        for condition = 1:3
            for criterion = 1:nCriteria
                
                %Estimated parameters - conf criteria and metaNoise
                if strcmp(modelToFit,'standard')
                    metaNoise = 0;
                    [confCriteria,respProb] = estimateCriteria2([metaNoise],subject);
                else
                    metaNoise =  modelFit.params(subject);
                    [confCriteria,respProb] = estimateCriteria2([metaNoise],subject);
                end
                
                c_negative = confCriteria(1); %Left criterion
                c_positive = confCriteria(3); %Right criterion
                
                %Generate confidence criteria
                if strcmp(modelToFit,'lognormal')
                    confCriteria_negative_ = lognrnd(c_negative, metaNoise, [1,N]);
                    confCriteria_positive_ = lognrnd(c_negative, metaNoise, [1,N]);
                else
                    confCriteria_negative_ = normrnd(c_negative, metaNoise,[1,N]);
                    confCriteria_positive_ = normrnd(c_positive, metaNoise,[1,N]);
                end
                
                % Simulate model activations
                %Generate actual stimulus values - 0 or 1
                stimulus = [zeros(1,N/2), ones(1,N/2)];
                
                %Sensory response
                for stimID = 0:1
                    rsens(1+N/2*(stimID):N/2*(stimID+1)) = normrnd((2*stimID-1)*mu(condition)/2, sigma, 1, N/2); %r_sens
                end
                
                %Decision responses
                response = rsens > c; % 1 - right & 0 - left
                
                %Correct responses
                correct = response == stimulus;
                
                %Confidence responses
                confCriteria_negative = [c*ones([1,N]);-confCriteria_negative_+c;-Inf*ones([1,N])];
                confCriteria_negative(confCriteria_negative>c)=c;
                confCriteria_positive = [c*ones([1,N]);confCriteria_positive_+c;Inf*ones([1,N])];
                confCriteria_positive(confCriteria_positive<c)=c;
                
                confidence = zeros([1,N]);
                for confValue = 1:2
                    confidence(rsens>=c & rsens>=confCriteria_positive(confValue,:)& rsens<confCriteria_positive(confValue+1,:))=confValue;
                    confidence(rsens<c & rsens<confCriteria_negative(confValue,:)& rsens>=confCriteria_negative(confValue+1,:))=confValue;
                end
                
                %Predicted Metacognition Measures
                %phi
                R = corrcoef(correct,confidence);
                phi(subject,condition,criterion) = R(2);
                %Type-2 AUC
                [nR_S1, nR_S2] = trials2counts(stimulus,response,confidence,nratings,[]);
                [t2AUC(subject,condition,criterion)] = type2ag(nR_S1, nR_S2, 1);
                %meta-d' and mratio
                meta_MLE = type2_SDT_MLE(stimulus,response,confidence,2,[],1);
                mratio(subject,condition,criterion) = meta_MLE.M_ratio;
                meta_d(subject,condition,criterion) = meta_MLE.meta_da;
                
                % Mean confidence
                meanConf(subject,condition,criterion) = mean(confidence);
                
                %Predicted ROC functions
                %From confcriteria
                [FAR_,HR_] = zROC(stimulus,response,confidence,2,1);
                FAR(subject,nCriteria+1-criterion,condition) = FAR_(1); FAR(subject,nCriteria+1+criterion,condition) = FAR_(3);
                HR(subject,nCriteria+1-criterion,condition) = HR_(1); HR(subject,nCriteria+1+criterion,condition) = HR_(3);
                
            end
            %HR and FAR from decision criterion
            HR(subject,nCriteria+1,condition) = length(find(correct==1 & stimulus==1))./length(find(stimulus==1));
            FAR(subject,nCriteria+1,condition) = length(find(correct==0 & stimulus==0))./length(find(stimulus==0));
        end
    end
else
    load(['predicted_measures_',modelToFit,'.mat']);
    load(['zROC_predicted_',modelToFit,'.mat']);
    
end

%Save data
if save_predictions
    save(['analyses/predicted_measures_',modelToFit,'.mat'],'phi','t2AUC','mratio','meta_d');
    save(['analyses/zROC_predicted_',modelToFit,'.mat'],'HR','FAR');
end


%% Figure 6 - plot predicted metacognition measures against confidence criterion location for all the models


axes('box','off','tickdir','out','LineWidth',1.25,'FontSize',16); hold on;

%Load data
models = {'standard','normal','lognormal'};
model_names = {'Standard SDT','Gaussian meta noise','Lognormal meta noise'};

%Line colors for each contrast
colors = {[.8 .1 .1],[.4 .6 .4],[.8 .8 1]};

%y-limits for each measure
ylimits = {[.5,1.1],[.5,3],[.5,.75],[.05,.35]};


for model = 1:3
    load(['predicted_measures_',models{model},'.mat'])
    
    %Exclude subject 7 (response probabilities of higher conf criteria are
    %too low to generate predictions)
    phi(7,:,:) = NaN;
    t2AUC(7,:,:) = NaN;
    mratio(7,:,:) = NaN;
    meta_d(7,:,:) = NaN;
    
    criterion_vect = [51:1:99];
    subjects = length(mratio);
    
    measures = {'mratio','meta_d','t2AUC','phi'};
    names = {'meta-d''/d''','meta-d''','Type-2 AUC','phi'};
    
    numMeasures = length(names);
    
    for measure = 1:numMeasures
        y_ = eval(measures{measure});
        
        %Plot data
        subplot(3,4, (model-1)*numMeasures+measure); hold on;
        for condition = 1:3
            for criterion = 1:length(criterion_vect)
                y(criterion) = nanmean(y_(:,condition,criterion));
                y_sem(criterion) = nanstd(y_(:,condition,criterion))./sqrt(20);
            end
            %Mratio for standard model as dotted lines
            if model == 1 && measure == 1 && condition == 1
                h(condition) = plot(criterion_vect,y,'--','linewidth',2,'color',colors{condition}); hold on;
                shadedplot(criterion_vect,y-y_sem,y+y_sem,colors{condition})
            else
                h(condition) = plot(criterion_vect,y,'linewidth',2,'color',colors{condition}); hold on;
                shadedplot(criterion_vect,y-y_sem,y+y_sem,colors{condition})
            end
            
        end
        
        %x and y-limits
        xlim([50 100])
        ylim(ylimits{measure})
        
        %x and y-labels
        ylabel(names{measure},'fontsize',18)
        if model == 3
            xlabel('Confidence criterion location','fontsize',14)
        end
    end
end

l = legend([h],'Contrast 1','Contrast 2','Contrast 3');
set(l,'box','off','orientation','horizontal','fontsize',18)

%% Analyses

% Curve fitting
fig_count = 0;

for model = 1:3
    load(['predicted_measures_',models{model},'.mat'])
    for subject = subjects
        for condition = 1:3
            figure;
            fig_count = fig_count+1;
            for measure = 1:numMeasures
                y_ = eval(measures{measure});
                
                y = y_(subject,condition,:);
                if measure < 3
                    clear p p1 p2
                    subplot(2,2,measure); hold on;
                    [fitresult, gof] = PolynomialFit(criterion_vect,y,1) %Linear fits
                    %Save the slope and intercept values
                    a_model = fitresult.p1; %Slope
                    b = fitresult.p2; %Intercept
                    params = [a_model,b];
                else
                    clear p p1 p2
                    subplot(2,2,measure)
                    [fitresult, gof] = PolynomialFit(criterion_vect,y,2) %Quadratic fits
                    %Save the polynomial coefficients
                    a_model = fitresult.p1; %degree of curvature - p1 < 0 --> curves downwards
                    b = fitresult.p2; %location
                    c = fitresult.p3; %y-intercept
                    params = [a_model,b,c];
                end
                
                ylabel(names{measure},'fontsize',18)
                xlabel('Confidence criterion','fontsize',18)
                curvefits{subject,measure,condition} = params;
            end
            %saveFigure(num2str(fig_count));
            %delete gcf
        end
        
    end
    % Save fitting results
    save(['analyses/measures_curvefits_',models{model},'.mat'],'curvefits')
end

%% One sample t-test for a < 0

%Models
models = {'standard','normal','lognormal'};
model_names = {'Standard SDT','Gaussian meta noise','Lognormal meta noise'};

for model = 1:3
    load(['measures_curvefits_',models{model},'.mat'])
    for measure = 1:4
        for condition = 1:3
            for subject = subjects                
                a_model(subject,condition) = curvefits{subject,measure,condition}(1);
            end
        end
        %one-sample t-tests
        [h,p,~,tstats] = ttest(a_model);
        mean_a(measure,:,model) = mean(a_model);
        pval(measure,:,model) = p;
        tstat(measure,:,model) = tstats.tstat;
                
    end
end

%% Plot the predicted zROC function

%Select indices of the 6 criteria used for model fitting
criterion_vect = [0 .51:.01:.99 Inf];
conf_binSize = .5/6;
for rating = 1:5
    model_confCriteria(rating) = .5+(rating)*conf_binSize;
    modelCriteria_index(rating) = find(criterion_vect<model_confCriteria(rating), 1, 'last' );
end

axes('box','off','tickdir','out','LineWidth',1.25,'FontSize',16); hold on;
xlims = {[-4,2],[-4,2],[-3,2]};
ylims = {[-2,3],[-2,3],[-2,3]};

for model = 1:length(models)
    
    subplot(2,2,model)
    load(['zROC_predicted_',models{model},'.mat'])
    
    for condition = 1:3
        
        %zHR and zFAR
        zHR = norminv(HR);
        zFAR = norminv(FAR);
        
        %Average across subjects
        zHR_avg = nanmean(zHR(:,:,condition));
        zFAR_avg = nanmean(zFAR(:,:,condition));
        nPoints = length(zHR_avg);
        
        h1(condition) = plot(zFAR_avg, zHR_avg, '.','markersize',30,'color',colors{condition});hold on;
        %Plot decision criterion in a different color
        h2 = plot(zFAR_avg(nCriteria+1), zHR_avg(nCriteria+1),'*k','markersize',10,'linewidth',2);
        %Plot the criteria used for modelfitting in a separate color
        h3 = plot(zFAR_avg([modelCriteria_index, nCriteria+1 + modelCriteria_index]), zHR_avg([modelCriteria_index, nCriteria+1 + modelCriteria_index]),'x','color',[0 0 0],'markersize',7.5,'linewidth',2);
        
       
        %Plot the linear zROC plot implied by the decision criterion
        dprime = zHR_avg(nCriteria+1)-zFAR_avg(nCriteria+1);
        xlimits = [-4, 2.5];
        plot(xlimits, dprime+[xlimits], '--','color',[.5 .5 .5],'linewidth',1.5)
        
        xlim(xlims{condition});
        ylim(ylims{condition});
        xlabel({'zFAR';''},'fontsize',18);
        ylabel('zHR','fontsize',18);
        
    end
    if model == 1
        l = legend([h2,h3,h1],'Decision criterion','Confidence criteria used for model fitting','Confidence criteria - Contrast 1','Confidence criteria - Contrast 2','Confidence criteria - Contrast 3','location','best');
    end
    set(l,'fontsize',18)
    title([model_names{model}],'fontsize',18)
    %set(h1(2),'color',
end

%% Do curve fitting for each subject (after rotating the ROC plot by 45 degrees)

for model = 1:3
    load(['zROC_predicted_',models{model},'.mat'])
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
            a_model = fitresult.p1; %degree of curvature - p1 < 0 --> curves downwards
            b = fitresult.p2; %location
            c = fitresult.p3; %y-intercept
            params = [a_model,b,c];
            
            ylabel('zHR','fontsize',18)
            xlabel('zFAR','fontsize',18)
            curvefits{subject,condition} = p;
        end
    end
    save(['zROC_curvefits_',models{model},'.mat'],'curvefits')
end

%% Analysis of the curve parameter (a < 0) 
clear mean_a tstat pval
for model = 1:3
    
    load(['zROC_curvefits_',models{model},'.mat'])
    for subject = subjects
        for condition = 1:3
            a_model(subject,condition) = curvefits{subject,condition}(1);
        end
    end


%1-way rep measures ANOVA on a with stimulus contrast as factor
[p,t,stats] = onewayRepmeasuresANOVA(a_model,'Contrast');
%Save p and F values
F = t{2,6}; P = p(1);

%Pairwise t-tests
[pw_pval(:,model),pw_tstat(:,model)] = pairwiseComparisons(a_model)

%Means
mean_a(:,model) = mean(a_model);

%One sample t-tests
[h,p,~,tstats] = ttest(a_model)
os_tstat(:,model) = tstats.tstat;
os_pval(:,model) = p';

% Comparison of predicted curvature  with empirical estimates 
   %Data
    load('zROC_curvefits')
    for subject = subjects
        for condition = 1:3
            a_data(subject,condition) = curvefits{subject,condition}(1);
        end
    end
    
    %Paired t-tests
    for condition = 1:3
        [h,pw2_p,~,pw2_tstats] = ttest(a_data(:,condition),a_model(:,condition));
        pw2_tstat(condition,model) = pw2_tstats.tstat';
        pw2_pval(condition,model) = pw2_p';
        mean_difference(condition,model) = mean(a_data(:,condition)-a_model(:,condition))
    end

end



