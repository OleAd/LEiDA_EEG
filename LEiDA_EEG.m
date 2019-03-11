%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LEiDA_EEG
% This script demonstrates how the LEiDA-method can be applied to EEG-data.
% It is based on Joana Cabral's work, available at
% github.com/juanitacabral/LEiDA
% 
% First implementation
% Joana Cabral Oct 2017
% First use in
% Cabral, et al. 2017 Scientific Reports 7, no. 1(2017): 5135.
% 
% Adapted by Ole Adrian Heggli
% ole.heggli@clin.au.dk
% 
% REQUIREMENTS:
% This script needs the SPM functions spm_vol() and spm_slice_vol(), hence 
% SPM needs to be downloaded and added to path.
% This script also needs the circ_mean() function for the Circular
% Statistic Toolbox (Directional Statistics), either in same folder or on
% path.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Load data
% This example script assumes you have source reconstructed EEG-data, based
% on the parcellation used in by Colclough et. al. in
% "Colclough, Giles L., et al. "A symmetric multivariate leakage correction 
% for MEG connectomes." Neuroimage 117 (2015): 439-448."
% If your data uses a different parcellations, specify this is N_areas, and
% rework the LEiDA_EEG_plot_nodes function.
% 
% The EEG-data should either be collated or continous.
% A sample EEG-dataset is included here.
% This implementation of LEiDA compares occurence probabilities between two
% groups.
% 
% This sample dataset contains "eight" subjects, four in each group.
% The data is scrambled anonymous data, and is provided only for
% informational purposes.
% For an actual dataset, it may be necessary to run the calcaulations using HPC. 


load('exampleSourceEEG_8.mat')
% put the data in a cell array
exampleData={p1, p2, p3, p4, p5, p6, p7, p8};
clear p1 p2 p3 p4 p5 p6 p7 p8
% Number of subjects
numSubjects=8;
% Create a vector describing which groups the participants belong to
group=[1,2,1,2,1,2,1,2];
% 



%% Step one, calculating eigenvectors

FS=250;             % Samplerate of the EEG signal
N_areas=39;         % Brain regions in the data
window_size=250;    % Number of samples to calculate the circular mean over
frequency='alpha';  % Select frequency band.

% Run the eigenvector calculation
eigenvectors=[];
design=[];
for s=1:size(exampleData,2)
    [thisEigen, thisDesign]=LEiDA_EEG_eigenvectors(exampleData{s}, group(s),...
        FS, N_areas, window_size, frequency);
    eigenvectors{s}=thisEigen;
    design{s}=thisDesign;
end

% 

%% Step two, k-means clustering

% First, collect all eigenvectors
collEigenvectors=[];
for s=1:numSubjects
    collEigenvectors=[collEigenvectors;eigenvectors{s}];
end
clear eigenvectors

% Select the range of k you want to cluster for. Normally anywhere between
% 2 and 20.
minK=2;
maxK=6;
rangeK=minK:maxK;

gpu=1; % Set gpu to 1 if you want to use gpu-acceleration, 0 if not.
% Beware that if you do a parfor with gpu-acceleration, you'll likely run
% into memory issues.

% Do K-means clustering
kMeans_results={};
for k=1:length(rangeK)
    thisK=rangeK(k);
    thisKmeans=LEiDA_EEG_kmeans(collEigenvectors, thisK, gpu);
    kMeans_results{k}=thisKmeans;
end



% 


%% Step three, selecting the best k
% First calculate probabilities
P=zeros(numSubjects,maxK-minK+1,maxK);

subjectIndicator=[];
for s=1:numSubjects
    currDesign=design{s};
    subjectIndicator=[subjectIndicator, repelem(s, length(currDesign))];
end

% Get the results per subject.
for k=1:length(rangeK)
        for s=1:numSubjects
            thisDesign=design{s};
            % Select the time points representing this subject              
            T=subjectIndicator==s;
            % This finds the time belonging to the subject
            Ctime=kMeans_results{k}.IDX(T);
            for c=1:rangeK(k)
                P(s,k,c)=mean(Ctime==c);
            end
    end
end

% Do statistics
% First define which participants belong to which group
Group1subj=[1,3,5,7];
Group2subj=[2,4,6,8];

disp('Test significance between groups')
P_pval = zeros(maxK-minK+1,maxK);
pvals=[];
for k=1:length(rangeK)
    disp(['Now running statistics for ' num2str(rangeK(k)) ' clusters'])
    for c=1:rangeK(k)
        % Compare Probabilities
        a=squeeze(P(Group1subj,k,c))';  % Vector containing prob of c in Group 1
        b=squeeze(P(Group2subj,k,c))';  % Vector containing prob of c in Group 2
        if ismember(a,0)
            disp('Found zero possibility.')
%             This means that one of the groups has zero possibility of
%             being in this state. That is possible, but not plausible,
%             depending on your data.
        end
        if ismember(b,0)
            disp('Found zero possibility.')
        end
        if sum([a,b])==0
           disp('Found a non-probable state')
           warning('Proceed with caution, probably an error in your design')
           pvals(c,:)=[NaN,NaN];
        else
           stats=permutation_htest2_np([a,b],[ones(1,numel(a)) 2*ones(1,numel(b))],1000,0.05,'ttest2');
           pvals(c,:)=stats.pvals;
        end
        P_pval(k,c)=min(stats.pvals);  
    end
   
end


% Plot and inspect to find k
P_Prob=P_pval;
sigC=zeros(1,length(rangeK));
for k=1:length(rangeK)   
    [Min_p_value(k), sigC(k)]=min(P_Prob(k,P_Prob(k,:)>0));    
end
figure

semilogy(rangeK(1)-1:rangeK(end)+1,0.05./(rangeK(1)-1:rangeK(end)+1).*ones(1,length(rangeK)+2),'g--','LineWidth',2)
hold on

for k=1:length(rangeK) 
    for c=1:rangeK(k)
        if P_Prob(k,c)<=0.05 && P_Prob(k,c)>(0.05/rangeK(k))
            semilogy(rangeK(k),P_Prob(k,c),'+r');
        end
        if P_Prob(k,c)>=0.05
            semilogy(rangeK(k),P_Prob(k,c),'.k','MarkerSize',10);
        end
        if P_Prob(k,c)<=(0.05/rangeK(k)) && P_Prob(k,c)>(0.05/sum(rangeK))
            semilogy(rangeK(k),P_Prob(k,c),'*g','MarkerSize',8);
        end
        if P_Prob(k,c)<=(0.05/sum(rangeK))
            semilogy(rangeK(k),P_Prob(k,c),'*b');
        end
    end
end
ylabel('p-value')
xlabel('K')
set(gca,'XTick', rangeK(1:2:end))
ylim([1e-4 1])
xlim([rangeK(1)-1 rangeK(end)+1])
title('P-values')
box off

% Choose your K based on visual inspection. Green line is corrected
% threshold.



% 
%% Select K
K=6; % Input your choosen K here.
Best_Clusters=kMeans_results{rangeK==K};
k=find(rangeK==K);

thisPvalues=P_Prob(k,1:K);
thisSignif=find(thisPvalues<0.05/K);

% Get the K patterns
V=Best_Clusters.C;
[~, N]=size(Best_Clusters.C);


% 

%% Step four, plotting data

% This plots the nodes and connections on the cortex, as well as the
% occurence probabilities.
figure;
colormap(jet)
counter=1;
columns=ceil(K/2);
rows=4;
for c=1:K
    subplot(columns,rows,counter) 
    counter=counter+1;
    Vc=V(c,:);
    LEiDA_EEG_plot_nodes(Vc)
    subplot(columns,rows,counter)
    counter=counter+1;
    Group1=squeeze(P(Group1subj,k,c));
    Group2=squeeze(P(Group2subj,k,c));
    bar([mean(Group1) mean(Group2)], 'EdgeColor','w','FaceColor',[.5 .5 .5])
    hold on
    title({['State #' num2str(c)]})
    errorbar([mean(Group1) mean(Group2)],[std(Group1)/sqrt(numel(Group1)) std(Group2)/sqrt(numel(Group2))],'LineStyle','none','Color','k')
    set(gca,'XTick',1:2,'XTickLabel',{'Group 1','Group 2'})
    ylim([0, max(max(mean(Group1), mean(Group2)))+max(std(Group1), std(Group2))*2]);
    box off
    
end



% 