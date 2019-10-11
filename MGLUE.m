%%
% Code to implement generalized likelihood uncertainty estimation (GLUE)
% Developed by: Ebrahim Ahmadisharaf
% Last updated: 5 Oct 2019
%%
tic;
clear all
clc;
%% Code settings
N = 2; % Shaping factor
BehavSelect = 1; % Approach to select behavioral simulations: 1: Cutoff threshold; 2: Top percentage
Behavtoppct = 0.1; % Percent of top simulations that are selected behabioral; Only needed when BehavSelect = 2

%% Read ensemble paramters
FileName = 'MGLUE_Inputs.xlsx'; % Spreadsheet containing inputs

%% Read inputs
% Read likelihood functions
[~,LFnam,~] = xlsread(FileName,'LFnam'); % Names of the likelihood functions
LFval = xlsread(FileName,'LFval'); % Values of the likelihood functions
LFcutoff = xlsread(FileName,'LFcutoff'); % Cutoff threshold of the likelihood functions; ; Only needed when BehavSelect = 1

% Read parameters
[num,parnam,raw] = xlsread(FileName,'parnam'); % Names of the parameters
parval = xlsread(FileName,'parval'); % Values of the parameters
parpriorp = xlsread(FileName,'parpriorp'); % Prior probabilities of the parameters

%% Get number of simulations and parameters
n = length(parval); % Number of burn-in simulations
for i=1:n
    simid(i) = i; % Vector showing the simulation IDs
end
npar = length(parnam);

%% Calculate cumulative prior probabilities
for i=1:length(parnam)
    x(:,1) = parval(:,i);
    x(:,2) = parpriorp(:,i);
    sortpriorp = sortrows(x); % Sorting parameter ascendingly based on parameter values
    priorpar(:,(i-1) * 2 + 1) = sortpriorp(:,1); % Sorted paramater values in all the simulations
    priorpar(:,(i-1) * 2 + 2) = cumsum(sortpriorp(:,2)); % Prior cumulative probability of parameters
    clear x
end

%% Select behavioral simulations
disp('It is assumed that the greater the value of likelihood function, the better the model performnace');
ilf = input('Enter the index of the likelihood function (in the order of the input spreadsheet worksheet LFnam): '); % Index of the primary likelihood function
if BehavSelect == 1
    j = 1; % Index that counts the behavioral simiulations
    for i=1:n
        if LFval(i,ilf) >= LFcutoff(ilf)
            LFBehavioral(j,1) = LFval(i,ilf);
            LFBehavioral(j,2) = i; % Index of behavioral simulation
            j = j+1;
        end
    end
else
    x(:,1) = LFval(:,ilf); % Values of likelihoodd function
    x(:,2) = transpose(simid); % Simulation IDs
    sortLF = sortrows(x); % Sorting simulations ascendingly based on the likelihood function
    LFBehavioral = sortLF(round(n-Behavtoppct*n+1):n,:);
end
clear x

disp(['There are ',num2str(length(LFBehavioral)),' behavioral simulations']);
%% Calculate likelihood and posterior probability
for j=1:length(LFBehavioral)
    L(j) = LFBehavioral(j,1).^N; % Likelihood
end
normL(:) = L/sum(L(:)); % Normalized likelihood
for i=1:length(parnam)
    for j=1:length(LFBehavioral)
        x(j) = normL(j)./parpriorp(LFBehavioral(j,2),i);
    end
    parpostp(:,i) = x/sum(x(:)); % Posterior probability
    clear x
end

%% Cumulative posterior probability calculation
for i=1:length(parnam)
    for j=1:length(LFBehavioral)
        x(j,1) = parval(LFBehavioral(j,2),i);
    end
    x(:,2) = parpostp(:,i);
    sortpostp = sortrows(x); % Sorting behavioral parameters ascendingly based on parameter values
    postcump(:,(i-1) * 2 + 1) = sortpostp(:,1); % Sorted paramater values in the behavioral simulations
    postcump(:,(i-1) * 2 + 2) = cumsum(sortpostp(:,2)); % Posterior cumulative probability of parameters
    clear x
end

%% Plotting posterior distributions of the parameters
for i=1:length(parnam)
    figure(i);
    plotpost = plot(postcump(:,(i-1) * 2 + 1),postcump(:,(i - 1) * 2 + 2),'LineWidth',4); % Posterior probabilities
    hold on
    plotpost = plot(priorpar(:,(i-1)* 2 + 1),priorpar(:,(i - 1) * 2 + 2),'LineWidth',4); % Prior probabilities
    xlabel(parnam(i),'fontweight','bold','fontsize',18);
    ylabel('Cumulative probability','fontweight','bold','fontsize',18);
%     legend({'Posterior distribution','Prior distribution'},'Location','southeast');
    ylim([0 1]); % Probability ranges from 0 to 1
    set(gca, 'FontName', 'Times New Roman');
    set(gca,'fontsize',18);
    set(gca,'fontweight','bold');
    set(gca,'LineWidth',2);
    parfig = [pwd,filesep,str2mat(parnam(i)),'.tif'];
    saveas(gca,parfig);
end

%%
toc;