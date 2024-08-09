load("SNr_new.mat");


% N is the number of neurons
N=34;
% S is the number of stimuli conditions (F1 frequencies in Romo's task)
Num_dspsize=4;
Num_conditions=2; % TP/TA

S=2; %Total number of stimuli conditions Find TP/TA
S_second=4; % display size

% D is the number of decisions (D=2)
D=2; %Accept or reject

% T is the number of time-points (note that all the trials should have the
T=1600; % Number of bins

Trial_total=240; %240 trial for each session
Max_trial_condition=39;


% Extracting data
% Bin names
column_names = cell(1, 1600);

for i = 1:1600
    column_names{i} = ['bin', num2str(i)];
end
Subject_number=3;
Temp = find(table.Subject==Subject_number);
Subject_data=table(Temp,:);

%extracting firing rate
Subject_data_firing= Subject_data{:, column_names(1:1600)};



% Find the unique numbers of sessions
Subject_sessions=Subject_data{:,57};

unique_sessions = unique(Subject_sessions);
num_sessions = length(unique_sessions);

E=Max_trial_condition;


%% Preparing data for pca
% We have already PSTH data so there is no 
% need to average among data points
firingRatesAverage=Subject_data_firing;
Trial_num=13258;



% Remove any rows with any NaN values
firingRatesAverage_clean = firingRatesAverage(~any(isnan(firingRatesAverage), 2), :);

combinedParams = {{1, [1 3]}, {2, [2 3]}, {3}, {[1 2], [1 2 3]}};
margNames = {'Stimulus', 'Decision', 'Condition-independent', 'S/D Interaction'};
margColours = [23 100 171; 187 20 25; 150 150 150; 114 97 171]/256;


% firingRates array has [N S D T E] size; herewe ignore the 1st dimension 
% (neurons), i.e. we have the following parameters:
%    1 - stimulus 
%    2 - decision
%    3 - time
% There are three pairwise interactions:
%    [1 3] - stimulus/time interaction 
% I think we should use this pairwise 
% As explained in the eLife paper, we group stimulus with stimulus/time interaction etc.:

%% Step 1: PCA of the dataset

X = firingRatesAverage_clean(:,:);
X = bsxfun(@minus, X, mean(X,2));

[W,~,~] = svd(X, 'econ');
W = W(:,1:20);

% minimal plotting
dpca_plot(firingRatesAverage_clean, W, W, @dpca_plot_default);

%% compute explain variance
% computing explained variance
explVar = dpca_explainedVariance(firingRatesAverage_clean, W, W);
time_axis = linspace(-0.4, 1.0, size(firingRatesAverage_clean,2));
timeEvents=0;

% a bit more informative plotting
dpca_plot(firingRatesAverage_clean, W, W, @dpca_plot_default, ...
    'explainedVar', explVar, ...
        'time', time_axis,                        ...
    'timeEvents', timeEvents,               ...
    'marginalizationNames', margNames, ...
    'marginalizationColours', margColours);


%% Step 2: PCA in each marginalization separately

dpca_perMarginalization(firingRatesAverage_clean, @dpca_plot_default);



%% Step 3: dPCA without regularization and ignoring noise covariance

% This is the core function.
% W is the decoder, V is the encoder (ordered by explained variance),
% whichMarg is an array that tells you which component comes from which
% marginalization

tic
[W,V,whichMarg] = dpca(firingRatesAverage_clean, 20);
toc

explVar = dpca_explainedVariance(firingRatesAverage_clean, W, V);

dpca_plot(firingRatesAverage_clean, W, V, @dpca_plot_default, ...
    'explainedVar', explVar, ...
    'marginalizationNames', margNames, ...
    'marginalizationColours', margColours, ...
    'whichMarg', whichMarg,                 ...
    'time', time_axis,                        ...
    'timeEvents', timeEvents,               ...
    'timeMarginalization', 3, ...
    'legendSubplot', 16);

%% Step 4: dPCA with regularization

% This function takes some minutes to run. It will save the computations 
% in a .mat file with a given name. Once computed, you can simply load 
% lambdas out of this file:
%   load('tmp_optimalLambdas.mat', 'optimalLambda')

% Please note that this now includes noise covariance matrix Cnoise which
% tends to provide substantial regularization by itself (even with lambda set
% to zero).

%check false/true
optimalLambda = dpca_optimizeLambda(firingRatesAverage, firingRates, trialNum, ...
    'simultaneous', false , ...
    'numRep', 2, ...  % increase this number to ~10 for better accuracy
    'filename', 'tmp_optimalLambdas.mat');

Cnoise = dpca_getNoiseCovariance(firingRatesAverage, ...
    firingRates, trialNum, 'simultaneous', ifSimultaneousRecording);

[W,V,whichMarg] = dpca(firingRatesAverage, 20, ...
    'combinedParams', combinedParams, ...
    'lambda', optimalLambda, ...
    'Cnoise', Cnoise);

explVar = dpca_explainedVariance(firingRatesAverage, W, V, ...
    'combinedParams', combinedParams);

dpca_plot(firingRatesAverage, W, V, @dpca_plot_default, ...
    'explainedVar', explVar, ...
    'marginalizationNames', margNames, ...
    'marginalizationColours', margColours, ...
    'whichMarg', whichMarg,                 ...
    'time', time,                        ...
    'timeEvents', timeEvents,               ...
    'timeMarginalization', 3,           ...
    'legendSubplot', 16);











