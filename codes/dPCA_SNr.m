%Applying dPCA to my data
%This is for SNr data

%% Preparing data

load("SNr_new.mat")

%We started with monkey P
subject3=table(1:13259,:);
Unique_neurons=unique(subject3.iUnit);


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




Main_data=zeros(28,2,1,1600);




for i=1:(length(Unique_neurons)-1)

    NeuronNumber=Unique_neurons(i);
    Temp = find(table.iUnit==NeuronNumber);
    Neuron = table(Temp,:);
    TPOnly = find(Neuron.EventValue==3);
    TAOnly = find(Neuron.EventValue==4);
    Neuron_TPOnly=Neuron(TPOnly,:);
    Neuron_TAOnly=Neuron(TAOnly,:);
    Neuron_TPOnly_bins = Neuron_TPOnly{:, column_names(1:1600)};
    Neuron_TAOnly_bins = Neuron_TAOnly{:, column_names(1:1600)};
    Neuron_TPOnly_bins_mean=nanmean(Neuron_TPOnly_bins,1);
    Neuron_TAOnly_bins_mean=nanmean(Neuron_TAOnly_bins,1);
    %TP_TA_merged=[Neuron_TPOnly_bins_mean, Neuron_TAOnly_bins_mean];
    Main_data(i,1,1,:)=Neuron_TPOnly_bins_mean;
    Main_data(i,2,1,:)=Neuron_TAOnly_bins_mean;

    %Main_data(i,:)=TP_TA_merged;


end    


%% Applying method

% firingRatesAverage: N x S x D x T
firingRatesAverage=Main_data(:,:,:,:);
N=28; % Usefull Number of neurons
S=2; % Number of conditions
D=1; % Number of Decisions
T=1600; % Number of time steps (Bins)



% firingRates array has [N S D T E] size; herewe ignore the 1st dimension 
% (neurons), i.e. we have the following parameters:
%    1 - stimulus 
%    2 - decision
%    3 - time
% There are three pairwise interactions:
%    [1 3] - stimulus/time interaction
%    [2 3] - decision/time interaction
%    [1 2] - stimulus/decision interaction
% And one three-way interaction:
%    [1 2 3] - rest
% As explained in the eLife paper, we group stimulus with stimulus/time interaction etc.:

combinedParams = {{1, [1 3]}, {2, [2 3]}, {3}, {[1 2], [1 2 3]}};
margNames = {'Stimulus', 'Decision', 'Condition-independent', 'S/D Interaction'};
margColours = [23 100 171; 187 20 25; 150 150 150; 114 97 171]/256;

% Time events of interest (e.g. stimulus onset/offset, cues etc.)
timeEvents = 0;
time_axis = linspace(-0.4, 1.0, size(firingRatesAverage,2));

%% Step 1: PCA of the dataset

X = firingRatesAverage(:,:);
X = bsxfun(@minus, X, mean(X,2));

[W,~,~] = svd(X, 'econ');
W = W(:,1:20);


















% Z=X'*W;
% 
% Z=Z';
% 
% Compnent1_TP=Z(1,1:1600);
% Compnent1_TA=Z(1,1601:3200);
% 
% 
% Compnent2_TP=Z(1,1:1600);
% Compnent2_TA=Z(1,1601:3200);
% 
% plot(1:1600,Compnent1_TA')
% hold on
% plot(1:1600,Compnent1_TP)

% minimal plotting
dpca_plot(firingRatesAverage, W, W, @dpca_plot_default);

%% compute explain variance
% computing explained variance
explVar = dpca_explainedVariance(firingRatesAverage, W, W);


% a bit more informative plotting
dpca_plot(firingRatesAverage, W, W, @dpca_plot_default, ...
    'explainedVar', explVar, ...
        'time', time_axis,                        ...
    'timeEvents', timeEvents,               ...
    'marginalizationNames', margNames, ...
    'marginalizationColours', margColours);

%% 3 components with Marginalization

dpca_perMarginalization(firingRatesAverage, @dpca_plot_default);

%% Step 3: dPCA without regularization and ignoring noise covariance

% This is the core function.
% W is the decoder, V is the encoder (ordered by explained variance),
% whichMarg is an array that tells you which component comes from which
% marginalization

tic
[W,V,whichMarg] = dpca(firingRatesAverage, 20);
toc



%% Test

Xt=firingRatesAverage(:,:);
XfullCen = reshape(Xt, size(firingRatesAverage));











%% 

explVar = dpca_explainedVariance(firingRatesAverage, W, V);

dpca_plot(firingRatesAverage, W, V, @dpca_plot_default, ...
    'explainedVar', explVar, ...
    'marginalizationNames', margNames, ...
    'marginalizationColours', margColours, ...
    'whichMarg', whichMarg,                 ...
    'time', time_axis,                        ...
    'timeEvents', timeEvents,               ...
    'timeMarginalization', 3, ...
    'legendSubplot', 16);



%% Apply Simple PCA

% Load your data
% dataMatrix = ...; % Your 129x1600 matrix
% labels = ...; % Your 129x1 label array (3 for TP, 4 for TA)
% Apply PCA
[coeff, score, ~, ~, explained] = pca(firingRatesAverage);

% Transform the data to 3 principal components
transformedData = score(:, 1:3);


hold on
bar(explained)
plot(1:numel(explained), cumsum(explained), 'o-', 'MarkerFaceColor', 'r')
yyaxis right
h = gca;
h.YAxis(2).Limits = [0 100];
h.YAxis(2).Color = h.YAxis(1).Color;
h.YAxis(2).TickLabel = strcat(h.YAxis(2).TickLabel, '%');

%% 

 
cumexplained = cumsum(explained);
cumunexplained = 100 - cumexplained;
plot(1:27, cumunexplained, 'x-');
grid on;
xlabel('Number of factors');
ylabel('Unexplained variance')

%% 


% Create 3D scatter plot
figure;
hold on;
for i = 1:size(transformedData, 1)
    if labels(i) == 3
        % TP class
        scatter3(transformedData(i, 1), transformedData(i, 2), transformedData(i, 3), 'r');
    elseif labels(i) == 4
        % TA class
        scatter3(transformedData(i, 1), transformedData(i, 2), transformedData(i, 3), 'b');
    end
end
xlabel('1st Principal Component');
ylabel('2nd Principal Component');
zlabel('3rd Principal Component');
title('3D PCA Plot');
grid on;
hold off;

view(-30,10)

% Change the view to the angle that provides the best separation
%view([azimuth,elevation]);

% Optionally, you can add rotation controls to interact with the plot
rotate3d on;


%% PCA with 2 components


% Transform the data to 2 principal components
transformedData2 = score(:, 1:2);
figure;
hold on;
for i = 1:size(transformedData2, 1)
    if labels(i) == 3
        % TP class
        scatter(transformedData2(i, 1), transformedData2(i, 2), 'r');
    elseif labels(i) == 4
        % TA class
        scatter(transformedData2(i, 1), transformedData2(i, 2), 'b');
    end
end
xlabel('1st Principal Component');
ylabel('2nd Principal Component');
title('2D PCA Plot');
grid on;
hold off;



%% 
hold on
bar(explained)
plot(1:numel(explained), cumsum(explained), 'o-', 'MarkerFaceColor', 'r')
yyaxis right
h = gca;
h.YAxis(2).Limits = [0 100];
h.YAxis(2).Color = h.YAxis(1).Color;
h.YAxis(2).TickLabel = strcat(h.YAxis(2).TickLabel, '%');


%% 
cumexplained = cumsum(explained);
cumunexplained = 100 - cumexplained;
plot(1:110, cumunexplained, 'x-');
grid on;
xlabel('Number of factors');
ylabel('Unexplained variance')



