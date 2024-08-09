%PCA transform on data
% Assuming firing_rates is a matrix of size (n_trials x 1600) containing firing rates for each trial
% Assuming labels is a cell array of size (n_trials x 1) containing labels for each trial (TP or TA)
%% Defining full matrix


% Define the column names array
load("SNr_new.mat")


NeuronNumber=20;
Temp = find(table.iUnit==NeuronNumber);
Neuron = table(Temp,:);


TP = find(Neuron.EventValue==3);
TA = find(Neuron.EventValue==4);

%filename = "data.mat";
%save(filename,"Neuron")
column_names = cell(1, 1600);

for i = 1:1600
    column_names{i} = ['bin', num2str(i)];
end
%% 

Eff=find(table.Search_Type==1);
Ineff=find(table.Search_Type==0);
DspSize3 = find(Neuron.DispSize==7);
Temp=intersect(DspSize3,Eff);
Neuron_Dspsz3=Neuron(DspSize3,:);




Neuron_Dspsz3_bins = Neuron_Dspsz3{:, column_names(1:1600)};




Neuron_TPTA_labels =Neuron_Dspsz3{:,"EventValue"};
%% Apply PCA

% Load your data
% dataMatrix = ...; % Your 129x1600 matrix
% labels = ...; % Your 129x1 label array (3 for TP, 4 for TA)
labels=Neuron_TPTA_labels;
% Apply PCA
[coeff, score, ~, ~, explained] = pca(Neuron_Dspsz3_bins);

% Transform the data to 3 principal components
transformedData = score(:, 1:3);

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



