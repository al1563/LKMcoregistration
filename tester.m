load data
addpath([pwd,'/SR']);
set(0,'DefaultFigureWindowStyle','docked'); % Figure is docked

options.N = 500;
options.nbin = 6; %6, 10
options.hsize = 250; %250, 150
%options.trainingindex = [1,2,6,9,12,15,16,17]; % Index of training data
options.trainingindex = 1:18;
options.alpha = .0001;
options.t = 4;

% % Train machine learning model
mlModel = LKM.trainModel(data, options);

avgtrainingparam = 0;
for j = options.trainingindex
   avgtrainingparam = avgtrainingparam + data{j}.gtparam;
end
avgtrainingparam = avgtrainingparam/5;
avgtrainingparam(4:5) = [0, 0];

%%
result = [];

for i = 15 % Data to test registration
    figure;
    
    % Test little difference from the groundtruth:
    initialparam = data{i}.gtparam + [.05 .05 .05 25 25 .2];%+ randn(1,6) .* [.1 .1 .3 5 5 .2]; 

    % How sensitive each parameter is 
    paramweight = [1 1 1 1 1 1];
    paramshift = [1 1 10 1000 1000 1];

    % % perform registration
    [T, sim] = LKM.register(data{i}.data3D, data{i}.data2D, options,...
                   mlModel,initialparam, 1, paramweight, paramshift);
               
    result = [result; data{i}.gtparam T];
end
