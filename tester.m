load data
addpath([pwd,'/SR']);

nbin = 8; 
hsize = 100;
N = 600;
% let N = 600 for about 1000 rows of training kernel per set

% Train machine learning model
mlModel = LKM.trainModel(data, N, nbin, hsize, trainindex);

% Which data to register?
i = 18;
% Test little difference from the groundtruth:
% initialparam = data{i}.gtparam + [.01 .1 .2 -30 -30 .5]; 
initialparam = [0 0 1.5 0 50 7]; 
% How sensitive each parameter is 
paramweight = [1 1 3 100 100 1.1];
paramshift = [.5 .5 .5 1000 500 1];
% 
% % perform registration
T = LKM.register(data{i}.data3D, data{i}.data2D, N, nbin, hsize,...
               mlModel{i},initialparam, ...
                true, paramweight, paramshift);
            