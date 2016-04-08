% load data
% addpath([pwd,'/SR']);
% 
% nbin = 10; 
% hsize = 150;
% % N = 500 for about 800 rows of training kernel per set
% N = 500;
% 
% % Train machine learning model
% mlModel = LKM.trainModel(data, N, nbin, hsize);

i = 15; % Data to test registration
% Test little difference from the groundtruth:
%initialparam = data{i}.gtparam + [.01 .1 .2 -30 -30 .5]; 
initialparam = [0 0 .5 -50 100 6];
% How sensitive each parameter is 
paramweight = [1 1 50 100 100 1.1];
paramshift = [.5 .5 1 10 10 1];

% % perform registration
T = LKM.register(data{i}.data3D, data{i}.data2D, N, nbin, hsize,...
               mlModel,initialparam, 1, paramweight, paramshift);
            