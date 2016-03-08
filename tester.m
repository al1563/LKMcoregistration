load data

nbin = 5;
hsize = 30;
N = 300;
%let N = 300 for about 1000 rows of training per set

%Train
mlModel1 = LKM.trainModel(data, N, nbin, hsize, trainindex);
initialparam = [0 0 1.7  .1 .1  7];
% 
% % perform registration
T = LKM.register(data{1}.data3D, data{1}.data2D, k, N, intv,...
               mlModel,initialparam, true);