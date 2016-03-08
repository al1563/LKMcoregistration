load data

nbin = 8;
hsize = 100;
N = 300;
%let N = 300 for about 1000 rows of training per set

%Train
mlModel = LKM.trainModel(data, N, nbin, hsize, trainindex);

i = 14;
initialparam = data{i}.gtparam + [.01 .1 .2 -30 -30 .5];
paramweight = [1 1 3 100 100 1.1];
paramscale = [.5 .5 .5 10 10 1];
initialparam = initialparam ./ paramweight;
initialparam = initialparam + paramscale;
% 
% % perform registration
T = LKM.register(data{i}.data3D, data{i}.data2D, N, nbin, hsize,...
               mlModel{2}.hbinsize{2,1},initialparam, ...
                true, paramweight, paramscale);