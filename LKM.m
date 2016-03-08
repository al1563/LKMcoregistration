classdef LKM
    %LKM Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods(Static)
%%        
        function [kernel, hist] = computeKernels(pts, nbin, hsize, graph)
        % pts - 2D, each row is one point and each column is a parameter
        % nbin - number of bin for each local histogram
        % hsize - how big total histogram be
        % graph - boolean. Graph the histogram and kernel?
        
        if(nargin < 4)
            graph = 0; % False by default
            if(nargin < 3)
                hsize = 100; % Arbitrarily chosen
            end
        else
            set(0,'DefaultFigureWindowStyle','docked'); % Figure is docked
        end
        
        Npts = size(pts);
        intv.size = hsize/nbin;
        intv.al = (-nbin/2*intv.size):intv.size:(nbin/2*intv.size);
        edge.x = repmat(pts(:,1),1,nbin+1) + repmat(intv.al,Npts(1),1);
        edge.y = repmat(pts(:,2),1,nbin+1) + repmat(intv.al,Npts(1),1);
        
        % Compute histograms and kernel
        kernel = [];
        for k = 1:Npts(1)
            hist{k} = hist2(pts(:,1),pts(:,2),edge.x(k,:),edge.y(k,:));
            %nbin+3 x nbin+3 matrix. Outer edges give out of bound values
            hist{k} = hist{k}(2:nbin+1,2:nbin+1); %Take local information
            hist{k} = hist{k}./sum(sum(hist{k})); %Normalize by numpoints
            kernel = [kernel; hist{k}(:)'];
        end
        
        if(graph)
            figure;
            for j= 1:100
                subplot(10,10,j);
                imshow(flip(hist{j}));
                caxis([0 max(max(hist{j}))]);
                title([num2str(j)]);
            end
            colormap default;
        end
        
        end
%%
        function [sim, pts2D] = similarity(pts3D_t, pts2D, mlModel, k, intv)
            % compute the similarity of point clouds using precomputed kernels
            kernel3D = LKM.computeKernels(pts3D_t, k, intv);
            kernel2 = LKM.computeKernels(pts2D,k, intv);
            
            % compute the distance from each point in 2D to nearest
            % neighbor in 3D transformation
            [nearestNeighbors, nnDist] = knnsearch(pts2D, pts3D_t);

            % pair up the kernels of nearest neighbors
            % note some bKernels will be paired with the same aKernel
            %pts2D = data2D(nearestNeighbors,:);
            kernel2D = kernel2(nearestNeighbors,:);
            kernels = [kernel3D kernel2D nnDist];

            kernels = normc(kernels);
            % apply machine learning model to 'kernels'
            Ktest = constructKernel(kernels, mlModel.training, mlModel.opts);
            Yhat = Ktest * mlModel.eigvector;

            sim = -mean(Yhat);
        end
%%
        function mlModel = trainModel(data, N, nbin, hsize, trainindex, a)
            % N: size of point sample
            % nbin: number of bins in local histogram
            % hsize: size of histogram
            % data is a cell array with fields data2D, data3D, center3D,
            % center2D, rotation, scale, and shift

            training = [];
            labels = [];
            traininga = [];
            labelsa = [];

            for i = [a, trainindex]
                % Gets rsid of noise...
                data3D = data{i}.data3D;
                data2D = data{i}.data2D;
                erode2D = Erode(data2D, 40, 50);
                erode3D = Erode(data3D, 10, 50);
                
                % Get N random points in data3D
                R = randperm(length(erode3D));
                pts3 = erode3D(R(1:2*N),:);
                R = randperm(length(erode2D));
                pts2 = erode2D(R(1:2*N),:);
                
                % transform 3D points by their groundtruths
                % param = [thetax thetay thetaz shiftx shifty scale]
                gtrot = rotToAxis(data{i}.rotation);
                gtparam = [gtrot data{i}.shift data{i}.scale];
                pts3_t = TransformPoint3D2D(gtparam, pts3);
                
                % Compute kernel for 3D data
                % Each row is a pt. Each col is dist to nearest neighbor
                [kernel3D, hist3D] = LKM.computeKernels(pts3_t, nbin, hsize); 
                [kernel2, hist2D] = LKM.computeKernels(pts2, nbin, hsize); 
                
                % Create training set of pairs of MATCHING points.
                [nearestNeighbors, distances] = knnsearch(pts2, pts3_t);
                kernel2D = kernel2(nearestNeighbors,:);
                matching = [kernel3D kernel2D distances]; 
                % Threshold points by distance
                matching = matching(matching(:,size(matching,2))<20,...
                            1:size(matching,2)-1);

                % create training set of pairs of NONMATCHING points.
                nonmatching = LKM.nonmatching(data, i, erode2D, ...
                                erode3D, N, nbin, hsize);
            
                if(i==a)
                    traininga = [traininga; matching; nonmatching];
                    labelsa = [labelsa; ones(size(matching,1), 1);...
                                zeros(size(nonmatching,1), 1)];
                else
                    training = [training; matching; nonmatching];
                    labels = [labels; ones(size(matching,1), 1);...
                                zeros(size(nonmatching,1), 1)];
                end

            end
            
            training = normc(training);
            traininga = normc(traininga);

            % perform machine learning on the training set
            opts            = [];
            opts.ReguAlpha  = 0.1;  % [0.0001, 0.1]
            opts.ReguType   = 'Ridge';
            opts.gnd        = labels;   % groundtruth (flair vector length)
            opts.KernelType = 'Gaussian';
            opts.t          = 6;    % [4, 10]  

            K = constructKernel(training, training, opts);
            [mlModel.eigvector, ~] = KSR(opts, labels, K);
            mlModel.opts = opts;
            mlModel.training = training;
            mlModel.labels = labels;
            mlModel.traininga = traininga;
            mlModel.labelsa = labelsa;
        end
%%
        function T = register(data3D, data2D, k, N, intv, mlModel, param, display)
            
            % A and B are 3Ddata and 2Ddata respectively.
            % k is the number of nearest-neighbor points to use in each kernel
            % 
            % T is a 2x3 (?) matrix representing an affine transform which optimally
            % maps A to B
            
            %%%%%%%%%%
            % initialize transform
            %%%%%%%%%%

            %param = [thetax thetay thetaz shiftx shifty scale]
            %param = [0 0 2.4 50 50 6];
            param = [param(1) param(2) param(3) param(4) param(5) param(6)];
            
            erode2D = Erode(data2D, 40, 50);
            erode3D = Erode(data3D, 10, 50);
                
            % Get N random points
            R = randperm(length(erode3D));
            pts3D = erode3D(R(1:2*N),:);
            pts3D_t = TransformPoint3D2D(param, pts3D);
            
            R = randperm(length(erode2D));
            pts2D = erode2D(R(1:2*N),:);
            %Already centered at origin for the time being
            
            %%%%%%%%%%
            % compute kernels
            %%%%%%%%%%

            kernel3D = LKM.computeKernels(pts3D_t, k, intv);
            kernel2D = LKM.computeKernels(pts2D, k, intv); %unnecessary

            %%%%%%%%%%
            % set up display
            %%%%%%%%%%

            if display
                subplot(1, 2, 1)
                displayPoints(pts3D_t,pts2D);
                set(gca,'FontSize',16)
                title('Initial position')
                drawnow
                subplot(1, 2, 2)
                set(gca,'FontSize',16)
            end

            %%%%%%%%%%
            % optimize transform with respect to local kernel matching
            %%%%%%%%%%

            % using anonymous function LKsim in order to use multiple parameters in
            % fminsearch
            simshow = @(x, varargin) x{varargin{:}};

            similarity = @(t) LKM.similarity(...
                TransformPoint3D2D([t(1) t(2) t(3) t(4) t(5) t(6)], pts3D), ...
                pts2D, mlModel,k, intv) + 10*abs(t(6)-6);
            show = @(t, sim) displayPoints(TransformPoint3D2D(t, pts3D),pts2D, sim);
            
            if display
                LKsim = @(t) simshow({show(t, similarity(t))}, 1);
            else
                LKsim = similarity;
            end
            
            % maybe a numerical gradient descent would be faster than fminsearch
            T = fminsearch(LKsim, param);

        end
%%
        function  nonmatch = nonmatching(data, i, erode2D, erode3D, N, nbin, hsize)
            % nonmatching(data, i, erode2D, erode3D, N, nbin, hsize)
            % data = dataset. i = which dataset? 
            % N number of points, nbin - number of bins of histogram
            % hsize - how big histogram be
            
            %Get 2D data
            R = randperm(length(erode2D));
            pts2DRand = erode2D(R(1:2*N),:);
            nonkernel2 = LKM.computeKernels(pts2DRand, nbin, hsize);

            % Four set of 3D data with different transforms
            R = randperm(length(erode3D)); pts3DRand = erode3D(R(1:2*N),:);
            gtrot = rotToAxis(data{i}.rotation);
            gtparam = [gtrot data{i}.shift data{i}.scale];
            pts3DRand_t = TransformPoint3D2D(gtparam, pts3DRand);
            R = (randperm(floor(N/2)));
            kernel2D1 = nonkernel2(R(1:floor(N/2)),:);
            kernel3D1 = LKM.computeKernels(pts3DRand_t, nbin, hsize);
            kernel3D1 = kernel3D1(R(1:floor(N/2)),:);
            d1 = sum((pts3DRand_t-pts2DRand) .^ 2, 2) .^ 0.5;
            d1 = d1(R(1:floor(N/2)),:);

            R = randperm(length(erode3D));
            pts3DRand = erode3D(R(1:2*N),:);
            gtparam = [gtrot data{i}.shift-30 data{i}.scale];
            pts3DRand_t2 = TransformPoint3D2D(gtparam, pts3DRand);
            kernel3D2 = LKM.computeKernels(pts3DRand_t2, nbin, hsize);
            R = (randperm(floor(N/2)));
            kernel3D2 = kernel3D2(R(1:floor(N/2)),:);
            [NNparam2, d2] = knnsearch(pts2DRand, ...
                            pts3DRand_t2(R(1:floor(N/2)),:));
            kernel2D2 = nonkernel2(NNparam2,:);
            
            R = randperm(length(erode3D));
            pts3DRand = erode3D(R(1:2*N),:);
            gtparam = [gtrot(1:2)-.1 gtrot(3)-1.3 data{i}.shift data{i}.scale];
            pts3DRand_t3 = TransformPoint3D2D(gtparam, pts3DRand);
            kernel3D3 = LKM.computeKernels(pts3DRand_t3, nbin, hsize);
            R = (randperm(floor(N/2)));
            kernel3D3 = kernel3D3(R(1:floor(N/2)),:);
            [NNparam3, d3] = knnsearch(pts2DRand, ...
                            pts3DRand_t3(R(1:floor(N/2)),:));
            kernel2D3 = nonkernel2(NNparam3,:);
            
            
            R = randperm(length(erode3D));
            pts3DRand = erode3D(R(1:2*N),:);
            gtparam = [gtrot(1:2) gtrot(3) data{i}.shift data{i}.scale*.5];
            pts3DRand_t4 = TransformPoint3D2D(gtparam, pts3DRand);
            kernel3D4 = LKM.computeKernels(pts3DRand_t4, nbin, hsize);
            R = (randperm(floor(N/2)));
            kernel3D4 = kernel3D4(R(1:floor(N/2)),:);
            [NNparam4, d4] = knnsearch(pts2DRand, ...
                            pts3DRand_t4(R(1:floor(N/2)),:));
            kernel2D4 = nonkernel2(NNparam4,:);

            % Compute kernels
            
            nonkernel3 = [kernel3D4; kernel3D3; kernel3D2; kernel3D1];
            nonkernel2 = [kernel2D4; kernel2D3; kernel2D2; kernel2D1];

            %distances = sum((pts3DRand_t-pts2DRand) .^ 2, 2) .^ 0.5;
            %distances = [d1; d2; d3; d4];
            nonmatch = [nonkernel3, nonkernel2];

        end
    end
    
end

