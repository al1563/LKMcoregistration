classdef LKM
    %LKM Local Kernel Matching functions
    
    properties
    end
    
    methods(Static)
%%        
        function [kernel, hist] = computeKernels(pts, nbin, hsize, graph)
        % pts - Nx2 matrix, each row is a pt and each col is a parameter
        % nbin - number of bin for each local histogram
        % hsize - length of side of histogram
        % graph - boolean. Graph some histograms? (for visualization)
        
        if(nargin < 4)
            graph = 0; % False by default
            if(nargin < 3)
                hsize = 100; % Arbitrarily chosen
                if(nargin < 2)
                    nbin = 8; % Arbitrarily chosen
                end
            end
        end
        
        Npts = size(pts);
        intv.size = hsize/nbin;
        intv.al = (-nbin/2*intv.size):intv.size:(nbin/2*intv.size);
        edge.x = repmat(pts(:,1),1,nbin+1) + repmat(intv.al,Npts(1),1);
        edge.y = repmat(pts(:,2),1,nbin+1) + repmat(intv.al,Npts(1),1);
        
        % Compute histograms and kernel
        kernel = [];
        for k = 1:Npts(1)
            hist{k} = zeros(length(edge.x(k,:))-1, length(edge.y(k,:))-1);
            
            for x_idx = 1:length(edge.x(k,:))-1
                for y_idx = 1:length(edge.y(k,:))-1
                    hist{k}(x_idx, y_idx) = ...
                     sum(pts(:,1) >= edge.x(k,x_idx) & pts(:,1) < edge.x(k,x_idx+1) & ...
                     pts(:,2) >= edge.y(k,y_idx) & pts(:,2) < edge.y(k,y_idx+1));
                end
            end
            % Note: hist takes local info (no bins for out of bound pts)
            
            hist{k} = hist{k}./sum(sum(hist{k})); % Normlze by tot # of pts
            kernel = [kernel; hist{k}(:)'];
        end
        
        % Graph 100 histograms for visualization
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
        function [sim, pts2] = similarity(pts3_t, pts2, mlModel, nbin, hsize)
            % Compute the similarity of point clouds using precomputed kernels
            % pts3D_t is the transformed and projected 3D point set, Nx2
            % pts2D is the 2D point set, Nx2 matrix
            % mlModel: pretrained model with eigvector and training kernel
            % nbin, hsize: Histogram parameters
            
            % Compute Kernels
            kernel3D = LKM.computeKernels(pts3_t, nbin, hsize);
            kernel2 = LKM.computeKernels(pts2, nbin, hsize);
            
             % Normalize the transformed 3D data and 2D data
             
%              pts3_t = pts3_t - repmat(mean(pts3_t),length(pts3_t),1);
%              hsize = hsize/max(max(abs(pts3_t)));
%              pts3_t(:,1) = pts3_t(:,1) / max(abs(pts3_t(:,1)));
%              pts3_t(:,2) = pts3_t(:,2) / max(abs(pts3_t(:,2)));
%                 
%              pts2 = pts2 - repmat(mean(pts2),length(pts2),1);
%              pts2(:,1) = pts2(:,1) / max(abs(pts2(:,1)));
%              pts2(:,2) = pts2(:,2) / max(abs(pts2(:,2)));
%                 
            % Nearest neighbor in 2D pts for each pt in transformed 3D pts
            [nearestNeighbors, nnDist] = knnsearch(pts2, pts3_t);
            
            % Pair up kernels. Note: Some 2D kernels will be repeated.
            kernel2D = kernel2(nearestNeighbors,:);
            kernels = [kernel3D kernel2D nnDist];
            
            minT = repmat(mlModel.minC,size(kernels,1),1);
            maxT = repmat(mlModel.maxC,size(kernels,1),1);
            kernels = (kernels - minT) ./ (maxT-minT);
            
            % Apply machine learning model
            Ktest = constructKernel(kernels, mlModel.training, mlModel.opts);
            Yhat = Ktest * mlModel.eigvector;
%             Yhat(Yhat > 1) = 1;
%             Yhat(Yhat < 0) = 0;
            
            % Similarity measure to minimize
            sim = -mean(Yhat);
        end
%%
        function mlModel = trainModel(data, options)
            % data: a cell array with fields data2D, data3D, center3D,
            %       center2D, rotation, scale, shift, and gtparam 
            % options: a structure containing
            %       hsize: size of sides of histogram
            %       N: size of point sample, 
            %       nbin: number of bins in local histogram
            %       trainingindex: array of index in data for training
            %       alpha, t: mlmodel parameters
            % TODO: add distance to threshold matching, erosion param
                
            training = [];
            labels = [];
            N = options.N;
            load matchedData;
            load mismatchedData;
            
            % Obtain matching and nonmatching training kernels
            for i = options.trainingindex
                i
                % Get rid of noise
                data3D = data{i}.data3D;
                data2D = data{i}.data2D;
                erode2D = Erode(data2D, 40, 50);
                erode3D = Erode(data3D, 10, 50);
                
                % Get N random points in 3D data
                R = randperm(length(erode3D));
                pts3 = erode3D(R(1:N),:);
                R = randperm(length(erode2D));
                pts2 = erode2D(R(1:N),:);
                
                % Transform 3D points by their groundtruths
                % param = [thetax thetay thetaz shiftx shifty scale]
                gtrot = rotToAxis(data{i}.rotation);
                gtparam = [gtrot data{i}.shift data{i}.scale];
                pts3_t = TransformPoint3D2D(gtparam, pts3);
                
                % Normalize the transformed 3D data and 2D data
%                 pts3_t = pts3_t - repmat(mean(pts3_t),length(pts3_t),1);
%                 hsize = hsize/max(max(abs(pts3_t)));
%                 pts3_t(:,1) = pts3_t(:,1) / max(abs(pts3_t(:,1)));
%                 pts3_t(:,2) = pts3_t(:,2) / max(abs(pts3_t(:,2)));
%                 
%                 pts2 = pts2 - repmat(mean(pts2),length(pts2),1);
%                 pts2(:,1) = pts2(:,1) / max(abs(pts2(:,1)));
%                 pts2(:,2) = pts2(:,2) / max(abs(pts2(:,2)));
                
                % Compute kernel for transformed 3D data and 2D data
                [kernel3D, hist3D] = LKM.computeKernels(pts3_t, options.nbin, options.hsize); 
                [kernel2, hist2D] = LKM.computeKernels(pts2, options.nbin, options.hsize); 
                
                % Create training set of pairs of MATCHING points.
                [nearestNeighbors, distances] = knnsearch(pts2, pts3_t);
                kernel2D = kernel2(nearestNeighbors,:);
                matching = [kernel3D kernel2D distances]; 
                % Threshold points by distance
                matching = matching( matching(:,size(matching,2))<20 , : );
                
                % more training for MATCHING points
                for j = 1:16
                    j
                    erodeMatched = Erode(matchedData{i,j},10,50);
                    R = randperm(length(erodeMatched));
                    if(length(erodeMatched) < N/50)
                        match3D = erodeMatched(R,:);
                    else
                        match3D = erodeMatched(R(1:N/50),:);
                    end;
                    [kernel3D,~] = LKM.computeKernels(match3D, options.nbin, options.hsize);
                    [nearestNeighbors, distances] = knnsearch(pts2, match3D);
                    kernel2D = kernel2(nearestNeighbors,:);
                    matching = [matching; kernel3D kernel2D distances]; 
                end

                % create training set of pairs of NONMATCHING points.
                nonmatching = LKM.nonmatching(data, i, erode2D, ...
                               erode3D, options);
                    % nonmatching = nonmatching(size(matching,1),:);

                % FUN NONMATCHING NESS
                for j = [1:16, 55:63]
                    j
                    erodeMismatched = Erode(mismatchedData{i,j},10,50);
                    R = randperm(length(erodeMismatched));
                    if(length(erodeMismatched) < N/50)
                        nonmatch3D = erodeMismatched(R,:);
                    else
                        nonmatch3D = erodeMismatched(R(1:N/50),:);
                    end;
                    [nonkernel3D,~] = LKM.computeKernels(nonmatch3D, options.nbin, options.hsize);
                    [nearestNeighbors, distances] = knnsearch(pts2, nonmatch3D);
                    nonkernel2D = kernel2(nearestNeighbors,:);
                    nonmatching = [nonmatching; nonkernel3D nonkernel2D distances]; 
                end

                training = [training; matching; nonmatching];
                labels = [labels; ones(size(matching,1), 1);...
                            zeros(size(nonmatching,1), 1)];

            end
            
            % Normalize columnwise
            minC = repmat(min(training),size(training,1),1); 
            maxC = repmat(max(training),size(training,1),1);
            training = (training - minC) ./ (maxC-minC);

            % Train machine learning model
            opts            = [];
            opts.ReguAlpha  = options.alpha;  % [0.0001, 0.1]
            opts.ReguType   = 'Ridge';
            opts.gnd        = labels;   % groundtruth (flair vector length)
            opts.KernelType = 'Gaussian';
            opts.t          = options.t;    % [4, 6, 10]  

            K = constructKernel(training, training, opts);
            [mlModel.eigvector, ~] = KSR(opts, labels, K);
            
            mlModel.opts = opts;
            mlModel.minC = minC(1,:);
            mlModel.maxC = maxC(1,:);
            mlModel.training = training;
            mlModel.labels = labels;
        end
%%
        function [T, sim] = register(data3D, data2D, options, mlModel, param, display, w, s)
            % data3D and data2D takes in the 3D (Nx3) and 2D(Nx2) pts
            % N - number of points to sample
            % options : nbin, hsize: histogram parameters
            % param is initialization of the transform
            % display - registration visualization
            % w and s: weight and shift to modify fminsearch
            if(nargin < 9)
                w = [1 1 1 1 1 1];
                s = [0 0 0 0 0 0];
            end
            
            nbin = options.nbin; hsize = options.hsize; 
            N = options.N;
            
            %%%%%%%%%%
            % Prepare data
            %%%%%%%%%%
  
            erode2D = Erode(data2D, 40, 50);
            erode3D = Erode(data3D, 10, 50);
                
            % Get N random points
            R = randperm(length(erode3D));
            pts3D = erode3D(R(1:N),:);
            pts3D_t = TransformPoint3D2D((param), pts3D);
            
            R = randperm(length(erode2D));
            pts2D = erode2D(R(1:N),:);
            
            %%%%%%%%%%
            % Compute kernels
            %%%%%%%%%%

            kernel3D = LKM.computeKernels(pts3D_t, nbin, hsize);
            kernel2D = LKM.computeKernels(pts2D, nbin, hsize); 

            %%%%%%%%%%
            % Set up display
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
            % Optimize transform with respect to local kernel matching
            %%%%%%%%%%
           
            % using anonymous function LKsim in order to use multiple parameters in
            % fminsearch
            simshow = @(x, varargin) x{varargin{:}};

            similarity = @(t) LKM.similarity(...
                TransformPoint3D2D(t, pts3D), ...
                pts2D, mlModel,nbin, hsize);
            show = @(t, sim) displayPoints(TransformPoint3D2D(t, pts3D),pts2D, sim);
            
            if display
                LKsim = @(t) simshow({show((t-s).*w, similarity((t-s).*w))}, 1);
            else
                LKsim = similarity;
            end
            
            % Better search function?
%             options = optimoptions(@fminunc,'Algorithm','quasi-newton');
%             T = fminunc(LKsim, param, options);
            param = param ./w + s; % Apply weighting
            T = fminsearch(LKsim, param, optimset('TolX',1e-2,'TolFun',1e-3));
            T = (T-s).*w;
            [sim, ~] = LKM.similarity(TransformPoint3D2D(T,pts3D), pts2D, mlModel, nbin, hsize)

        end
%%
        function  nonmatch = nonmatching(data, i, erode2D, erode3D, options)
            % Rather specific function to create nonmatching kernels
            % nonmatching(data, i, erode2D, erode3D, N, nbin, hsize)
            % data = dataset. i = which dataset? 
            % options:
            %   N number of points, nbin number of bins of histogram
            %   hsize how big histogram be
            N = options.N; nbin = options.nbin; hsize = options.hsize;
            
            %Get 2D data
            R = randperm(length(erode2D));
            pts2DRand = erode2D(R(1:N),:);
            nonkernel2 = LKM.computeKernels(pts2DRand, nbin, hsize);
            gtparam = data{i}.gtparam;

            % Four set of 3D data with different transforms
            R = randperm(length(erode3D)); 
            pts3DRand = erode3D(R(1:N),:);
            pts3DRand_t = TransformPoint3D2D(gtparam, pts3DRand);
            R = (randperm(floor(N/4)));
            kernel2D1 = nonkernel2(R(1:floor(N/4)),:);
            kernel3D1 = LKM.computeKernels(pts3DRand_t, nbin, hsize);
            kernel3D1 = kernel3D1(R(1:floor(N/4)),:);
            d1 = sum((pts3DRand_t-pts2DRand) .^ 2, 2) .^ 0.5;
            d1 = d1(R(1:floor(N/4)),:);

            R = randperm(length(erode3D));
            pts3DRand = erode3D(R(1:N),:);
            gtrot = gtparam(1:3);
            gtparam = [gtrot data{i}.shift-30 data{i}.scale];
            pts3DRand_t2 = TransformPoint3D2D(gtparam, pts3DRand);
            kernel3D2 = LKM.computeKernels(pts3DRand_t2, nbin, hsize);
            R = (randperm(floor(N/4)));
            kernel3D2 = kernel3D2(R(1:floor(N/4)),:);
            [NNparam2, d2] = knnsearch(pts2DRand, ...
                            pts3DRand_t2(R(1:floor(N/4)),:));
            kernel2D2 = nonkernel2(NNparam2,:);
            
            R = randperm(length(erode3D));
            pts3DRand = erode3D(R(1:N),:);
            gtparam = [gtrot(1:2)-.1 gtrot(3)-1.3 data{i}.shift data{i}.scale];
            pts3DRand_t3 = TransformPoint3D2D(gtparam, pts3DRand);
            kernel3D3 = LKM.computeKernels(pts3DRand_t3, nbin, hsize);
            R = (randperm(floor(N/4)));
            kernel3D3 = kernel3D3(R(1:floor(N/4)),:);
            [NNparam3, d3] = knnsearch(pts2DRand, ...
                            pts3DRand_t3(R(1:floor(N/4)),:));
            kernel2D3 = nonkernel2(NNparam3,:);
            
            
            R = randperm(length(erode3D));
            pts3DRand = erode3D(R(1:N),:);
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
            distances = [d1; d2; d3; d4];
            nonmatch = [nonkernel3, nonkernel2 distances];

        end
    end
    
end

