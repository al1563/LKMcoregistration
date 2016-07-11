%Generate noise to data

%% Load the data first
load data;
numImages = 18;

%% Generate mismatched data by adding noise
numMismatchedCopies = 1;
numMismatchedData = numMismatchedCopies*63;
mismatchedData = cell(numImages, numMismatchedData);

for i = 1:numImages
    data3D = data{i}.data3D;
    
    for j=1:numMismatchedData;
        mismatchedParam= data{i}.gtparam + getMismatchedNoise(ceil(j/numMismatchedCopies));
        if(mismatchedParam(6) <= 0)
            mismatchedParam(6) = 0.1;
        end;
        
        mismatchedData{i, j} = TransformPoint3D2D(mismatchedParam,data3D);
    end;
end

%% Generate matched data by adding noise
numMatchedCopies = 1;
numMatchedData = numMatchedCopies*63;
matchedData = cell(numImages, numMatchedData);

for i = 1:numImages
    data3D = data{i}.data3D;
    
    for j=1:numMatchedData;
        matchedParam= data{i}.gtparam + getMatchedNoise(ceil(j/numMatchedCopies));
        if(matchedParam(6) <= 0)
            matchedParam(6) = 0.1;
        end;
        
        matchedData{i, j} = TransformPoint3D2D(matchedParam,data3D);
    end;
end

%% Plot the mismatched data
for i = 1:4
    figure;
    for j = 16*(i-1)+1:16*i
        if j==64
            break;
        end;
        subplot(4,4,mod(j,16)+1);
        plot(data{1}.data2D(:,1),data{1}.data2D(:,2),'r.');
        hold on;
        data3 = mismatchedData{1,j};
        plot(data3(:,1),data3(:,2),'b+');
    end
    title('Mismatched');
end

%% Plot the matched data
for i = 1:4
    figure;
    for j = 16*(i-1)+1:16*i
        if j==64
            break;
        end;
        subplot(4,4,mod(j,16)+1);
        plot(data{1}.data2D(:,1),data{1}.data2D(:,2),'r.');
        hold on;
        data3 = matchedData{1,j};
        plot(data3(:,1),data3(:,2),'b+');
    end
    title('Matched');
end

%%  Save the data
save('mismatchedData', 'mismatchedData');
save('matchedData', 'matchedData');
