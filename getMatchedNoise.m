function noise = getMatchedNoise(subsetIndex)
%ADDNOISE Summary of this function goes here
%   Detailed explanation goes here

noise = zeros(1,6); 

xShiftParam = 10;
yShiftParam = 10;

% The first bit corresponds to scaling
if mod(subsetIndex,2)==1
    noise(6) = 0.2*rand()-0.1;
end;

% The second bit is for shifting in y direction
subsetIndex = floor(bitsra(subsetIndex,1));

if mod(subsetIndex,2)==1
    noise(5) = 2*yShiftParam*(1+rand())-yShiftParam;
end;

% THe third bit is for shifting in x direction
subsetIndex = floor(bitsra(subsetIndex,1));

if mod(subsetIndex,2)==1
    noise(5) = 2*xShiftParam*(1+rand())-xShiftParam;
end;

% The fourth bit is for rotation about z axis
subsetIndex = floor(bitsra(subsetIndex,1));

if mod(subsetIndex,2)==1
    noise(3) = pi/12*rand()-pi/24;
end;

% The fifth bit is for rotation about y axis
subsetIndex = floor(bitsra(subsetIndex,1));

if mod(subsetIndex,2)==1
    noise(2) = pi/24*rand()-pi/48;
end;

% The sixth bit is for rotation about x axis
subsetIndex = floor(bitsra(subsetIndex,1));

if mod(subsetIndex,2)==1
    noise(1) = pi/24*rand()-pi/48;
end;

end

