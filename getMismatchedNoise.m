function noise = getMismatchedNoise(subsetIndex)
%ADDNOISE Summary of this function goes here
%   Detailed explanation goes here

noise = zeros(1,6); 

xShiftParam = 100;
yShiftParam = 100;

% The first bit corresponds to scaling
if mod(subsetIndex,2)==1
    sign = 1;
    if rand()>0.5
        sign = -1;
    end;
    noise(6) = sign*(3 + 3*rand());
end;

% The second bit is for shifting in y direction
subsetIndex = floor(bitsra(subsetIndex,1));

if mod(subsetIndex,2)==1
    sign = 1;
    if rand()>0.5
        sign = -1;
    end;
    noise(5) = sign*yShiftParam*(1+rand());
end;

% THe third bit is for shifting in x direction
subsetIndex = floor(bitsra(subsetIndex,1));

if mod(subsetIndex,2)==1
    sign = 1;
    if rand()>0.5
        sign = -1;
    end;
    noise(4) = sign*xShiftParam*(1+rand());
end;

% The fourth bit is for rotation about z axis
subsetIndex = floor(bitsra(subsetIndex,1));

if mod(subsetIndex,2)==1
    sign = 1;
    if rand()>0.5
        sign = -1;
    end;
    noise(3) = sign*(pi/6+5*pi/6*rand());
end;

% The fifth bit is for rotation about y axis
subsetIndex = floor(bitsra(subsetIndex,1));

if mod(subsetIndex,2)==1
    sign = 1;
    if rand()>0.5
        sign = -1;
    end;
    noise(2) = sign*(pi/6+5*pi/6*rand());
end;

% The sixth bit is for rotation about x axis
subsetIndex = floor(bitsra(subsetIndex,1));

if mod(subsetIndex,2)==1
    sign = 1;
    if rand()>0.5
        sign = -1;
    end;
    noise(1) = sign*(pi/6+5*pi/6*rand());
end;

end

