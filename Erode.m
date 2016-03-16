function erode3D = Erode(data3D, extent, threshold)

is2D = 0;

if(size(data3D,2) < 3)
    is2D = 1;
    data3D = [data3D; ones(1,size(data3D,2))];
end
        
%First let's manually set a range box
%range = zeros(1,2);
%range(1) = minval+maxval/extent?
range = [extent, extent, extent];
nearby = zeros(size(data3D,1),1);
erode3D = zeros(0,3);

if(size(data3D,2) == 2)
    data3D = [data3D, ones(size(data3D,1),1)];
end

%More efficient algorithm?
for i = 1:size(data3D,1)
    x_cmp = abs(data3D(:, 1)-data3D(i, 1)) < range(1);
    data3Dx = data3D(x_cmp,:);
    neighbors = abs(data3Dx(:, 2)-data3D(i, 2)) < range(2);

    if(~is2D)
     data3Dxy = data3Dx(neighbors,:);
     neighbors = abs(data3Dx(:, 3)-data3D(i, 3)) < range(3);
    end
    
    nearby(i) = nearby(i) + sum(neighbors);
end


highz = max(data3D);
diff = max(data3D) - min(data3D);
zcutoff = highz(3)-diff*.35;
zcutoff2 = highz(3) - diff*.25;

erode3D = [erode3D; data3D(nearby > threshold,:)];

if(is2D)
    erode3D = erode3D(:,1:2);
end

end