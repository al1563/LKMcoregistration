function norm = norm01(vector)
% Currently normalizes the elements of the input vector between 0 and 1
maxv = max(vector);
minv = min(vector);
norm = (vector - minv)/(maxv - minv);
end