% Transforms 3D points to a 2D set of points
%
% param = [thetax thetay thetaz shiftx shifty scale]
%
% P will first rotate around x, y, and z axis (mostly around z)
% Then P will project to 2D: Take (y, z) component
% Finally, P is scaled and shifted
% param = [0 0 0 0 0 1] is identity
%
function PT = TransformPoint3D2D(param,P)

% Rotation matrices
thetax = param(1); thetay = param(2); thetaz = param(3);
x_rot = [1 0 0;
         0 cos(thetax) -sin(thetax);
         0 sin(thetax) cos(thetax)];
y_rot = [cos(thetay) 0 sin(thetay);
         0 1 0;
         -sin(thetay) 0 cos(thetay)];
z_rot = [cos(thetaz) -sin(thetaz) 0;
         sin(thetaz) cos(thetaz) 0;
         0 0 1];

% Apply 3D rotation
P = z_rot * (y_rot * (x_rot * P'));
P = P';

% Project points
PT = [P(:,2) P(:,3)]; 
 
% Shift and Scale points
shiftx = param(4);
shifty = param(5);
scale = abs(param(6));

%center = mean(PT);
%PT = PT - repmat(center,size(PT,1),1);
M     = [scale 0 shiftx;
         0 scale shifty]; 
PT = (M * [PT'; ones(1,size(PT,1))] )';

end
