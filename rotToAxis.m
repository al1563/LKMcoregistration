function theta = rotToAxis(rot)

thetax = atan2(rot(3,2), rot(3,3));
thetay = atan2(-rot(3,1),sqrt(rot(3,2)^2 + rot(3,3)^2));
thetaz = atan2(rot(2,1), rot(1,1));

theta = [thetax thetay thetaz];

end