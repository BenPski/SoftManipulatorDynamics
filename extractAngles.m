function out = extractAngles(M)
    %https://d3cw3dd2w32x2b.cloudfront.net/wp-content/uploads/2012/07/euler-angles1.pdf
    theta1 = atan2(M(2,3),M(3,3));
    c2 = sqrt(M(1,1)^2+M(1,2)^2);
    theta2 = atan2(-M(1,3),c2);
    s1 = sin(theta1);
    c1 = cos(theta1);
    theta3 = atan2(s1*M(3,1)-c1*M(2,1),c1*M(2,2)-s1*M(3,2));
    out = [theta1;theta2;theta3];
end