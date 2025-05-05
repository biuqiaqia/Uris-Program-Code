function theta = incident_angle_cal(XS,XR,Xinter)

c = norm(XS-XR);
a = norm(XS-Xinter);
b = norm(XR-Xinter);
theta = acosd((a^2+b^2-c^2)/(2*a*b))/2;%incidence angle
