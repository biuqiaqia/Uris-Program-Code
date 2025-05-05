%% plane_eqn
function [a, b, c, d] = plane_eqn(pt1, pt2, pt3)
A = [pt1;pt2;pt3];
A(:,1) = 1;
B = [pt1;pt2;pt3];
B(:,2) = 1;
C = [pt1;pt2;pt3];
C(:,3) = 1;
D = [pt1;pt2;pt3];
D = D * -1;

a = det(A);
b = det(B);
c = det(C);
d = det(D);
end
