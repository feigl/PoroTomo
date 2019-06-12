function area = triangular_area(x,y,z)
% calculate the area of a 3-D triangle
%https://www.mathworks.com/matlabcentral/answers/14928-area-of-triangle
% x = rand(3,1);
% y = rand(3,1);
% z = rand(3,1);
% fill3(x,y,z,'r')
% x = reshape(x,1,3);
% y = reshape(y,1,3);
% z = reshape(z,1,3);
% figure;
% plot3(x,y,z,'k-o');
% hold on;
% fill3(x,y,z,'cyan');
% ons = [1 1 1];
% area= 0.5*sqrt(det([x;y;ons])^2 + det([y;z;ons])^2 + det([z;x;ons])^2);
A=[x(1),y(1),z(1)];
B=[x(2),y(2),z(2)];
C=[x(3),y(3),z(3)];
area = 0.5 * norm(cross(B-A,C-A));
return
end

