function [] = plot_ellipse( covariance,mean,color)


[eigenvec, eigenval] = eigs(covariance);
eigenval=sort(diag(eigenval),'descend');

% Get the index of the largest eigenvector: This is the evector associated
% to the highest evalue

max_evc = eigenvec(:, 1);

% Get the largest eigenvalue
max_evl = eigenval(2);
min_evl = eigenval(1);
min_evc= eigenvec(2,:);

% Calculate the angle between the x-axis and the largest eigenvector
angle = atan2(max_evc(2), max_evc(1));
 
% This angle is between -pi and pi.
% Let's shift it such that the angle is between 0 and 2pi
if(angle < 0)   
    angle = angle + 2*pi;
end

% Get the 90% confidence interval error ellipse
chisquare_val = chi2inv(0.9,2);
theta_grid = linspace(0,2*pi);
phi = angle;
X0  = mean(1);
Y0  = mean(2);
a   = sqrt(chisquare_val*max_evl);
b   = sqrt(chisquare_val*min_evl);

% the ellipse in x and y coordinates 
ellipse_x_r  = a*cos( theta_grid );
ellipse_y_r  = b*sin( theta_grid );

%Define a rotation matrix
R = [ cos(phi) sin(phi); -sin(phi) cos(phi) ];

%let's rotate the ellipse to some angle phi
r_ellipse = [ellipse_x_r;ellipse_y_r]' * R;

% Draw the error ellipse
%figure;
plot(r_ellipse(:,1) + X0,r_ellipse(:,2) + Y0,color)

end

