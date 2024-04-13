function [d] = emilyparencliticdeviation(m,b,x,y)

x_2 = (y-b)/m;
y_2 = m*x+b;

R_x = abs(x-x_2);
R_y = abs(y-y_2);

d=(R_x*R_y)/((R_x^2+R_y^2)^0.5);

end

