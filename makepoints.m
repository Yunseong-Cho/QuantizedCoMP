function [c_x, c_y] = makepoints(SP)

N = SP.Nu;
R_max = SP.dmax;     %Radius of Hexagon
R_min = SP.dmin;
%Define the vertexes of the hexagon. They for angles 0, 60, 120, 180, 240 and 300 withe origin.
%Vertexes
while(1)
    v_x_max = R_max * cos((0:6)*pi/3);
    v_y_max = R_max * sin((0:6)*pi/3);

    v_x_min = R_min * cos((0:6)*pi/3);
    v_y_min = R_min * sin((0:6)*pi/3);
    %The method used here is to generate many points in a square and choose N points that fall within the hexagon
    %Generate 4N random points with square that is 2R by 2R
    c_x = R_max-rand(1, 4*N)*2*R_max;
    c_y = R_max-rand(1, 4*N)*2*R_max;
    %There is a command in MATLAB inploygon. 
    %The command finds points within a polygon region.
    %get the points within the polygon
    IN_max = inpolygon(c_x, c_y, v_x_max, v_y_max);

    %drop nodes outside the hexagon
    c_x = c_x(IN_max);
    c_y = c_y(IN_max);

    IN_min = inpolygon(c_x, c_y, v_x_min, v_y_min);


    c_x = c_x(~IN_min);
    c_y = c_y(~IN_min);

    if length(c_x)>=N
        %randomly choose only N points
        idx = randperm(length(c_x));
        c_x = c_x(idx(1:N));
        c_y = c_y(idx(1:N));
        break
    end

end


% plot(c_x, c_y, 'r*');
% hold on;
% plot(v_x_max,v_y_max);
% plot(v_x_min,v_y_min);
% 
% axis square;
% hold off
end