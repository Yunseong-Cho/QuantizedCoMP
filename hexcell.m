function [coeff] = hexcell(SP)

Nu = SP.Nu; %Number of users
Nc = SP.Nc; 
r_cell = SP.dmax; %Radius of Hexagon
d0 = SP.d0; % meter
rho = SP.rho; % pathloss exponent
s = SP.s; % (dB) shadow fading standard deviation
LD = SP.LD;
NF = SP.NF;
BW = SP.BW;


CX = sqrt(3)*r_cell * [0, cos((0:Nc-2)*pi/3+pi/6)]; % x coordinates of centers
CY = sqrt(3)*r_cell * [0, sin((0:Nc-2)*pi/3+pi/6)]; % y coordinates of centers



points = zeros(Nu,Nc,2);


    for c = 1:Nc
        
        v_x = r_cell * cos((0:6)*pi/3) + CX(c); % vertax of hexagon centered at center
        v_y = r_cell * sin((0:6)*pi/3) + CY(c);
        
        %% Place Nu samples in the hexagon
        
        [p_x, p_y] = makepoints(SP);
        
        x = CX(c) + p_x;
        y = CY(c) + p_y;
        
        %% Coordinates of Nu points in a cell
        
        points(:,c,1) = x;
        points(:,c,2) = y;
        
%         hold on
%         plot(CX(c), CY(c), 'o');
%         plot(v_x, v_y);
%         plot(x, y, '.');
        
    end
    
    D = zeros(Nu,Nc,Nc);

    g = zeros(Nu, Nc, Nc); 
    coeff = zeros(Nu*Nc, Nc);


    for c = 1:Nc
       diff_x = points(:,:,1) - [CX(c)]; % distance is measured between corresponding cell center
       diff_y = points(:,:,2) - [CY(c)];

       D = sqrt(diff_x.^2 + diff_y.^2); % D(:,:,c)
       PL = 20*log10(4*3.14*d0/LD) + 10*rho*log10(D/d0) + s*randn(Nu,Nc); % 2.4GHz, 8dB Shadowing, n = 3, d0 = 100 m
       %PL = 128.1 + 37.6*log10(D/10^3) + 8;%*randn(Nu*Nc,2);

       Pnoise = -174 + 10*log10(BW) + NF;
       g = 10.^(-(PL + Pnoise)/10); % Large-scale fading 

       
       
       coeff(:,c) = g(:);
    end

end

 




