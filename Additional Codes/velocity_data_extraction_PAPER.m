%% Load velocity data from ANSYS and model chamber dimensions
% ANSYS data import
import_data = importdata("C:\Users\Nida Qayyum\Downloads\Chamber_Simplified_Velocity_2");
ANSYS_data = import_data.data;
xc = ANSYS_data(:,2); xc = xc*100; % Chamber x-coordinate from ANSYS mesh
yc = ANSYS_data(:,4); yc = yc*100; % Chamber y-coordinate from ANSYS mesh
zc = ANSYS_data(:,3); zc = zc*100; % Chamber z-coordinate from ANSYS mesh
xv = ANSYS_data(:,6); % x-velocity
yv = ANSYS_data(:,8); % y-velocity
zv = ANSYS_data(:,7); % z-velocity
magv = ANSYS_data(:,5); % Velocity magnitude
cell_vol = ANSYS_data(:,12); % ANSYS cell volume

% Chamber dimensions
ch_l = 14; % Chamber length (cm)
ch_w = 9; % Chamber width (cm)
ch_h = 2.7; % Chamber height (cm)

% Define data extraction nodes
dx = 0.1; dy = 0.1; dz = 0.1; 
nodes_y = round(ch_l/dy) + 1;
nodes_x = round(ch_w/dx) + 1;
nodes_z = round(ch_h/dz) + 1;

ph = -10000; % Placeholder value to be used during data extraction issues

%% Loop for data extraction into defined nodes
% LOOP STARTS
for z=1:nodes_z
    if z == 1 || z == nodes_z % Top and bottom walls of chamber
        for y=1:nodes_y
            for x = 1:nodes_x
                vx(y,x,z) = 0; % No slip boundary condition
                vy(y,x,z) = 0; % No slip boundary condition
                vz(y,x,z) = 0; % No slip boundary condition
                vmag(y,x,z) = 0;
                vol_cell(y,x,z) = 0;
            end
        end
    % Layer 1 (flow obstruction by plate lid)
    elseif z>1 && z<12 % Height of water bath (Layer 1)
        for y=1:nodes_y
            if y == 1 || y == nodes_y % Back and front walls of chamber
                for x=1:nodes_x
                    vx(y,x,z) = 0; % No slip boundary condition
                    vy(y,x,z) = 0; % No slip boundary condition
                    vz(y,x,z) = 0; % No slip boundary condition
                    vmag(y,x,z) = 0;
                    vol_cell(y,x,z) = 0;
                end
            elseif y>=11 && y<=131 % Space occupied by water bath
                for x=1:nodes_x
                    if x == 1 || x>=11 && x<=81 || x == nodes_x % (x=1 and x=nodes_x are the left and right chamber walls. The middle condition is the space occupied by the water bath)
                        vx(y,x,z) = 0; % No slip boundary condition/no flow
                        vy(y,x,z) = 0; % No slip boundary condition/no flow
                        vz(y,x,z) = 0; % No slip boundary condition/no flow
                        vmag(y,x,z) = 0;
                        vol_cell(y,x,z) = 0;
                    else % Sides of the water bath with flow
                       % Coordinates of defined mesh
                        y_coor = (y-1)*dy;
                        x_coor = (x-1)*dx;
                        z_coor = (z-1)*dz;
    
                        % Match coordinates between defined and ANSYS mesh
                        coor_x = find(abs(xc-x_coor)<0.1);
                        coor_y = find(abs(yc-y_coor)<0.1);
                        coor_z = find(abs(zc-z_coor)<0.1);
                        coor_rows1 = intersect(coor_x,coor_y,'stable'); % Identify match between desired and ANSYS coordinate in 2-D plane
                        coor_rows = intersect(coor_z,coor_rows1,'stable'); % Identify match between desired and ANSYS coordinate in 3-D plane

                        if double(isempty(coor_rows)) == 0 % Match identified between desired and ANSYS coordinate
                            clear point_dist
                            for i = 1:length(coor_rows)
                                point_dist(i) = abs(zc(coor_rows(i)) - z_coor) + abs(xc(coor_rows(i)) - x_coor) + abs(yc(coor_rows(i)) - y_coor);
                            end
                            [~, target_row] = min(point_dist);
                            actual_row = coor_rows(target_row); % Find ANSYS coordinate closet to desired coordinate
    
                            % Extract data from identified ANSYS coordinate which matches desired coordinate
                            vx(y,x,z) = xv(actual_row);
                            vy(y,x,z) = yv(actual_row);
                            vz(y,x,z) = zv(actual_row);
                            vmag(y,x,z) = magv(actual_row);
                            vol_cell(y,x,z) = cell_vol(actual_row);
                        else
                            % Set variables to placeholder value if a coordinate match cannot be found
                            vx(y,x,z) = ph;
                            vy(y,x,z) = ph;
                            vz(y,x,z) = ph;
                            vmag(y,x,z) = ph;
                            vol_cell(y,x,z) = ph;
                        end
                    end
                end
            else % Unobstructed flow zone
                for x=1:nodes_x
                    if x == 1 || x == nodes_x % Left and right walls of chamber
                        vx(y,x,z) = 0; % No slip boundary condition
                        vy(y,x,z) = 0; % No slip boundary condition
                        vz(y,x,z) = 0; % No slip boundary condition
                        vmag(y,x,z) = 0;
                        vol_cell(y,x,z) = 0;
                    else
                        % Coordinates of defined mesh
                        y_coor = (y-1)*dy;
                        x_coor = (x-1)*dx;
                        z_coor = (z-1)*dz;
    
                        % Match coordinates between defined and ANSYS mesh
                        coor_x = find(abs(xc-x_coor)<0.1);
                        coor_y = find(abs(yc-y_coor)<0.1);
                        coor_z = find(abs(zc-z_coor)<0.1);
                        coor_rows1 = intersect(coor_x,coor_y,'stable'); % Identify match between desired and ANSYS coordinate in 2-D plane
                        coor_rows = intersect(coor_z,coor_rows1,'stable'); % Identify match between desired and ANSYS coordinate in 3-D plane

                        if double(isempty(coor_rows)) == 0 % Match identified between desired and ANSYS coordinate
                            clear point_dist
                            for i = 1:length(coor_rows)
                                point_dist(i) = abs(zc(coor_rows(i)) - z_coor) + abs(xc(coor_rows(i)) - x_coor) + abs(yc(coor_rows(i)) - y_coor);
                            end
                            [~, target_row] = min(point_dist);
                            actual_row = coor_rows(target_row); % Find ANSYS coordinate closet to desired coordinate
    
                            % Extract data from identified ANSYS coordinate which matches desired coordinate
                            vx(y,x,z) = xv(actual_row);
                            vy(y,x,z) = yv(actual_row);
                            vz(y,x,z) = zv(actual_row);
                            vmag(y,x,z) = magv(actual_row);
                            vol_cell(y,x,z) = cell_vol(actual_row);
                        else
                            % Set variables to placeholder value if a coordinate match cannot be found
                            vx(y,x,z) = ph;
                            vy(y,x,z) = ph;
                            vz(y,x,z) = ph;
                            vmag(y,x,z) = ph;
                            vol_cell(y,x,z) = ph;
                        end
                    end
                end
            end
        end
    % Layer 2 (flow obstruction by wells)
    elseif z>12 && z<20 % Layer 2 separated to account for inlet and outlet
        for y=1:nodes_y
            if y == 1 % Back wall of chamber
                for x=1:nodes_x
                    if x<28 || x>34 % Excluding chamber inlet
                        vx(y,x,z) = 0; % No slip boundary condition
                        vy(y,x,z) = 0; % No slip boundary condition
                        vz(y,x,z) = 0; % No slip boundary condition
                        vmag(y,x,z) = 0;
                        vol_cell(y,x,z) = 0;
                    else
                        % Coordinates of defined mesh
                        y_coor = (y-1)*dy;
                        x_coor = (x-1)*dx;
                        z_coor = (z-1)*dz;
    
                        % Match coordinates between defined and ANSYS mesh
                        coor_x = find(abs(xc-x_coor)<0.1);
                        coor_y = find(abs(yc-y_coor)<0.1);
                        coor_z = find(abs(zc-z_coor)<0.1);
                        coor_rows1 = intersect(coor_x,coor_y,'stable'); % Identify match between desired and ANSYS coordinate in 2-D plane
                        coor_rows = intersect(coor_z,coor_rows1,'stable'); % Identify match between desired and ANSYS coordinate in 3-D plane

                        if double(isempty(coor_rows)) == 0 % Match identified between desired and ANSYS coordinate
                            clear point_dist
                            for i = 1:length(coor_rows)
                                point_dist(i) = abs(zc(coor_rows(i)) - z_coor) + abs(xc(coor_rows(i)) - x_coor) + abs(yc(coor_rows(i)) - y_coor);
                            end
                            [~, target_row] = min(point_dist);
                            actual_row = coor_rows(target_row); % Find ANSYS coordinate closet to desired coordinate
    
                            % Extract data from identified ANSYS coordinate which matches desired coordinate
                            vx(y,x,z) = xv(actual_row);
                            vy(y,x,z) = yv(actual_row);
                            vz(y,x,z) = zv(actual_row);
                            vmag(y,x,z) = magv(actual_row);
                            vol_cell(y,x,z) = cell_vol(actual_row);
                        else
                            % Set variables to placeholder value if a coordinate match cannot be found
                            vx(y,x,z) = ph;
                            vy(y,x,z) = ph;
                            vz(y,x,z) = ph;
                            vmag(y,x,z) = ph;
                            vol_cell(y,x,z) = ph;
                        end
                    end
                end
            elseif y>=11 && y<=41 || y>=51 && y<=81 || y>=101 && y<=131 % Region obstructed by wells
                for x=1:nodes_x
                    if x == 1 || x>=11 && x<=41 || x>=51 && x<=81 || x == nodes_x % Region obstructed by well and left and right walls of chamber
                        vx(y,x,z) = 0; % No slip boundary condition/no flow
                        vy(y,x,z) = 0; % No slip boundary condition/no flow
                        vz(y,x,z) = 0; % No slip boundary condition/no flow
                        vmag(y,x,z) = 0;
                        vol_cell(y,x,z) = 0;
                    else
                        % Coordinates of defined mesh
                        y_coor = (y-1)*dy;
                        x_coor = (x-1)*dx;
                        z_coor = (z-1)*dz;
    
                        % Match coordinates between defined and ANSYS mesh
                        coor_x = find(abs(xc-x_coor)<0.1);
                        coor_y = find(abs(yc-y_coor)<0.1);
                        coor_z = find(abs(zc-z_coor)<0.1);
                        coor_rows1 = intersect(coor_x,coor_y,'stable'); % Identify match between desired and ANSYS coordinate in 2-D plane
                        coor_rows = intersect(coor_z,coor_rows1,'stable'); % Identify match between desired and ANSYS coordinate in 3-D plane

                        if double(isempty(coor_rows)) == 0 % Match identified between desired and ANSYS coordinate
                            clear point_dist
                            for i = 1:length(coor_rows)
                                point_dist(i) = abs(zc(coor_rows(i)) - z_coor) + abs(xc(coor_rows(i)) - x_coor) + abs(yc(coor_rows(i)) - y_coor);
                            end
                            [~, target_row] = min(point_dist);
                            actual_row = coor_rows(target_row); % Find ANSYS coordinate closet to desired coordinate
    
                            % Extract data from identified ANSYS coordinate which matches desired coordinate
                            vx(y,x,z) = xv(actual_row);
                            vy(y,x,z) = yv(actual_row);
                            vz(y,x,z) = zv(actual_row);
                            vmag(y,x,z) = magv(actual_row);
                            vol_cell(y,x,z) = cell_vol(actual_row);
                        else
                            % Set variables to placeholder value if a coordinate match cannot be found
                            vx(y,x,z) = ph;
                            vy(y,x,z) = ph;
                            vz(y,x,z) = ph;
                            vmag(y,x,z) = ph;
                            vol_cell(y,x,z) = ph;
                        end
                    end
                end
            elseif y == nodes_y % Front wall of chamber
                for x=1:nodes_x
                    if x<58 || x>64 % Excluding chamber outlet
                        vx(y,x,z) = 0; % No slip boundary condition
                        vy(y,x,z) = 0; % No slip boundary condition
                        vz(y,x,z) = 0; % No slip boundary condition
                        vmag(y,x,z) = 0;
                        vol_cell(y,x,z) = 0;
                    else
                         % Coordinates of defined mesh
                        y_coor = (y-1)*dy;
                        x_coor = (x-1)*dx;
                        z_coor = (z-1)*dz;
    
                        % Match coordinates between defined and ANSYS mesh
                        coor_x = find(abs(xc-x_coor)<0.1);
                        coor_y = find(abs(yc-y_coor)<0.1);
                        coor_z = find(abs(zc-z_coor)<0.1);
                        coor_rows1 = intersect(coor_x,coor_y,'stable'); % Identify match between desired and ANSYS coordinate in 2-D plane
                        coor_rows = intersect(coor_z,coor_rows1,'stable'); % Identify match between desired and ANSYS coordinate in 3-D plane

                        if double(isempty(coor_rows)) == 0 % Match identified between desired and ANSYS coordinate
                            clear point_dist
                            for i = 1:length(coor_rows)
                                point_dist(i) = abs(zc(coor_rows(i)) - z_coor) + abs(xc(coor_rows(i)) - x_coor) + abs(yc(coor_rows(i)) - y_coor);
                            end
                            [~, target_row] = min(point_dist);
                            actual_row = coor_rows(target_row); % Find ANSYS coordinate closet to desired coordinate
    
                            % Extract data from identified ANSYS coordinate which matches desired coordinate
                            vx(y,x,z) = xv(actual_row);
                            vy(y,x,z) = yv(actual_row);
                            vz(y,x,z) = zv(actual_row);
                            vmag(y,x,z) = magv(actual_row);
                            vol_cell(y,x,z) = cell_vol(actual_row);
                        else
                            % Set variables to placeholder value if a coordinate match cannot be found
                            vx(y,x,z) = ph;
                            vy(y,x,z) = ph;
                            vz(y,x,z) = ph;
                            vmag(y,x,z) = ph;
                            vol_cell(y,x,z) = ph;
                        end
                    end
                end
            else
                for x=1:nodes_x % Unobstructed flow regions in Layer 2
                    if x == 1 || x == nodes_x % Left and right walls of chamber
                        vx(y,x,z) = 0; % No slip boundary condition
                        vy(y,x,z) = 0; % No slip boundary condition
                        vz(y,x,z) = 0; % No slip boundary condition
                        vmag(y,x,z) = 0;
                        vol_cell(y,x,z) = 0;
                    else
                         % Coordinates of defined mesh
                        y_coor = (y-1)*dy;
                        x_coor = (x-1)*dx;
                        z_coor = (z-1)*dz;
    
                        % Match coordinates between defined and ANSYS mesh
                        coor_x = find(abs(xc-x_coor)<0.1);
                        coor_y = find(abs(yc-y_coor)<0.1);
                        coor_z = find(abs(zc-z_coor)<0.1);
                        coor_rows1 = intersect(coor_x,coor_y,'stable'); % Identify match between desired and ANSYS coordinate in 2-D plane
                        coor_rows = intersect(coor_z,coor_rows1,'stable'); % Identify match between desired and ANSYS coordinate in 3-D plane

                        if double(isempty(coor_rows)) == 0 % Match identified between desired and ANSYS coordinate
                            clear point_dist
                            for i = 1:length(coor_rows)
                                point_dist(i) = abs(zc(coor_rows(i)) - z_coor) + abs(xc(coor_rows(i)) - x_coor) + abs(yc(coor_rows(i)) - y_coor);
                            end
                            [~, target_row] = min(point_dist);
                            actual_row = coor_rows(target_row); % Find ANSYS coordinate closet to desired coordinate
    
                            % Extract data from identified ANSYS coordinate which matches desired coordinate
                            vx(y,x,z) = xv(actual_row);
                            vy(y,x,z) = yv(actual_row);
                            vz(y,x,z) = zv(actual_row);
                            vmag(y,x,z) = magv(actual_row);
                            vol_cell(y,x,z) = cell_vol(actual_row);
                        else
                            % Set variables to placeholder value if a coordinate match cannot be found
                            vx(y,x,z) = ph;
                            vy(y,x,z) = ph;
                            vz(y,x,z) = ph;
                            vmag(y,x,z) = ph;
                            vol_cell(y,x,z) = ph;
                        end
                    end
                end
            end
        end
    elseif z == 12 || z>19 && z<25
        for y=1:nodes_y
            if y == 1 || y == nodes_y % Front and back wall of chamber
                for x=1:nodes_x
                    vx(y,x,z) = 0; % No slip boundary condition
                    vy(y,x,z) = 0; % No slip boundary condition
                    vz(y,x,z) = 0; % No slip boundary condition
                    vmag(y,x,z) = 0;
                    vol_cell(y,x,z) = 0;
                end
            elseif y>=11 && y<=41 || y>=51 && y<=81 || y>=101 && y<=131 % Region obstructed by wells
                for x=1:nodes_x
                    if x == 1 || x>=11 && x<=41 || x>=51 && x<=81 || x == nodes_x % Region obstructed by wells and left and right walls of chamber
                        vx(y,x,z) = 0; % No slip boundary condition/no flow
                        vy(y,x,z) = 0; % No slip boundary condition/no flow
                        vz(y,x,z) = 0; % No slip boundary condition/no flow
                        vmag(y,x,z) = 0;
                        vol_cell(y,x,z) = 0;
                    else
                        % Coordinates of defined mesh
                        y_coor = (y-1)*dy;
                        x_coor = (x-1)*dx;
                        z_coor = (z-1)*dz;
    
                        % Match coordinates between defined and ANSYS mesh
                        coor_x = find(abs(xc-x_coor)<0.1);
                        coor_y = find(abs(yc-y_coor)<0.1);
                        coor_z = find(abs(zc-z_coor)<0.1);
                        coor_rows1 = intersect(coor_x,coor_y,'stable'); % Identify match between desired and ANSYS coordinate in 2-D plane
                        coor_rows = intersect(coor_z,coor_rows1,'stable'); % Identify match between desired and ANSYS coordinate in 3-D plane

                        if double(isempty(coor_rows)) == 0 % Match identified between desired and ANSYS coordinate
                            clear point_dist
                            for i = 1:length(coor_rows)
                                point_dist(i) = abs(zc(coor_rows(i)) - z_coor) + abs(xc(coor_rows(i)) - x_coor) + abs(yc(coor_rows(i)) - y_coor);
                            end
                            [~, target_row] = min(point_dist);
                            actual_row = coor_rows(target_row); % Find ANSYS coordinate closet to desired coordinate
    
                            % Extract data from identified ANSYS coordinate which matches desired coordinate
                            vx(y,x,z) = xv(actual_row);
                            vy(y,x,z) = yv(actual_row);
                            vz(y,x,z) = zv(actual_row);
                            vmag(y,x,z) = magv(actual_row);
                            vol_cell(y,x,z) = cell_vol(actual_row);
                        else
                            % Set variables to placeholder value if a coordinate match cannot be found
                            vx(y,x,z) = ph;
                            vy(y,x,z) = ph;
                            vz(y,x,z) = ph;
                            vmag(y,x,z) = ph;
                            vol_cell(y,x,z) = ph;
                        end
                    end
                end
            else
                for x=1:nodes_x
                    if x == 1 || x == nodes_x % Left and right walls of chamber
                        vx(y,x,z) = 0; % No slip boundary condition
                        vy(y,x,z) = 0; % No slip boundary condition
                        vz(y,x,z) = 0; % No slip boundary condition
                        vmag(y,x,z) = 0;
                        vol_cell(y,x,z) = 0;
                    else
                        % Coordinates of defined mesh
                        y_coor = (y-1)*dy;
                        x_coor = (x-1)*dx;
                        z_coor = (z-1)*dz;
    
                        % Match coordinates between defined and ANSYS mesh
                        coor_x = find(abs(xc-x_coor)<0.1);
                        coor_y = find(abs(yc-y_coor)<0.1);
                        coor_z = find(abs(zc-z_coor)<0.1);
                        coor_rows1 = intersect(coor_x,coor_y,'stable'); % Identify match between desired and ANSYS coordinate in 2-D plane
                        coor_rows = intersect(coor_z,coor_rows1,'stable'); % Identify match between desired and ANSYS coordinate in 3-D plane

                        if double(isempty(coor_rows)) == 0 % Match identified between desired and ANSYS coordinate
                            clear point_dist
                            for i = 1:length(coor_rows)
                                point_dist(i) = abs(zc(coor_rows(i)) - z_coor) + abs(xc(coor_rows(i)) - x_coor) + abs(yc(coor_rows(i)) - y_coor);
                            end
                            [~, target_row] = min(point_dist);
                            actual_row = coor_rows(target_row); % Find ANSYS coordinate closet to desired coordinate
    
                            % Extract data from identified ANSYS coordinate which matches desired coordinate
                            vx(y,x,z) = xv(actual_row);
                            vy(y,x,z) = yv(actual_row);
                            vz(y,x,z) = zv(actual_row);
                            vmag(y,x,z) = magv(actual_row);
                            vol_cell(y,x,z) = cell_vol(actual_row);
                        else
                            % Set variables to placeholder value if a coordinate match cannot be found
                            vx(y,x,z) = ph;
                            vy(y,x,z) = ph;
                            vz(y,x,z) = ph;
                            vmag(y,x,z) = ph;
                            vol_cell(y,x,z) = ph;
                        end
                    end
                end
            end
        end
    % Layer 3 (no flow obstruction)
    else
        for y=1:nodes_y
            if y == 1 || y == nodes_y % Front and back walls of chamber
                for x=1:nodes_x
                    vx(y,x,z) = 0; % No slip boundary condition
                    vy(y,x,z) = 0; % No slip boundary condition
                    vz(y,x,z) = 0; % No slip boundary condition
                    vmag(y,x,z) = 0;
                    vol_cell(y,x,z) = 0;
                end
            else
                for x=1:nodes_x
                    if x == 1 || x == nodes_x % Left and right walls of chamber
                        vx(y,x,z) = 0; % No slip boundary condition
                        vy(y,x,z) = 0; % No slip boundary condition
                        vz(y,x,z) = 0; % No slip boundary condition
                        vmag(y,x,z) = 0;
                        vol_cell(y,x,z) = 0;
                    else
                         % Coordinates of defined mesh
                        y_coor = (y-1)*dy;
                        x_coor = (x-1)*dx;
                        z_coor = (z-1)*dz;
    
                        % Match coordinates between defined and ANSYS mesh
                        coor_x = find(abs(xc-x_coor)<0.1);
                        coor_y = find(abs(yc-y_coor)<0.1);
                        coor_z = find(abs(zc-z_coor)<0.1);
                        coor_rows1 = intersect(coor_x,coor_y,'stable'); % Identify match between desired and ANSYS coordinate in 2-D plane
                        coor_rows = intersect(coor_z,coor_rows1,'stable'); % Identify match between desired and ANSYS coordinate in 3-D plane

                        if double(isempty(coor_rows)) == 0 % Match identified between desired and ANSYS coordinate
                            clear point_dist
                            for i = 1:length(coor_rows)
                                point_dist(i) = abs(zc(coor_rows(i)) - z_coor) + abs(xc(coor_rows(i)) - x_coor) + abs(yc(coor_rows(i)) - y_coor);
                            end
                            [~, target_row] = min(point_dist);
                            actual_row = coor_rows(target_row); % Find ANSYS coordinate closet to desired coordinate
    
                            % Extract data from identified ANSYS coordinate which matches desired coordinate
                            vx(y,x,z) = xv(actual_row);
                            vy(y,x,z) = yv(actual_row);
                            vz(y,x,z) = zv(actual_row);
                            vmag(y,x,z) = magv(actual_row);
                            vol_cell(y,x,z) = cell_vol(actual_row);
                        else
                            % Set variables to placeholder value if a coordinate match cannot be found
                            vx(y,x,z) = ph;
                            vy(y,x,z) = ph;
                            vz(y,x,z) = ph;
                            vmag(y,x,z) = ph;
                            vol_cell(y,x,z) = ph;
                        end
                    end
                end
            end
        end
    end

    % After velocity profile is extracted use nearest points to fill in missing data
    Plane_data = squeeze(vx(:,:,z)); % Filled in data for current z slice of chamber

    [row, col, ~] = find(Plane_data ~= 0 & Plane_data ~= ph); % Find data points extracted from ANSYS
    Plane_data_e = [row col]; % Extracted data points in plane

    [row_m, col_m, ~] = find(Plane_data == ph); % Find missing data points
    Plane_data_m = [row_m col_m]; % Missing data points in plane

    size_m = size(Plane_data_m); 
    no_elements_m = size_m(1); % Total number of missing points

    for n=1:no_elements_m % Cycle through all of the missing points
        ref_point = Plane_data_m(n,:);
        distances = sqrt(sum((Plane_data_e - ref_point).^2, 2));
        [~, sorted_indices] = sort(distances);
        nearest_points = Plane_data_e(sorted_indices(1:3), :); % 3 nearest points to missing point
        vx(ref_point(1),ref_point(2),z) = mean([vx(nearest_points(1,1),nearest_points(1,2),z),vx(nearest_points(2,1),nearest_points(2,2),z),vx(nearest_points(3,1),nearest_points(3,2),z)]);
        vy(ref_point(1),ref_point(2),z) = mean([vy(nearest_points(1,1),nearest_points(1,2),z),vy(nearest_points(2,1),nearest_points(2,2),z),vy(nearest_points(3,1),nearest_points(3,2),z)]);
        vz(ref_point(1),ref_point(2),z) = mean([vz(nearest_points(1,1),nearest_points(1,2),z),vz(nearest_points(2,1),nearest_points(2,2),z),vz(nearest_points(3,1),nearest_points(3,2),z)]);
        vmag(ref_point(1),ref_point(2),z) = mean([vmag(nearest_points(1,1),nearest_points(1,2),z),vmag(nearest_points(2,1),nearest_points(2,2),z),vmag(nearest_points(3,1),nearest_points(3,2),z)]);
        vol_cell(ref_point(1),ref_point(2),z) = mean([vol_cell(nearest_points(1,1),nearest_points(1,2),z),vol_cell(nearest_points(2,1),nearest_points(2,2),z),vol_cell(nearest_points(3,1),nearest_points(3,2),z)]);
    end

end
% LOOP ENDS

%% Post-processing
% Convert velocity data to cm/s
vx = vx*100;
vy = vy*100;
vz = vz*100;
vmag = vmag*100;

%% Left and right flow splits near inlets for each layer
% Determine x negative and x positive split for Layer 2, Element (1,3)
for z = 12:24 % Cycle through each z-slice
    l = z-11;
    Inlet_x = vx(1:10,26:36,z); % x-velocity for 1 cm plane in Layer 2 around inlet
    Inlet_vol = vol_cell(1:10,26:36,z); % ANSYS cell volume for 1 cm plane in Layer 2 around inlet
    
    Inlet_p_x = find(Inlet_x>0); % Find x-velocity > 0
    Vol_tot_p_x(l) = sum(Inlet_vol(Inlet_p_x));
    avg_pos_x(l) = sum(Inlet_x(Inlet_p_x).*Inlet_vol(Inlet_p_x))/Vol_tot_p_x(l); % Average x-velocity > 0 for z-slice
    
    Inlet_n_x = find(Inlet_x<0); % Find x-velocity < 0
    Vol_tot_n_x(l) = sum(Inlet_vol(Inlet_n_x));
    avg_neg_x(l) = sum(Inlet_x(Inlet_n_x).*Inlet_vol(Inlet_n_x))/Vol_tot_n_x(l); % Average x-velocity < 0 for z-slice
end
x_p_avg = mean(avg_pos_x); % Average x-velocity > 0 for Layer 2, Element (1,3)
x_n_avg = mean(avg_neg_x); % Average x-velocity < 0 for Layer 2, Element (1,3)

% Determine x negative and x positive split for Layer 1, Element (1,3)
clear Vol_tot_p_x Vol_tot_n_x avg_pos_x avg_neg_x
for z = 2:11
    l = z-1;
    Inlet_x = vx(1:10,26:36,z); % x-velocity for 1 cm plane in Layer 1 around inlet
    Inlet_vol = vol_cell(1:10,26:36,z); % ANSYS cell volume for 1 cm plane in Layer 1 around inlet
    
    Inlet_p_x = find(Inlet_x>0); % Find x-velocity > 0
    Vol_tot_p_x(l) = sum(Inlet_vol(Inlet_p_x));
    avg_pos_x(l) = sum(Inlet_x(Inlet_p_x).*Inlet_vol(Inlet_p_x))/Vol_tot_p_x(l); % Average x-velocity > 0 for z-slice
    
    Inlet_n_x = find(Inlet_x<0); % Find x-velocity < 0
    Vol_tot_n_x(l) = sum(Inlet_vol(Inlet_n_x));
    avg_neg_x(l) = sum(Inlet_x(Inlet_n_x).*Inlet_vol(Inlet_n_x))/Vol_tot_n_x(l); % Average x-velocity < 0 for z-slice
end
x_p_avg_l1 = mean(avg_pos_x); % Average x-velocity > 0 for Layer 1, Element (1,3)
x_n_avg_l1 = mean(avg_neg_x); % Average x-velocity < 0 for Layer 1, Element (1,3)

% Determine x negative and x positive split for Layer 3, Element (1,3)
clear Vol_tot_p_x Vol_tot_n_x avg_pos_x avg_neg_x
for z = 25:nodes_z-1
    l = z-24;
    Inlet_x = vx(1:10,26:36,z); % x-velocity for 1 cm plane in Layer 3 around inlet
    Inlet_vol = vol_cell(1:10,26:36,z); % ANSYS cell volume for 1 cm plane in Layer 3 around inlet
    
    Inlet_p_x = find(Inlet_x>0); % Find x-velocity > 0
    Vol_tot_p_x(l) = sum(Inlet_vol(Inlet_p_x));
    avg_pos_x(l) = sum(Inlet_x(Inlet_p_x).*Inlet_vol(Inlet_p_x))/Vol_tot_p_x(l); % Average x-velocity > 0 for z-slice
    
    Inlet_n_x = find(Inlet_x<0); % Find x-velocity < 0
    Vol_tot_n_x(l) = sum(Inlet_vol(Inlet_n_x));
    avg_neg_x(l) = sum(Inlet_x(Inlet_n_x).*Inlet_vol(Inlet_n_x))/Vol_tot_n_x(l); % Average x-velocity < 0 for z-slice
end
x_p_avg_l3 = mean(avg_pos_x); % Average x-velocity > 0 for Layer 3, Element (1,3)
x_n_avg_l3 = mean(avg_neg_x); % Average x-velocity > 0 for Layer 3, Element (1,3)

% Determine x negative and x positive split for Layer 3, Element (2,3)
clear Vol_tot_p_x Vol_tot_n_x avg_pos_x avg_neg_x
for z = 25:nodes_z-1
    l = z-24;
    Inlet_x = vx(10:20,26:36,z); % x-velocity for 1 cm plane in Layer 3, Row 2 around inlet
    Inlet_vol = vol_cell(10:20,26:36,z); % ANSYS cell volume for 1 cm plane in Layer 3, Row 2 around inlet
    
    Inlet_p_x = find(Inlet_x>0); % Find x-velocity > 0
    Vol_tot_p_x(l) = sum(Inlet_vol(Inlet_p_x));
    avg_pos_x(l) = sum(Inlet_x(Inlet_p_x).*Inlet_vol(Inlet_p_x))/Vol_tot_p_x(l); % Average x-velocity > 0 for z-slice
    
    Inlet_n_x = find(Inlet_x<0); % Find x-velocity < 0
    Vol_tot_n_x(l) = sum(Inlet_vol(Inlet_n_x));
    avg_neg_x(l) = sum(Inlet_x(Inlet_n_x).*Inlet_vol(Inlet_n_x))/Vol_tot_n_x(l); % Average x-velocity < 0 for z-slice
end
x_p_avg_l3_2 = mean(avg_pos_x); % Average x-velocity > 0 for Layer 3, Element (2,3)
x_n_avg_l3_2 = mean(avg_neg_x); % Average x-velocity < 0 for Layer 3, Element (2,3)

%% x-y flow splits for elements of each layer
% Determine x-split of flow for all defined node points
vmag_por = zeros(nodes_y,nodes_x,nodes_z);
for z=1:nodes_z
    for y=1:nodes_y
        for x=1:nodes_x
            vmag_por(y,x,z) = abs(vx(y,x,z))/(abs(vy(y,x,z)) + abs(vx(y,x,z)));
        end
    end
end

% x-y flow split for elements in Layer 2
z = 12:24; % Cycle through each z-slice
for y=1:14
    if y == 1 % x-y split for Layer 2, Element (1,5) 
        y_cell_s = (y-1)*10 + 2;
        y_cell_e = y_cell_s + 8;
        for x=1:9
            if x == 5
                x_cell_s = (x-1)*10 + 1;
                x_cell_e = x_cell_s + 9;
                x_y_por_l2(y,x) = mean(mean(mean(vmag_por(y_cell_s:y_cell_e,x_cell_s:x_cell_e,z))));
            else
                x_y_por_l2(y,x) = 0;
            end
        end
    elseif y == 5   
        y_cell_s = (y-1)*10 + 2;
        y_cell_e = y_cell_s + 8;
        for x=1:9
            if x == 1 % x-y split for Layer 2, Element (5,1)
                x_cell_s = (x-1)*10 + 2;
                x_cell_e = x_cell_s + 8;
                x_y_por_l2(y,x) = mean(mean(mean(vmag_por(y_cell_s:y_cell_e,x_cell_s:x_cell_e,z))));
            elseif x == 5 % x-y split for Layer 2, Element (5,5)
                x_cell_s = (x-1)*10 + 1;
                x_cell_e = x_cell_s + 9;
                x_y_por_l2(y,x) = mean(mean(mean(vmag_por(y_cell_s:y_cell_e,x_cell_s:x_cell_e,z))));
            else
                x_y_por_l2(y,x) = 0;
            end
        end
    elseif y == 9 || y == 10 
        y_cell_s = (y-1)*10 + 2;
        y_cell_e = y_cell_s + 8;
        for x=1:9
            if x == 1 % x-y split for Layer 2, Elements (9,1) and (10,1)
                x_cell_s = (x-1)*10 + 2;
                x_cell_e = x_cell_s + 8;
                x_y_por_l2(y,x) = mean(mean(mean(vmag_por(y_cell_s:y_cell_e,x_cell_s:x_cell_e,z))));
            elseif x == 9 % Only 1 outlet
                x_y_por_l2(y,x) = 0;
            else
                x_cell_s = (x-1)*10 + 1;
                x_cell_e = x_cell_s + 9;
                x_y_por_l2(y,x) = mean(mean(mean(vmag_por(y_cell_s:y_cell_e,x_cell_s:x_cell_e,z))));
            end
        end
    else
        for x=1:9
            x_y_por_l2(y,x) = 0;
        end
    end
end

% x-y flow split for elements in Layer 3
z = 25:nodes_z-1;
for y=1:14
    if y == 1
        y_cell_s = (y-1)*10 + 2;
        y_cell_e = y_cell_s + 8;
        for x=1:9
            if x == 1 || x == 9
                x_cell_s = (x-1)*10 + 2;
                x_cell_e = x_cell_s + 8;
                x_y_por_l3(y,x) = mean(mean(mean(vmag_por(y_cell_s:y_cell_e,x_cell_s:x_cell_e,z))));
            else
                x_cell_s = (x-1)*10 + 1;
                x_cell_e = x_cell_s + 9;
                x_y_por_l3(y,x) = mean(mean(mean(vmag_por(y_cell_s:y_cell_e,x_cell_s:x_cell_e,z))));
            end
        end
    elseif y == 14 % No flow in y-direction
        x_y_por_l3(y,x) = 0;
    else
        y_cell_s = (y-1)*10 + 1;
        y_cell_e = y_cell_s + 9;
        for x=1:9
            if x == 1 || x == 9
                x_cell_s = (x-1)*10 + 2;
                x_cell_e = x_cell_s + 8;
                x_y_por_l3(y,x) = mean(mean(mean(vmag_por(y_cell_s:y_cell_e,x_cell_s:x_cell_e,z))));
            else
                x_cell_s = (x-1)*10 + 1;
                x_cell_e = x_cell_s + 9;
                x_y_por_l3(y,x) = mean(mean(mean(vmag_por(y_cell_s:y_cell_e,x_cell_s:x_cell_e,z))));
            end
        end
    end
end