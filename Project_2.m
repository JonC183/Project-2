%% Part A

startup_rvc;

clear
clc
close all

% Load File

filename = '4.jpg'

image = imread(filename);
image_hsv = rgb2hsv(image);

% Load from Camera
% 
% camList = webcamlist;
% cam = webcam(1);
% 
% preview(cam)

%%% Don't Start Until Ready
prompt = 'Start';
inputCommand = input(prompt)
% 
% 
% while true
%   image = snapshot(cam)
    %%% CODE FOR EXTRACTING ESSENTIAL OBJECTS HERE %%%
% end
% 

rgb = prism(6);
hsv = rgb2hsv(rgb);

figure(1);
imshow(image_hsv);

% Mask for Purple Circles
h_purple = [0.75 0.79];
s_purple = [0.54 0.58];
v_purple = [0.4 0.6];

mask_purple = createMaskAndShow(image_hsv,h_purple,s_purple,v_purple);
centers_purple = findCenters(mask_purple);

hold on
plot(centers_purple(:,1),centers_purple(:,2),'*r');
hold off

% Mask for Yellow Corners
h_yellow = [0.07 0.11];
s_yellow = [0.5 0.65];
v_yellow = [0.7 0.9];

mask_yellow = createMaskAndShow(image_hsv,h_yellow,s_yellow,v_yellow);
centers_yellow = findCenters(mask_yellow);

hold on
plot(centers_yellow(:,1),centers_yellow(:,2),'*r');
hold off

%%%% Conversion for world frame %%%%%
point1 = [-250, 75];
point2 = [-250, -525];
point3 = [-900, 75];
point4 = [-900, -525];

world = [point1 ; point3 ; point2 ; point4];

outputFrameWorld = [abs(900-250) abs(525+75)];

tform_world = fitgeotrans(centers_purple,world,'projective');
board_trans_world = imwarp(image,tform_world,'OutputView',imref2d(outputFrameWorld));

%%%% Conversion for Image %%%%
point1 = [0 0];
point2 = [380 0];
point3 = [0 590];
point4 = [380 590];

world_img = [point1 ; point3 ; point2 ; point4];

outputFrameImg = [590 380];

tform_img = fitgeotrans(centers_yellow,world_img,'projective');

board_trans_img = imwarp(image_hsv,tform_img,'OutputView',imref2d(outputFrameImg));
board_trans_img_rgb = imwarp(image,tform_img,'OutputView',imref2d(outputFrameImg));

figure(1);
imshow(board_trans_img);

% Get Coordinates of Corners from Image
board_corners_img = centers_yellow;

% Convert to World Coordinates
board_corners_world = transformPointsForward(tform_world,board_corners_img);

% Produce Transform for Transformed Image to World
% trans_T_world = trans_T_original*original_T_world
% trans_T_original = original_T_trans'
% tform_img_to_world = tform_img.T'*tform_world

% Test to see inverse of tform_img
P_img = [0 0];

P_original = transformPointsInverse(tform_img,P_img);
P_world = transformImgToWorld(tform_img,tform_world,P_img);

figure(4)
imshow(image)
hold on
plot(P_original(1),P_original(2),'g*','MarkerSize',30)

num_rows = 5;
num_cols = 8;

cols_size = round(outputFrameImg(1)/num_cols);
rows_size = round(outputFrameImg(2)/num_rows);

square_center_img = {};
square_kernel_img = {};
square_center_world = {};

figure(1)
for row = 1:num_rows
    for col = 1:num_cols
        point = [(rows_size*row - round(rows_size/2))...
            (cols_size*col - round(cols_size/2))];

        square_center_img{row,col} = point;
        
        hold on
        % Plot the point
        plot(square_center_img{row,col}(1),square_center_img{row,col}(2),'r*');

        hold off
        square_center_world{row,col} = transformImgToWorld(tform_img,tform_world,point);
    end
end

%%%% Detect Red Pucks %%%%
h_red = [0.96 0.98];
s_red = [0.6 0.85];
v_red = [0.6 0.8];

mask_red = createMaskAndShow(board_trans_img,h_red,s_red,v_red);

[red_puck_img_coord,red_puck_cell_coord,red_puck_world_coord] = ...
    findColouredPuck(mask_red,square_center_img,2,tform_img,tform_world,'red');

hold on
plot(red_puck_img_coord(:,1),red_puck_img_coord(:,2),'r*','MarkerSize',30)

%%%% Detect Blue Pucks %%%%
h_blue = [0.6 0.7];
s_blue = [0.7 0.9];
v_blue = [0.6 0.7];

mask_blue = createMaskAndShow(board_trans_img,h_blue,s_blue,v_blue)

[blue_puck_img_coord,blue_puck_cell_coord,blue_puck_world_coord] = ...
    findColouredPuck(mask_blue,square_center_img,2,tform_img,tform_world,'blue');

hold on
plot(blue_puck_img_coord(:,1),blue_puck_img_coord(:,2),'b*','MarkerSize',30)


%%%% Detect Green Puck %%%%
h_green = [0.3 0.4];
s_green = [0.8 1];
v_green = [0.5 0.7];

mask_green = createMaskAndShow(board_trans_img,h_green,s_green,v_green);

[green_puck_img_coord,green_puck_cell_coord,green_puck_world_coord] = ...
    findColouredPuck(mask_green,square_center_img,2,tform_img,tform_world,'green');

hold on
plot(green_puck_img_coord(:,1),green_puck_img_coord(:,2),'g*','MarkerSize',30)

display(green_puck_cell_coord)

%%%%% Create Grid %%%%%
grid = [ones(1,10) ; [ones(5,1) zeros(5,8) ones(5,1)] ; ones(1,10)];

grid = putInGrid(grid,red_puck_cell_coord,1);
grid = putInGrid(grid,blue_puck_cell_coord,2);
grid = putInGrid(grid,green_puck_cell_coord,3);

% grid = putInGrid(grid,[5 7],1)
grid = putInGrid(grid,[1,5],4); % Goal Position


display(grid')

%% Move Robot %%%%%

% % TCP Host and Port settings
host = '127.0.0.1'; % THIS IP ADDRESS MUST BE USED FOR THE VIRTUAL BOX VM
% host = '192.168.230.128'; % THIS IP ADDRESS MUST BE USED FOR THE VMWARE
% host = '192.168.0.100'; % THIS IP ADDRESS MUST BE USED FOR THE REAL ROBOT
rtdeport = 30003;
% vacuumport = 63352;

% Calling the constructor of rtde to setup tcp connction
rtde = rtde(host,rtdeport);

home = [-588.5,-133, 371, 2.2214, -2.2214, 0.00];

rtde.movej(home);

% Calling the constructor of vacuum to setup tcp connction
% vacuum = vacuum(host,vacuumport);

% corners = [];

% for i = 1:length(board_corners_world)
%     corners = cat(1,corners,[board_corners_world(i,1), board_corners_world(i,2), 20, 2.2214, -2.2214, 0.00])
% end

% for i = 1:length(blue_puck_world_coord)
%     puck_blue = cat(1,corners,[board_corners_world(i,1), board_corners_world(i,2), 20, 2.2214, -2.2214, 0.00])
% end
% [[4,3,2,1]',corners]



% pose1 = rtde.movej(pts(1,:));
% pose2 = rtde.movej(pts(2,:));
% pose3 = rtde.movej(pts(3,:));
% pose4 = rtde.movej(pts(4,:));

% poses = [pose1;pose2;pose3;pose4];
% poses = [pose2;pose3;pose4];


% Pick up Green Puck
puck_green = [green_puck_world_coord(1), green_puck_world_coord(2), 20, 2.2214, -2.2214, 0.00];

pose1 = rtde.movej(puck_green);

% Move to bottom right

pt = square_center_world{5,8}

bottom_right = [pt(1), pt(2), 50, 2.2214, -2.2214, 0.00];

pt = square_center_world{1,8}

bottom_left = [pt(1), pt(2), 50, 2.2214, -2.2214, 0.00];

pt = square_center_world{5,1}

top_right = [pt(1), pt(2), 50, 2.2214, -2.2214, 0.00];

pt = square_center_world{1,1}

top_left = [pt(1), pt(2), 50, 2.2214, -2.2214, 0.00];

pose1 = rtde.movel(bottom_right);
pose2 = rtde.movel(bottom_left);
pose3 = rtde.movel(top_left);
pose4 = rtde.movel(top_right);

poses = [pose1;pose2;pose3;pose4];
% poses = [pose2;pose3];


rtde.drawPath(poses);
% yaxis([-0.9 0])
% xaxis([-0.4 0.1])

XMIN = -0.9;
XMAX = -0.1;
YMIN = -0.6;
YMAX = 0.2;

axis([XMIN XMAX YMIN YMAX])

%% Part B - BUG 2

% clear
% clc
% close all

startup_rvc;
% 
% grid = [
% 1 1 1 1 1 1 1 1 1 1;
% 1 3 1 4 0 0 0 0 0 1;
% 1 0 0 1 2 0 0 0 0 1;
% 1 0 0 0 0 0 2 0 0 1;
% 1 0 0 0 0 1 2 0 0 1;
% 1 0 0 0 0 0 1 0 0 1;
% 1 1 1 1 1 1 1 1 1 1];

for i = 1:(length(grid(1,:))-1)
    for j = 1:(length(grid(:,1))-1)
        if (grid(j,i) == 3)
            start = [j,i]
            grid(j,i) = 0;
        end

        if (grid(j,i) == 4)
            goal = [j,i]
            grid(j,i) = 0;
        end
    end
end

bug = myBug2(grid');
bug.plot();

bug_path = bug.query(start,goal,'animate')

axis xy
bug.plot
bug.plot_mline
hold on
plot(bug_path(:, 1), bug_path(:, 2),'LineWidth', 4, 'Color', 'g');

%%

function [puck_img_coord,puck_cell_coord,puck_world_coord] = findColouredPuck(mask,centers,width,tform_img,tform_world,colour)
    
    puck_img_coord = [];
    puck_cell_coord = [];
    puck_world_coord = [];
    
    for row = 1:length(centers(:,1))
        for col = 1:length(centers(1,:))
            
            point = centers{row,col};

            lower_x = point(2) - width;
            upper_x = point(2) + width;
            lower_y = point(1) - width;
            upper_y = point(1) + width;
    
            kernel = mask(lower_x:upper_x,lower_y:upper_y);
            if (sum(kernel(:)) == length(kernel(:)))
                puck_img_coord = cat(1,puck_img_coord,centers{row,col});
                puck_cell_coord = cat(1,puck_cell_coord,[row,col]);
                world_point = transformImgToWorld(tform_img,tform_world,point);
                puck_world_coord = cat(1,puck_world_coord,world_point);
            elseif (sum(kernel(:)) < length(kernel(:)) && sum(kernel(:)) > 0) % There is part of a puck in the image
                input(strcat('Please move ',colour,' into cell'));
            end
        end
    end

end

function grid = putInGrid(grid,puck_cell_coord,val)

if (length(puck_cell_coord(:,1)) == 1)
    cell_x = puck_cell_coord(1) + 1;
    cell_y = puck_cell_coord(2) + 1;
    
    grid(cell_x,cell_y) = val;
    return;
end

for i = 1:length(puck_cell_coord)
    cell_x = puck_cell_coord(i,1) + 1;
    cell_y = puck_cell_coord(i,2) + 1;
    
    grid(cell_x,cell_y) = val;
end

end

function mask = createMaskAndShow(im_hsv,h,s,v)
    
    mask =  (im_hsv(:,:,1) <= max(h))&(im_hsv(:,:,1) > min(h))&...
        (im_hsv(:,:,2) <= max(s))&(im_hsv(:,:,2) > min(s))&...
        (im_hsv(:,:,3) <= max(v))&(im_hsv(:,:,3) > min(v));
    se = strel('disk',7);
    mask = imclose(mask,se);
    mask = bwareaopen(mask,100);

    figure;
    imshow(mask)
end

function centers = findCenters(mask)
    
    % Find Circles and Return Centers
    blobs = regionprops(mask,'Centroid')
    
    centers = [blobs(1).Centroid ; blobs(2).Centroid ; blobs(3).Centroid ...
        ; blobs(4).Centroid]
    
    % Order Coordinates
    centers = sortrows(centers,1);
    centers = [sortrows(centers(1:2,:),2);sortrows(centers(3:4,:),2);]

end
