%% Part A

startup_rvc;

clear
clc
close all

% Load File

% filenames = strcat({'3','4','5'},'.jpg')
% filenames = strcat({'Pink_edges'},'.jpg')

% Load from Camera
% 
camList = webcamlist;
cam = webcam(1);
 
preview(cam)

%
images = {};
images_hsv = {};
% image_folder = 'webcamImages/';
image_folder = '';

% for i = 1:length(filenames)
%     images{i} = imread(filenames{i});
%     images_hsv{i} = rgb2hsv(images{i});
% 
%     webcam_filename = strcat(image_folder,string(i),'.mat')
%     save(webcam_filename,"images");
%     figure(1);
%     imshow(images_hsv{i})
% end

%%% Don't Start Until Ready
% prompt = 'Start';
% inputCommand = input(prompt)
% 
% i = 1;
for i = 1:30
    images{i} = snapshot(cam);
    images_hsv{i} = rgb2hsv(images{i});

    figure(1);
    imshow(images_hsv{i})
    % CODE FOR EXTRACTING ESSENTIAL OBJECTS HERE %%%
end

webcam_filename = strcat(image_folder,'test.mat')
save(webcam_filename,"images");


%% Grab Images Of Same Board


clear all
clc
close all

load('test.mat');

% Grab Images of the Index
img_idx = [5:30];

for i = 1:length(img_idx)
    images_hsv{i} = rgb2hsv(images{img_idx(i)}(:,200:1400,:));
    % images_hsv{i} = rgb2hsv(images{i+1});
    
    % webcam_filename = strcat(image_folder,string(i),'.mat')
    % save(webcam_filename,"images");

    figure(1);
    imshow(images_hsv{i})
    % CODE FOR EXTRACTING ESSENTIAL OBJECTS HERE %%%
end

figure(1);
imshow(images_hsv{end});

% Create Masks

% Mask for Purple Circles
h_purple = [0.76 0.80];
s_purple = [0.46 0.7];
v_purple = [0.3 0.6];

% Mask for Yellow Corners
h_yellow = [0.04 0.15];
s_yellow = [0.5 1];
v_yellow = [0.65 1];

% % Mask for Pink Corners
% h_pink = [0.9 0.95];
% s_pink = [0.67 0.88];
% v_pink = [0.8 1];

% Mask for Orange Corners
h_orange = [0 0.1];
s_orange = [0.5 0.75];
v_orange = [0.8 1];



rgb = prism(6);
hsv = rgb2hsv(rgb);

centers_purple = {};
centers_pink = {};
mask_purple = {};
mask_orange = {};

purple_mask_figure = 2;
orange_mask_figure = 3;


for idx = 1:length(images_hsv)
    % % Acquire a single image.
    % image = snapshot(cam);
    % image_hsv = rgb2hsv(image);
    image_hsv = images_hsv{idx};

    % Display the image.
    mask_purple{idx} = createMaskAndShow(image_hsv,h_purple,s_purple,v_purple, ...
        purple_mask_figure,'purple');
    centers_purple{idx} = findCenters(mask_purple{idx});
%     
%     mask_pink{idx} = createMaskAndShow(image_hsv,h_pink,s_pink,v_pink, ...
%         pink_mask_figure,'pink');
%     centers_pink{idx} = findCenters(mask_pink{idx});

    
    mask_pink{idx} = createMaskAndShow(image_hsv,h_orange,s_orange,v_orange, ...
        orange_mask_figure,'pink');
    centers_pink{idx} = findCenters(mask_pink{idx});
    
    % mask_yellow{idx} = createMaskAndShow(image_hsv,h_yellow,s_yellow,v_yellow, ...
    %     yellow_mask_figure);
    % centers_yellow{idx} = findCenters(mask_yellow{idx});
    
end

sum_centers_purple = zeros(4,2);
% sum_centers_yellow = zeros(4,2);
sum_centers_pink = zeros(4,2);

sum_centers_purple = zeros(4,2);
% sum_centers_yellow = zeros(4,2);
sum_centers_pink = zeros(4,2);

outliers = zeros(length(images_hsv),1);

% Remove Outliers
for j = 1:length(centers_purple)
    for pt = 1:4
        points{pt}(j,:) = centers_purple{j}(pt,:)
    end
end

% Bug with this method since if points are aronud same
% area it will still count
for pt = 1:4
    [outlier,L,U,C] = isoutlier(points{1,pt})
    display([points{1,pt} outlier])
    outlier = outlier(:,1) + outlier(:,2)
    outliers = outliers + outlier;
end

% Averages at start value
% average_centers_purple = centers_purple{1};
% average_centers_pink = centers_pink{1};

n_sum = 0;

for idx = 1:length(images_hsv)
    figure(2);
    hold on
    plot(centers_purple{idx}(:,1),centers_purple{idx}(:,2),'*r');
    % Get Sum of all Purple Centers
%     if (centers_purple{idx}(:,1) - centers_purple{idx - 1,:})
%         
%     end
    
    % Get Rid of Outliers
    if ((outliers(idx) == 0) && (length(centers_purple{idx}(:,1)) == 4) && (length(centers_pink{idx}(:,1)) == 4))
        sum_centers_purple(:,1) = sum_centers_purple(:,1) + centers_purple{idx}(:,1);
        sum_centers_purple(:,2) = sum_centers_purple(:,2) + centers_purple{idx}(:,2);
        n_sum = n_sum + 1;
    end
    % figure(3);
    % hold on
    % plot(centers_yellow{idx}(:,1),centers_yellow{idx}(:,2),'*r-');
    % % Get Sum of all Yellow Centers
    % sum_centers_yellow(:,1) = sum_centers_yellow(:,1) + centers_yellow{idx}(:,1);
    % sum_centers_yellow(:,2) = sum_centers_yellow(:,2) + centers_yellow{idx}(:,2);
    % 

    figure(3);
    imshow(mask_pink{idx});
    title('Pink');
    hold on
    plot(centers_pink{idx}(:,1),centers_pink{idx}(:,2),'*r-');
    % Get Sum of outliers
    if ((outliers(idx) == 0) && (length(centers_purple{idx}(:,1)) == 4) && (length(centers_pink{idx}(:,1)) == 4))
        sum_centers_pink(:,1) = sum_centers_pink(:,1) + centers_pink{idx}(:,1);
        sum_centers_pink(:,2) = sum_centers_pink(:,2) + centers_pink{idx}(:,2);
    end

    
end

% Get Averages
average_centers_purple = sum_centers_purple./n_sum;
% average_centers_yellow = sum_centers_yellow./length(images_hsv)
% average_centers_pink = sum_centers_pink./n_sum;
average_centers_pink = centers_pink{1,11};

% Plot Averages of Each
figure(2);
hold on
plot(average_centers_purple(:,1),average_centers_purple(:,2),'*g-');
% 
% figure(3);
% hold on
% plot(average_centers_yellow(:,1),average_centers_yellow(:,2),'*g-');

figure(3);
hold on
plot(average_centers_pink(:,1),average_centers_pink(:,2),'*g-');

% Conversion for world frame %%%%%

% variables_folder = 'H:\MTRN4320\GitHub\Project-2\savedVariables\';
variables_folder = '';

% Mask for Red Pucks
h_red = [0.95 0.99];
s_red = [0.6 0.98];
v_red = [0.55 1];

% Mask for Blue Pucks
h_blue = [0.6 0.7];
s_blue = [0.3 1];
v_blue = [0.55 1];

% Mask for Green Pucks
h_green = [0.3 0.4];
s_green = [0.8 1];
v_green = [0.4 0.7];

point1 = [-250, 75];
point2 = [-250, -525];
point3 = [-900, 75];
point4 = [-900, -525];

world = [point1 ; point3 ; point2 ; point4];

outputFrameWorld = [abs(900-250) abs(525+75)];

tform_world = fitgeotrans(average_centers_purple,world,'projective');
board_trans_world = imwarp(images{end},tform_world,'OutputView',imref2d(outputFrameWorld));

%%%% Conversion for Image %%%%
point1 = [0 0];
point2 = [380 0];
point3 = [0 590];
point4 = [380 590];

world_img = [point1 ; point3 ; point2 ; point4];

outputFrameImg = [590 380];

tform_img = fitgeotrans(average_centers_pink,world_img,'projective');

board_trans_img = imwarp(images_hsv{end},tform_img,'OutputView',imref2d(outputFrameImg));
board_trans_img_rgb = imwarp(images{end},tform_img,'OutputView',imref2d(outputFrameImg));

%%% SAVE TRANS BOARD IMAGE %%%
webcam_filename = strcat(variables_folder,'board_trans.mat')
save(webcam_filename,"board_trans_img","board_trans_img_rgb");

figure(4);
imshow(board_trans_img);

% Get Coordinates of Corners from Image
board_corners_img = average_centers_pink;

% Convert to World Coordinates
board_corners_world = transformPointsForward(tform_world,board_corners_img);

%%% SAVE BOARD CORNERS %%%
webcam_filename = strcat(variables_folder,'board_corners.mat')
save(webcam_filename,"board_corners_img","board_corners_world");

% Produce Transform for Transformed Image to World
% trans_T_world = trans_T_original*original_T_world
% trans_T_original = original_T_trans'
% tform_img_to_world = tform_img.T'*tform_world

% Test to see inverse of tform_img
P_img = [0 0];

P_original = transformPointsInverse(tform_img,P_img);
P_world = transformImgToWorld(tform_img,tform_world,P_img);

figure(1)
hold on
plot(P_original(1),P_original(2),'g*','MarkerSize',30)

num_rows = 5;
num_cols = 8;

cols_size = round(outputFrameImg(1)/num_cols);
rows_size = round(outputFrameImg(2)/num_rows);

square_center_img = {zeros(num_rows,num_cols)};
square_kernel_img = {};
square_center_world = {};

figure(4)
for row = 1:num_rows
    for col = 1:num_cols
        point = [(rows_size*row - round(rows_size/2))...
            (cols_size*(num_cols+1 - col) - round(cols_size/2))];
        % Reversed order along cols to all bottom left square to be cell
        % [1,1]

        square_center_img{row,col} = point;
        
        hold on
        % Plot the point
        plot(square_center_img{row,col}(1),square_center_img{row,col}(2),'r*');

        hold off
        square_center_world{row,col} = transformImgToWorld(tform_img,tform_world,point);
    end
end

%%% SAVE SQUARE CENTERS VARIABLES %%%
variables_filename = strcat(variables_folder,'Square_Centers.mat')
save(variables_filename,"square_center_img","square_center_world");


mask_red = createMaskAndShow(board_trans_img,h_red,s_red,v_red,5,'red');
mask_blue = createMaskAndShow(board_trans_img,h_blue,s_blue,v_blue,6,'blue');
mask_green = createMaskAndShow(board_trans_img,h_green,s_green,v_green,7,'green');

%%%% Detect Red Pucks %%%%
[red_puck_img_coord,red_puck_cell_coord,red_puck_world_coord] = ...
    findColouredPuck(mask_red,square_center_img,2,tform_img,tform_world,'red');

%%% SAVE PUCK VARIABLES %%%
variables_filename = strcat(variables_folder,'red_pucks.mat')
save(variables_filename,"red_puck_img_coord","red_puck_cell_coord","red_puck_world_coord");

figure(5)
hold on
plot(red_puck_img_coord(:,1),red_puck_img_coord(:,2),'r*','MarkerSize',30)

%%%% Detect Blue Pucks %%%%
[blue_puck_img_coord,blue_puck_cell_coord,blue_puck_world_coord] = ...
    findColouredPuck(mask_blue,square_center_img,2,tform_img,tform_world,'blue');

%%% SAVE PUCK VARIABLES %%%
variables_filename = strcat(variables_folder,'blue_pucks.mat')
save(variables_filename,"blue_puck_img_coord","blue_puck_cell_coord","blue_puck_world_coord");

figure(6)
hold on
plot(blue_puck_img_coord(:,1),blue_puck_img_coord(:,2),'b*','MarkerSize',30)

%%%% Detect Green Puck %%%%
[green_puck_img_coord,green_puck_cell_coord,green_puck_world_coord] = ...
    findColouredPuck(mask_green,square_center_img,2,tform_img,tform_world,'green');

figure(7)
hold on
plot(green_puck_img_coord(:,1),green_puck_img_coord(:,2),'g*','MarkerSize',30)

%%% SAVE PUCK VARIABLES %%%
variables_filename = strcat(variables_folder,'green_pucks.mat')
save(variables_filename,"green_puck_img_coord","green_puck_cell_coord","green_puck_world_coord");

display("Green puck at;")
display(green_puck_cell_coord)

%% Create Grid %%%%%
grid = [ones(1,10) ; [ones(5,1) zeros(5,8) ones(5,1)] ; ones(1,10)];

grid = putInGrid(grid,red_puck_cell_coord,1);
grid = putInGrid(grid,blue_puck_cell_coord,2);
grid = putInGrid(grid,green_puck_cell_coord,3);

% grid = putInGrid(grid,[5 7],1)

x_cell = 5;
y_cell = 6;
grid = putInGrid(grid,[x_cell,y_cell],4); % Goal Position

% flip grid 90 degrees
% grid = rot90(grid,-1)

% Plot Goal on Grid
figure(4);
hold on
plot(square_center_img{5,5}(1),square_center_img{5,5}(2),'r*','markersize',20)

%% Part B - BUG 2

% variables_folder = 'savedVariables\';
variables_folder = '';

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

grid = [
1 1 1 1 1 1 1 1 1 1
1 0 0 1 0 0 4 0 0 1
1 3 0 1 0 0 0 0 0 1
1 0 0 2 0 2 0 0 0 1
1 0 0 1 0 1 0 0 0 1
1 0 0 0 0 2 0 0 0 1
1 1 1 1 1 1 1 1 1 1];

% flip grid 90 degrees
% grid = rot90(grid,-1)

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

bug = Bug2(grid');

figure;
bug.plot();

bug_path = bug.query(start,goal,'animate')

axis xy
bug.plot
bug.plot_mline
hold on
plot(bug_path(:, 1), bug_path(:, 2),'LineWidth', 4, 'Color', 'g');

% Convert bug_path to real coordinates
bug_path(:,1) = bug_path(:,1) - 1
bug_path(:,2) = bug_path(:,2) - 1

%%% SAVE PUCK VARIABLES %%%
variables_filename = strcat(variables_folder,'bug2_path.mat')
save(variables_filename,"bug_path");

%% Move Robot %%%%%

clear all
clc
close all

% % TCP Host and Port settings
% host = '127.0.0.1'; % THIS IP ADDRESS MUST BE USED FOR THE VIRTUAL BOX VM
% host = '192.168.230.128'; % THIS IP ADDRESS MUST BE USED FOR THE VMWARE
host = '192.168.0.100'; % THIS IP ADDRESS MUST BE USED FOR THE REAL ROBOT
rtdeport = 30003;

vacuumport = 63352;

% Calling the constructor of rtde to setup tcp connction
rtde = rtde(host,rtdeport);

% Calling the constructor of vacuum to setup tcp connction
vacuum = vacuum(host,vacuumport);

% load('\\savedVariables\blue_pucks.mat')
% load('savedVariables\blue_pucks.mat');
% load('savedVariables\red_pucks.mat');
% load('savedVariables\green_pucks.mat');
% load('savedVariables\board_corners.mat');
% load('savedVariables\Square_Centers.mat');
% load('savedVariables\board_trans.mat');
% load('savedVariables\bug2_path.mat');

load('blue_pucks.mat');
load('red_pucks.mat');
load('green_pucks.mat');
load('board_corners.mat');
load('Square_Centers.mat');
load('board_trans.mat');
load('bug2_path.mat');

figure(1);
imshow(board_trans_img);
hold on

for i = 1:length(bug_path)
    plot(square_center_img{bug_path(i, 1)}(1),square_center_img{bug_path(i, 1)}(2),'g-');
end

home = [-588.5,-133, 371, 2.2214, -2.2214, 0.00];

puck_height = 6;
raised_height = 50;

lowered = [0,0,puck_height, 2.2214, -2.2214, 0.00];
raised = [0,0,raised_height, 2.2214, -2.2214, 0.00];

rtde.movej(home);



% pose1 = rtde.movej(pts(1,:));
% pose2 = rtde.movej(pts(2,:));
% pose3 = rtde.movej(pts(3,:));
% pose4 = rtde.movej(pts(4,:));

% poses = [pose1;pose2;pose3;pose4];
% poses = [pose2;pose3;pose4];

pt = square_center_world{2,4}

move_pt = [pt(1), pt(2), raised(3:end)];

Prompt = 'Going to Board Corners'
input(Prompt)

% Move to bottom right

pt = square_center_world{1,1}

bottom_right = [pt(1), pt(2), 50, 2.2214, -2.2214, 0.00];

pt = square_center_world{5,1}

bottom_left = [pt(1), pt(2), 50, 2.2214, -2.2214, 0.00];

pt = square_center_world{5,8}

top_right = [pt(1), pt(2), 50, 2.2214, -2.2214, 0.00];

pt = square_center_world{1,8}

top_left = [pt(1), pt(2), 50, 2.2214, -2.2214, 0.00];

pose1 = rtde.movel(bottom_right);
pose2 = rtde.movel(bottom_left);
pose3 = rtde.movel(top_left);
pose4 = rtde.movel(top_right);

% Move to blue obstacles

poses = [pose1;pose2;pose3;pose4];
% poses = [pose2;pose3];

figure(3)
rtde.drawPath(poses);
% yaxis([-0.9 0])
% xaxis([-0.4 0.1])

XMIN = -0.9;
XMAX = -0.1;
YMIN = -0.6;
YMAX = 0.2;

axis([XMIN XMAX YMIN YMAX])

%%%%% Move Green Puck Using Bug2 %%%%%
%%%% Pick up Green Puck %%%%
puck_green_lowered = [green_puck_world_coord(1), green_puck_world_coord(2), lowered(3:end)];
puck_green_raised = [green_puck_world_coord(1), green_puck_world_coord(2), raised(3:end)];

pose1 = rtde.movel(puck_green_raised);

Prompt = 'At Green Puck Raised'
input(Prompt)

pose2 = rtde.movel(puck_green_lowered);

vacuum.grip();
Prompt = 'At Green Puck Lowered'
input(Prompt)

pose3 = rtde.movel(puck_green_raised)

Prompt = 'At Green Puck Raised'
input(Prompt)

poses = [pose1;pose2;pose3]

figure(2);
rtde.drawPath(poses);
% yaxis([-0.9 0])
% xaxis([-0.4 0.1])

XMIN = -0.9;
XMAX = -0.1;
YMIN = -0.6;
YMAX = 0.2;

axis([XMIN XMAX YMIN YMAX])

poses = [];

for i = 1:length(bug_path)
    pt = square_center_world{bug_path(i,1),bug_path(i,2)}
    move_green = [pt(1), pt(2), raised(3:end)];
    rtde.movel(move_green)
end

% figure(1);
% rtde.drawPath(poses);
% yaxis([-0.9 0])
% xaxis([-0.4 0.1])

XMIN = -0.9;
XMAX = -0.1;
YMIN = -0.6;
YMAX = 0.2;

axis([XMIN XMAX YMIN YMAX])

final_green_lowered = [pt(1), pt(2), lowered(3:end)];

pose2 = rtde.movel(final_green_lowered);

vacuum.release();

pause(3)

final_green_raised = [pt(1), pt(2), raised(3:end)];
pose2 = rtde.movel(final_green_raised);

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

function mask = createMaskAndShow(im_hsv,h,s,v,idx,colour)
    
    mask =  (im_hsv(:,:,1) <= max(h))&(im_hsv(:,:,1) > min(h))&...
        (im_hsv(:,:,2) <= max(s))&(im_hsv(:,:,2) > min(s))&...
        (im_hsv(:,:,3) <= max(v))&(im_hsv(:,:,3) > min(v));
    se = strel('disk',3);
    mask = imclose(mask,se);
    mask = bwareaopen(mask,30);

    figure(idx);
    imshow(mask);
    title(colour);
end

function centers = findCenters(mask)
    
    % Find Circles and Return Centers
    blobs = regionprops(mask,'Centroid');
    
    centers = [blobs(1).Centroid ; blobs(2).Centroid ; blobs(3).Centroid ...
        ; blobs(4).Centroid];
    
    % Order Coordinates
    centers = sortrows(centers,1);
    centers = [sortrows(centers(1:2,:),2);sortrows(centers(3:4,:),2);];

end

function P_world = transformImgToWorld(tform_img,tform_world,P_img)
    P_original = transformPointsInverse(tform_img,P_img)
    P_world = transformPointsForward(tform_world,P_original)
end
