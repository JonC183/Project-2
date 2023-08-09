%% Part A

startup_rvc;

clear
clc
close all

% Load File

filename = '4.jpg'

blank_board = imread(filename);
blank_board_hsv = rgb2hsv(blank_board);

rgb = prism(6);
hsv = rgb2hsv(rgb);

figure(1);
imshow(blank_board_hsv);

mask_yellow =  (blank_board_hsv(:,:,1) <= 0.1)&(blank_board_hsv(:,:,1) > 0.07)&...
        (blank_board_hsv(:,:,2) <= 0.65)&(blank_board_hsv(:,:,2) > 0.50)&...
        (blank_board_hsv(:,:,3) <= 0.9)&(blank_board_hsv(:,:,3) > 0.7);

mask_purple =  (blank_board_hsv(:,:,1) <= 0.79)&(blank_board_hsv(:,:,1) > 0.75)&...
        (blank_board_hsv(:,:,2) <= 0.58)&(blank_board_hsv(:,:,2) > 0.54)&...
        (blank_board_hsv(:,:,3) <= 0.6)&(blank_board_hsv(:,:,3) > 0.4);

se = strel('disk',7,4);
mask_purple = imclose(mask_purple,se);
mask_purple = bwareaopen(mask_purple,100);
figure(2);
imshow(mask_purple)

% Find Circles and Return Centers
[centers_purple,radii_purple] = imfindcircles(mask_purple,[5 15]);
hold on
plot(centers_purple(:,1),centers_purple(:,2),'*r');

X = centers_purple(:,1);
Y = centers_purple(:,2);

% Conversion for world frame
point1 = [-250, 75];
point2 = [-250, -525];
point3 = [-900, 75];
point4 = [-900, -525];

% world = [point3; point4 ; point1; point2]
world = [point4; point3 ; point2; point1]

outputFrame = [abs(900-250) abs(525-75) 3];

tform = fitgeotrans([X Y],world,'projective');

board_trans_world = imwarp(blank_board,tform,'OutputView',imref2d(outputFrame));

% Mask for Yellow Corners
mask_yellow = imclose(mask_yellow,se);
mask_yellow = bwareaopen(mask_yellow,100);
figure(3);
imshow(mask_yellow)

% Find Points for Yellow Corners
% [centers_yellow,radii_yellow] = imfindcircles(mask_yellow,[3 10]);
blobs = regionprops(mask_yellow,'Centroid')

centers_yellow = [blobs(1).Centroid ; blobs(2).Centroid ; blobs(3).Centroid ...
    ; blobs(4).Centroid]

hold on
plot(centers_yellow(:,1),centers_yellow(:,2),'*r')

% Display as Image
point1 = [0, 0];
point2 = [abs(900-250), 0];
point3 = [0, abs(525-75)];
point4 = [abs(900-250), abs(525-75)];

% world_img = [point3; point4 ; point1; point2]
world_img = [0 0 ; 0 590 ; 380 0 ; 380 590];
% world_img = [point4; point1 ; point2; point3]

% outputFrameImg = [abs(900-250) abs(525-75)  3];
outputFrameImg = [590 380];

tform_img = fitgeotrans(centers_yellow,world_img,'projective');

board_trans_img = imwarp(blank_board,tform_img,'OutputView',imref2d(outputFrameImg));
figure(1);
imshow(board_trans_img);

% Get Coordinates of Corners from Image
board_corners_img = centers_yellow

% Convert to World Coordinates
board_corners_world = transformPointsForward(tform,board_corners_img)

%%%%% Move Robot %%%%%

% % TCP Host and Port settings
host = '127.0.0.1'; % THIS IP ADDRESS MUST BE USED FOR THE VIRTUAL BOX VM
% host = '192.168.230.128'; % THIS IP ADDRESS MUST BE USED FOR THE VMWARE
% host = '192.168.0.100'; % THIS IP ADDRESS MUST BE USED FOR THE REAL ROBOT
port = 30003;
% 

% Calling the constructor of rtde to setup tcp connction
rtde = rtde(host,port);

pts = [];

for i = 1:length(board_corners_world)
    pts = cat(1,pts,[board_corners_world(i,1), board_corners_world(i,2), 20, 2.2214, -2.2214, 0.00])
end

[[4,3,2,1]',pts]

home = [-588.53, -133.30, 371.91, 2.2214, -2.2214, 0.00];


% pose1 = rtde.movej(pts(1,:));
pose2 = rtde.movej(pts(2,:));
pose3 = rtde.movej(pts(3,:));
pose4 = rtde.movej(pts(4,:));

% poses = [pose1;pose2;pose3;pose4];
poses = [pose2;pose3;pose4];

rtde.drawPath(poses);

%% Part B - BUG 2

clear
clc
close all

startup_rvc;

grid = [
1 1 1 1 1 1 1 1 1 1;
1 3 1 4 0 0 0 0 0 1;
1 0 0 1 2 2 0 0 0 1;
1 0 0 0 0 0 0 0 0 1;
1 0 0 0 0 1 0 0 0 1;
1 0 0 0 0 2 1 0 0 1;
1 1 1 1 1 1 1 1 1 1];

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
bug.plot();

bug_path = bug.query(start,goal,'animate')

axis xy
bug.plot
bug.plot_mline
hold on
plot(bug_path(:, 1), bug_path(:, 2),'LineWidth', 4, 'Color', 'g');

