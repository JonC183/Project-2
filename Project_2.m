%% Part A

startup_rvc;

clear
clc
close all

% Load File

filename = '1.jpg'

blank_board = imread(filename);
blank_board_hsv = rgb2hsv(blank_board);

rgb = prism(6);
hsv = rgb2hsv(rgb);

figure;
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
figure;
imshow(mask_purple)

mask_yellow = imclose(mask_yellow,se);
mask_yellow = bwareaopen(mask_yellow,100);
figure;
idisp(mask_yellow)

[centers_purple,radii_purple] = imfindcircles(mask_purple,[5 15]);
figure(2)
hold on
plot(centers_purple(:,1),centers_purple(:,2),'*r');

X = centers_purple(:,1);
Y = centers_purple(:,2);

point1 = [0, 0];
point2 = [abs(900-250), 0];
point3 = [0, abs(525-75)];
point4 = [abs(900-250), abs(525-75)];

[centers_yellow,radii_yellow] = imfindcircles(mask_yellow,[3 10]);

figure(3)
hold on
plot(centers_yellow(:,1),centers_yellow(:,2),'*r')

% world = [point3; point4 ; point1; point2]
world = [point1; point3 ; point4; point2]

outputFrame = [abs(900-250) abs(525-75) 3];


tform = fitgeotrans([X Y],world,'projective');
figure;

board_trans_img = imwarp(blank_board,tform,OutputView =imref2d(outputFrame));
imshow(board_trans_img);


point1 = [-250, 75];
point2 = [-250, -525];
point3 = [-900, 75];
point4 = [-900, -525];

world = [point4; point1 ; point2; point3]
% outputFrame = [380 590 3]
outputFrame = [590 380 3]

tform = fitgeotrans([X Y],world,'projective');
figure;
board_trans_world = imwarp(blank_board,tform,OutputView =imref2d(outputFrame))
imshow(board_trans_world)

board_corners_img = centers_yellow
board_corners_world = transformPointsForward(tform,board_corners_img)

%% Move Robot


clear
close all
clc

% % TCP Host and Port settings
host = '127.0.0.1'; % THIS IP ADDRESS MUST BE USED FOR THE VIRTUAL BOX VM
%host = '192.168.230.128'; % THIS IP ADDRESS MUST BE USED FOR THE VMWARE
 %host = '192.168.0.100'; % THIS IP ADDRESS MUST BE USED FOR THE REAL ROBOT
port = 30003;
% 

% Calling the constructor of rtde to setup tcp connction
rtde = rtde(host,port);


pts = [board_corners_world(1,:), board_corners_world(2,:), 0, 2.2214, -2.2214, 0.00];


%% Part B - BUG 2

startup_rvc;

clear
clc
close all

grid = [
1 1 1 1 1 1 1 1 1 1;
1 3 1 4 0 0 0 0 0 1;
1 0 0 1 0 2 0 0 0 1;
1 0 0 0 0 2 0 0 0 1;
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

