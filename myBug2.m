%BUG2 Bug navigation class
%
% A concrete subclass of the abstract Navigation class that implements the bug2 
% navigation algorithm.  This is a simple automaton that performs local 
% planning, that is, it can only sense the immediate presence of an obstacle.
%
% Methods::
%   Bug2        Constructor
%   query       Find a path from start to goal
%   plot        Display the obstacle map
%   display     Display state/parameters in human readable form
%   char        Convert to string
%
% Example::
%         load map1             % load the map
%         bug = Bug2(map);      % create navigation object
%         start = [20,10]; 
%         goal = [50,35];
%         bug.query(start, goal);   % animate path
%
% Reference::
% -  Dynamic path planning for a mobile automaton with limited information on the environment,,
%    V. Lumelsky and A. Stepanov, 
%    IEEE Transactions on Automatic Control, vol. 31, pp. 1058-1063, Nov. 1986.
% -  Robotics, Vision & Control, Sec 5.1.2,
%    Peter Corke, Springer, 2011.
%  
% See also Navigation, DXform, Dstar, PRM.



% Copyright (C) 1993-2017, by Peter I. Corke
%
% This file is part of The Robotics Toolbox for MATLAB (RTB).
% 
% RTB is free software: you can redistribute it and/or modify
% it under the terms of the GNU Lesser General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% RTB is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU Lesser General Public License for more details.
% 
% You should have received a copy of the GNU Leser General Public License
% along with RTB.  If not, see <http://www.gnu.org/licenses/>.
%
% http://www.petercorke.com

classdef myBug2 < Navigation

    properties(Access=protected)
        H       % hit points
        j       % number of hit points
        mline   % line from starting position to goal
        step    % state, in step 1 or step 2 of algorithm
        edge    % edge list
        k       % edge index
        heading % states heading of robot
    end

    methods

        function bug = myBug2(varargin)
            %Bug2.Bug2 Construct a Bug2 navigation object 
            %
            % B = Bug2(MAP, OPTIONS) is a bug2 navigation object, and MAP is an occupancy grid,
            % a representation of a planar world as a matrix whose elements are 0 (free
            % space) or 1 (occupied).
            %
            % Options::
            % 'goal',G      Specify the goal point (1x2)
            % 'inflate',K   Inflate all obstacles by K cells.
            %
            % See also Navigation.Navigation.

            % invoke the superclass constructor
            bug = bug@Navigation(varargin{:});

            bug.H = [];
            bug.j = 1;
            bug.step = 1;
            bug.heading = [];

        end

        function pp = query(bug, start, goal, varargin)
            %Bug2.query  Find a path
            %
            % B.query(START, GOAL, OPTIONS) is the path (Nx2) from START (1x2) to GOAL
            % (1x2).  Row are the coordinates of successive points along the path.  If
            % either START or GOAL is [] the grid map is displayed and the user is
            % prompted to select a point by clicking on the plot.
            %
            % Options::
            %  'animate'   show a simulation of the robot moving along the path
            %  'movie',M   create a movie
            %  'current'   show the current position position as a black circle
            %
            % Notes::
            % - START and GOAL are given as X,Y coordinates in the grid map, not as
            %   MATLAB row and column coordinates.
            % - START and GOAL are tested to ensure they lie in free space.
            % - The Bug2 algorithm is completely reactive so there is no planning
            %   method.
            % - If the bug does a lot of back tracking it's hard to see the current
            %   position, use the 'current' option.
            % - For the movie option if M contains an extension a movie file with that
            %   extension is created.  Otherwise a folder will be created containing
            %   individual frames.
            %
            % See also Animate.
         
            opt.animate = false;
            opt.movie = [];
            opt.current = false;
            
            opt = tb_optparse(opt, varargin);
            
            if ~isempty(opt.movie)
                anim = Animate(opt.movie);
                opt.animate = true;
            end
       
            % make sure start and goal are set and valid
            bug.start = []; bug.goal = [];
            bug.checkquery(start, goal);
            
            % compute the m-line
            %  create homogeneous representation of the line
            %  line*[x y 1]' = 0
            bug.mline = homline(bug.start(1), bug.start(2), ...
                bug.goal(1), bug.goal(2));
            bug.mline = bug.mline / norm(bug.mline(1:2));
            
            display(bug.mline)

            if opt.animate
                bug.plot();
                
                bug.plot_mline();
            end
            
            % iterate using the next() method until we reach the goal
            robot = bug.start(:);
            bug.step = 1;
            path = bug.start(:);

            dirs = [1 0; 0 1; -1 0; 0 -1; -1 1; -1 -1; 1 -1; 1 1];

            while true
                if opt.animate
                    plot(robot(1), robot(2), 'g.', 'MarkerSize', 12);
                    if opt.current
                        h = plot(robot(1), robot(2), 'ko', 'MarkerSize', 8);
                    end
                    drawnow
                    if ~isempty(opt.movie)
                        anim.add();
                    end
                    if opt.current
                        delete(h)
                    end
                end

                % move to next point on path
                robot = bug.next(robot,dirs);

                % are we there yet?
                if isempty(robot)
                    % yes, exit the loop
                    break
                else
                    % no, append it to the path
                    path = [path robot(:)];
                    % [path ; bug.heading']
                end
            end
            
            if ~isempty(opt.movie)
                anim.close();
            end

            % only return the path if required
            if nargout > 0
                pp = path';
            end
        end
        
        function plot_mline(bug, ls)
            
                % parameters of the M-line, direct from initial position to goal
                % as a vector mline, such that [robot 1]*mline = 0
                
                if nargin < 2
                    ls = 'k--';
                end
                dims = axis;
                xmin = dims(1); xmax = dims(2);
                ymin = dims(3); ymax = dims(4);
                
                hold on
                if bug.mline(2) == 0
                    % handle the case that the line is vertical
                    plot([bug.start(1) bug.start(1)], [ymin ymax], 'k--');
                else
                    x = [xmin xmax]';
                    y = -[x [1;1]] * [bug.mline(1); bug.mline(3)] / bug.mline(2);
                    plot(x, y, ls);
                end
        end
        
        function n = next(bug, robot, dirs)
            
            % implement the main state machine for bug2
            n = [];
            robot = robot(:);
            % these are coordinates (x,y)
          
            if bug.step == 1
                % Step 1.  Move along the M-line toward the goal

                if colnorm(bug.goal - robot) == 0 % are we there yet?
                    return
                end

                % motion on line toward goal
                d = bug.goal-robot;
                if abs(d(1)) > abs(d(2))
                    % line slope less than 45 deg
                    dx = sign(d(1));
                    L = bug.mline;
                    y = -( (robot(1)+dx)*L(1) + L(3) ) / L(2);
                    % dy = round(y - robot(2));
                    dy = 0;
                else
                    % line slope greater than 45 deg
                    dy = sign(d(2));
                    L = bug.mline;
                    x = -( (robot(2)+dy)*L(2) + L(3) ) / L(1);
                    % dx = round(x - robot(1));
                    dx = 0;
                end
                

                % detect if next step is an obstacle
                if bug.isoccupied(robot + [dx; dy])
                    bug.message('(%d,%d) obstacle!', n);
                    bug.H(bug.j,:) = robot; % define hit point
                    bug.step = 2;
                    bug.heading = [round(dx); round(dy)];
                else
                    n = robot + [dx; dy];
                    bug.heading = [round(dx); round(dy)];
                end
            end % step 1

            if bug.step == 2
                % Step 2.  Move around the obstacle until we reach a point
                % on the M-line closer than when we started.
                if colnorm(bug.goal-robot) == 0 % are we there yet?
                    return
                end

                E = 1;
                S = 2;
                W = 3;
                N = 4;
                % SW = 5;
                % NW = 6;
                % NE = 7;
                % SE = 8;
                
                % Convert NW,SW,SE,NE to N,E,S,W
                % Only applies if heading is diagonal which will only 
                % take place during movement along mline
                % If this is the case always goes back to moving E
                if ~(isequal(bug.heading, dirs(N,:)') || isequal(bug.heading, dirs(E,:)') || ...
                    isequal(bug.heading, dirs(S,:)')|| isequal(bug.heading, dirs(W,:)'))
                    bug.heading = dirs(E,:)';
                end
                
                leftWallDir = bug.leftWallDirection(dirs);
                frontWallDir = bug.heading;

                leftWallCoord = robot + leftWallDir; 
                frontWallCoord = robot + frontWallDir;
    
                if bug.detectWall(robot,leftWallCoord) 
                    if bug.detectWall(robot,frontWallCoord) %%% There is a wall in front and to the left so turn right
                        bug.turnHeadingRight(dirs);
                    end
                    %%% Else stay on path and keep heading the same
                elseif ~bug.detectWall(robot,leftWallCoord)  
                        bug.turnHeadingLeft(dirs);            
                end 
                
                % Test new position to see if there is an obstacle
                new_position = robot + bug.heading;

                if bug.occgridnav(new_position(2),new_position(1)) == 0 
                    n = new_position;
                else % Robot stays put until finds a direction where it can move
                    n = robot;
                end                                                                                      

                % are we on the M-line now ?
                if abs( [n' 1]*bug.mline') <= 0.5
                    bug.message('(%d,%d) moving along the M-line', n);
                    % back to moving along the M-line
                    bug.j = bug.j + 1;
                    bug.step = 1;
                    return;
                end



                % no, keep going around
                bug.message('(%d,%d) keep moving around obstacle', n)
                bug.k = bug.k+1;
            end % step 2
        end % next
        
        function plan(bug)
            error('RTB:Bug2:badcall', 'This class has no plan method');
        end

        function wallDetected = detectWall(bug,robot,wallCoord)
            im = bug.occgridnav;

            if im(wallCoord(2),wallCoord(1)) == 1 % Wall Detect!! :)))
                wallDetected = 1;
                return;
            elseif im(wallCoord(2),wallCoord(1)) == 0 % Wall not Detected :((
                wallDetected = 0;
                return;
            else    % Your at edge of map
                wallDetected = -1;
                return;
            end
        end

        function turnHeadingLeft(bug,dirs)
            
            % Define N,E,W,S
            E = 1;
            S = 2;
            W = 3;
            N = 4;

            if isequal(bug.heading, dirs(N,:)')
                bug.heading = dirs(W,:)';
            elseif isequal(bug.heading, dirs(E,:)')
                bug.heading  = dirs(N,:)';
            elseif isequal(bug.heading, dirs(S,:)')
                bug.heading = dirs(E,:)';
            elseif isequal(bug.heading, dirs(W,:)')
                bug.heading  = dirs(S,:)';
            end
        end

        function turnHeadingRight(bug,dirs)
            
            % Define N,E,W,S
            E = 1;
            S = 2;
            W = 3;
            N = 4;

            if isequal(bug.heading, dirs(N,:)')
                bug.heading = dirs(E,:)';
            elseif isequal(bug.heading, dirs(E,:)')
                bug.heading  = dirs(S,:)';
            elseif isequal(bug.heading, dirs(S,:)')
                bug.heading = dirs(W,:)';
            elseif isequal(bug.heading, dirs(W,:)')
                bug.heading  = dirs(N,:)';
            end
        end

        function leftWallDir = leftWallDirection(bug,dirs)
            % Define N,E,W,S
            E = 1;
            S = 2;
            W = 3;
            N = 4;

            if isequal(bug.heading, dirs(N,:)')
                leftWallDir = dirs(W,:)';
            elseif isequal(bug.heading, dirs(E,:)')
                leftWallDir  = dirs(N,:)';
            elseif isequal(bug.heading, dirs(S,:)')
                leftWallDir = dirs(E,:)';
            elseif isequal(bug.heading, dirs(W,:)')
                leftWallDir  = dirs(S,:)';
            end
        end


    end % methods
end % classdef
