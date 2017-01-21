%Assumptions:
%   -There are no objects/obstacles in the way of the fish(for the
%   beginning simulation) 
%   -All fish follow three rules: cohesion rule (fish moves toward the mean
%   position of its "closest neighbors"), separation rule (fish does not get
%   closer than some minimum distance to any neighbor, alignment rule(fish
%   heads in the mean direction to which its "closest neighbor heads"
%   -Separation rule takes priority over these rules 
%   -Fish are not allowed to collide with each other 
%   -Birthing fish and death of fish are not considered
%   -All fish have the same mass
%   -Assume no force is acting on the fish (perfect vaccum)
%Agent:only fish
%Modeling position of fish schooling
%Randomly generate fish positions 

%Setting width and height of environment 
upperHeight = 50;
lowerHeight = 0;
upperWidth = 50;
lowerWidth = 0; 

%number of fish wanted
numFish = 20;

%separation radius of fish
r = 2; 

%range of "closest friends" for cohesion
friendRange = 10;

%Time to wait between each frame(only if you don't want to press a key
%to continue)
wait = 0.2;
%Array to hold x and y position of fish
fishPos = zeros(2,numFish);
%Array to hold velocity in the x direction and velocity in the y direction
v = ones(2,numFish);
%The maximum velocity a fish is allowed 
max = 1;
%The minimum velocity a fish is allowed 
min = -1;
%Array to store maximum velocity 
vMax = max * ones(2,numFish);
%Array to store minimum velocity 
vMin = min * ones(2,numFish);
%The percentage to move fish closer to center of mass of its friends 
percentageToMoveFish = 100;
%Initializes all numFish fish 
for i = 1:numFish
    %Generates random starting velocity for fish that will never be above
    %velocity max 
    v(1,i) = rand(1); 
    v(2,i) = rand(1);
    display(v);
    %Initial random x position
    x(1) = randi([lowerWidth + r, upperWidth - r]); 
    %Initial random y position 
    y(1) = randi([lowerHeight + r, upperHeight - r]);
    
    %Create numFish fish using struct
    fish(i) = struct('ID', i, 'xPos', x(1), 'yPos', y(1), 'velocity', v(1));
    
    %Array of x and y position of each fish
    fishPos(1,i) = fish(i).xPos;
    fishPos(2,i) = fish(i).yPos;
    
end
%Anonymous functions to compute sqaure of a number and distance 
sqr = @(x) x .* x;
distance = @(a, b, c, d) sqrt(sqr(a - b) + sqr(c - d));
%Array to store separation values needed for each fish
separation = zeros(2,numFish);
%Array to store cohesion values needed for each fish
cohesion = zeros(2,numFish);
%Array to store alignment values neede for each fish
align = zeros(2,numFish);

%timestep
dt = 1;
%simulation length
simLength = 200; 
%number of iterations
numIterations = simLength/dt; 

for loop = 1:numIterations
    
    %Makes sure the velocity of the fish is not greater than the max
    %velocity and that the velocity is not less than the min velocity
    for eachRow = 1:2
        for eachFish = 1:numFish
            if v(eachRow,eachFish) > vMax(eachRow,eachFish)
                v(eachRow,eachFish) = vMax(eachRow,eachFish);
            elseif v(eachRow,eachFish) < vMin(eachRow,eachFish)
                v(eachRow,eachFish) = vMin(eachRow,eachFish);
            end
        end
    end
    
    %Separation
    for fish1 = 1:numFish %fish that is trying to find fish friends
        for fish2 = 1:numFish %all other fish
            if fish2 ~= fish1 %makes sure that it is not itself
                if(abs(distance(fishPos(1,fish2),fishPos(1,fish1),...
                        fishPos(2,fish2),fishPos(2,fish1))) < r)
                    %If the distance between fish1 and fish2 is less than 
                    %the separation radius then subtract the position of 
                    %the two fish and add to separation array so it can 
                    %later be added to velocity which will 
                    separation(1,fish1) = separation(1,fish1) +... 
                    abs(fishPos(1,fish1) - fishPos(1,fish2));
                
                    separation(2,fish1) = separation(2,fish1) +...
                        abs(fishPos(2,fish1) - fishPos(2,fish2));
                end
            end
        end
    end
    
    %Cohesion
    %Counter for number of friends
    numFriends = 0;
    for fish1 = 1:numFish %Fish that is trying to find friends
        friendCenter = zeros(2,fish1); %Reinitialize the array
        for fish2 = 1:numFish %Fish within friendRange
            if fish2 ~= fish1 %If it is not comparing itself
                if(abs(distance(fishPos(1,fish2),fishPos(1,fish1),...
                        fishPos(2,fish2),fishPos(2,fish1))) <= friendRange)
                    %Store x and y position of friends in array
                    friendCenter(1,fish1) = friendCenter(1,fish1)...
                        + fishPos(1,fish2);
                    friendCenter(2,fish1) = friendCenter(2,fish1)...
                        + fishPos(2,fish2);
                    %Increment counter
                    numFriends = numFriends + 1;
                else
                    %If not a friend, set to 0
                    friendCenter(:,fish1) = 0;
                end
            end
        end
        if numFriends > 0
            %Get mean position of friends
            friendCenter(1,fish1) = abs(sum(friendCenter(1,fish1)))/numFriends;
            friendCenter(2,fish1) = abs(sum(friendCenter(2,fish1)))/numFriends;
            
            %Moves fish closer to mean position of friends
            cohesion(1,fish1) = (friendCenter(1,fish1)...
                - fishPos(1,fish2))/percentageToMoveFish;
            cohesion(2,fish1) = (friendCenter(2,fish1)...
                - fishPos(2,fish2))/percentageToMoveFish;
        end
    end
    
    %Alignment 
    %Counter for fish friends
    numFriends = 0;
    for fish1 = 1:numFish %fish that is trying to find friends
        friendVelocity = zeros(2,fish1); %Reinitialize array
        for fish2 = 1:numFish %All other fish
            if fish2 ~= fish1 %Makes sure that it is not itslef
                if(abs(distance(fishPos(1,fish2),fishPos(1,fish1),...
                        fishPos(2,fish2),fishPos(2,fish1))) <= friendRange)
                    %Add x and y velocity of friends to array
                    friendVelocity(1,fish1) = friendVelocity(1,fish1)...
                        + v(1,fish2);
                    friendVelocity(2,fish1) = friendVelocity(2,fish1)...
                        + v(2,fish2);
                    %Increment counter
                    numFriends = numFriends + 1;
                else
                    %If not a friend, set to 0
                    friendVelocity(:,fish1) = 0;
                end
            end
        end
        if numFriends > 0
            %Get mean velocity of fish friends 
            friendVelocity(1,fish1) = sum(friendVelocity(1,fish1))/numFriends;
            friendVelocity(2,fish1) = sum(friendVelocity(2,fish1))/numFriends;
            %Move fish velocity closer to mean velocity 
            align(1,fish1) = friendVelocity(1,fish1) - v(1,fish2);
            align(2,fish1) = friendVelocity(2,fish1) - v(2,fish2);
        end
    end
    
    %Update velocity
    v = v + separation + cohesion + align; 
    
    %Rechecks to make sure velocity is within velocity bounds 
    for eachRow = 1:2
        for eachFish = 1:numFish
            if v(eachRow,eachFish) > vMax(eachRow,eachFish)
                
                v(eachRow,eachFish) = vMax(eachRow,eachFish);
            
            elseif v(eachRow,eachFish) < vMin(eachRow,eachFish)
                
                v(eachRow,eachFish) = vMin(eachRow,eachFish);
                
            end
        end
    end
    
    %Update position of fish
    fishPos = fishPos + v;
    
    %Wraps fish around so it is periodic
    for pos = 1:2
        for eachFish = 1:numFish
            if fishPos(1,eachFish) >= upperHeight
                
                fishPos(1,eachFish) = lowerHeight + r;
                
            elseif fishPos(1,eachFish) <= lowerHeight
                
                fishPos(1,eachFish) = upperHeight - r;
                
            elseif fishPos(2,eachFish) >= upperWidth
                
                fishPos(2,eachFish) = lowerWidth + r;
                
            elseif fishPos(2,eachFish) <= lowerWidth
                
                fishPos(2,eachFish) = upperWidth - r;
                
            end
        end
    end
    %Getting the updated x and y positions of numFish fish
    fishPosX{loop} = fishPos(1,:);
    fishPosY{loop} = fishPos(2,:);         

    %Graphing the school of fish
     clf;
     string = ['Iteration: ', num2str(loop)];
     title(string)
     xlabel('Width of Ocean Area')
     ylabel('Height of Ocean Area')
     hold on;

    scatter(fishPosX{loop}, fishPosY{loop});
    axis([lowerHeight upperHeight lowerWidth upperWidth]);
    % wait to go on to next image
    %fprintf('Waiting for any key to be pressed\n');
    %w = waitforbuttonpress;
    pause(wait);
end









            
