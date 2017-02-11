%Kristin_Diep_withPressure
function separation = withPressure(fishPos, v)

global numFish friendRange pressure 
sqr = @(x) x .* x;
distance = @(a, b, c, d) sqrt(sqr(a - b) + sqr(c - d));
separation = zeros(2,numFish);

%Calculates separation for each fish with pressure
for fish1 = 1:numFish %Fish that is trying to find friends
    for fish2 = 1:numFish %All other fish
        if fish2 ~= fish1 %Makes sure it isn't comparing itself
            %If fish are moving in the positive x direction
            if(v > 0)
                %fish friend is to the left of main fish
                if(fishPos(2,fish2) > fishPos(2,fish1))
                    %Within friend range and checks if right
                    %sensitivity of main fish is greater than left
                    %sensitivity of friend fish
                    if abs(distance(fishPos(1,fish2),fishPos(1,fish1),...
                            fishPos(2,fish2),fishPos(2,fish1))) < friendRange...
                            && pressure(1,fish1) > pressure(2,fish2)
                        
                        %if friend fish is too close, then get
                        %difference and add to y position of fish
                        %friend
                        if sqrt(sqr(abs(fishPos(2,fish2) - fishPos(2,fish1))))...
                                < pressure(1,fish1)
                            
                            separation(:,fish1) = (fishPos(:,fish2) -...
                                fishPos(:,fish1)) + fishPos(:,fish2);
                        end
                        
                    %Within friend range and checks if right
                    %sensitivity of fish friend is greater than left
                    %sensitivity of main fish
                    elseif abs(distance(fishPos(1,fish2),fishPos(1,fish1),...
                            fishPos(2,fish2),fishPos(2,fish1))) < friendRange...
                            && pressure(1,fish2) > pressure(2,fish1)
                        
                        %if main fish is too close, then get
                        %difference and add to y position of main
                        %fish
                        if sqrt(sqr(abs((fishPos(2,fish1) - fishPos(2,fish2)))))...
                                < pressure(1,fish2)
                            
                            separation(:,fish1) = (fishPos(:,fish2) -...
                                fishPos(:,fish1)) + fishPos(:,fish1);
                        end
                    end
                    
                %fish friend is to the right of main fish
                elseif(fishPos(2,fish2) < fishPos(2,fish1))
                    %Within friend range and checks if right
                    %sensitivity of main fish is greater than left
                    %sensitivity of friend fish
                    if abs(distance(fishPos(1,fish2),fishPos(1,fish1),...
                            fishPos(2,fish2),fishPos(2,fish1))) < friendRange...
                            && pressure(1,fish1) > pressure(2,fish2)
                        
                        %if friend fish is too close, then get
                        %difference and add to y position of fish
                        %friend
                        if sqrt(sqr(abs((fishPos(2,fish1) - fishPos(2,fish2)))))...
                                < pressure(1,fish1)
                            
                            separation(:,fish1) = (fishPos(:,fish2) -...
                                fishPos(:,fish1)) + fishPos(:,fish2);
                        end
                        
                    %Within friend range and checks if right
                    %sensitivity of fish friend is greater than left
                    %sensitivity of main fish
                    elseif abs(distance(fishPos(1,fish2),fishPos(1,fish1),...
                            fishPos(2,fish2),fishPos(2,fish1))) < friendRange...
                            && pressure(1,fish2) > pressure(2,fish1)
                        
                        %if main fish is too close, then get
                        %difference and add to y position of main
                        %fish
                        if sqrt(sqr(abs((fishPos(2,fish1) - fishPos(2,fish2)))))...
                                 < pressure(2,fish2)
                            
                            separation(:,fish1) = (fishPos(2,fish2) -...
                                fishPos(2,fish1)) + fishPos(2,fish1);
                        end
                    end
                end
            end
            
        elseif v < 0
            %fish friend is to the right of main fish
            if(fishPos(2,fish2) > fishPos(2,fish1))
                %Within friend range and checks if right
                %sensitivity of main fish is greater than left
                %sensitivity of friend fish
                if abs(distance(fishPos(1,fish2),fishPos(1,fish1),...
                            fishPos(2,fish2),fishPos(2,fish1))) < friendRange...
                            && pressure(1,fish1) > pressure(2,fish2)
                        
                        %if friend fish is too close, then get
                        %difference and add to y position of fish
                        %friend
                        if sqrt(sqr(abs((fishPos(2,fish1) - fishPos(2,fish2)))))...
                                < pressure(1,fish1)
                            
                            separation(:,fish1) = (fishPos(2,fish2) -...
                                fishPos(2,fish1)) + fishPos(2,fish2);
                        end
                        
                        %Within friend range and checks if left
                        %sensitivity of fish friend is greater than right
                        %sensitivity of main fish
                    elseif abs(distance(fishPos(1,fish2),fishPos(1,fish1),...
                            fishPos(2,fish2),fishPos(2,fish1))) < friendRange...
                            && pressure(2,fish2) > pressure(1,fish1)
                        
                        %if main fish is too close, then get
                        %difference and add to y position of main
                        %fish
                        if sqrt(sqr(abs((fishPos(2,fish1) - fishPos(2,fish2)))))...
                                < pressure(2,fish2)
                            
                            separation(:,fish1) = (fishPos(2,fish2) -...
                                fishPos(2,fish1)) + fishPos(2,fish1);
                        end
                end
                %fish friend is to the left of main fish
                if(fishPos(2,fish1) > fishPos(2,fish2))
                    %Within friend range and checks if left
                    %sensitivity of main fish is greater than right
                    %sensitivity of friend fish
                    if abs(distance(fishPos(1,fish2),fishPos(1,fish1),...
                            fishPos(2,fish2),fishPos(2,fish1))) < friendRange...
                            && pressure(2,fish1) > pressure(1,fish2)
                        
                        %if friend fish is too close, then get
                        %difference and add to y position of fish
                        %friend
                        if sqrt(sqr(abs(fishPos(2,fish2) - fishPos(2,fish1))))...
                                < pressure(2,fish1)
                            
                            separation(:,fish1) = (fishPos(2,fish2) -...
                                fishPos(2,fish1)) + fishPos(2,fish2);
                        end
                        
                        %Within friend range and checks if right
                        %sensitivity of fish friend is greater than left
                        %sensitivity of main fish
                    elseif abs(distance(fishPos(1,fish2),fishPos(1,fish1),...
                            fishPos(2,fish2),fishPos(2,fish1))) < friendRange...
                            && pressure(1,fish2) > pressure(2,fish1)
                        
                        %if main fish is too close, then get
                        %difference and add to y position of main
                        %fish
                        if sqrt(sqr(abs((fishPos(2,fish1) - fishPos(2,fish2)))))...
                                < pressure(1,fish2)
                            
                            separation(:,fish1) = (fishPos(2,fish2) -...
                                fishPos(2,fish1)) + fishPos(2,fish1);
                        end
                    end
                end
            end
        end
    end
end
end