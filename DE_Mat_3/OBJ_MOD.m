%clear all
	%clc
    
function y = OBJ_MOD (X)

% Input and parameters.
        %c = 0.707;
        %B = [150 90; -100 -120; -80 130; 140 -70; 60 120; -90 -130];
        
 % Beacon Nodes
    B1 = [150 90 50];
    B2 = [-100 -120 50];
    B3 = [-80 130 50];
    B4 = [140 -70 50];
    B5 = [60 120 50];
    B6 = [-90 -130 50];
    
  % Sensor Coordinate
    S = [80 40 0];
    
  %Euclidean distance from beacon nodes to known sensor node
    d1 = pdist2(B1,S,'euclidean'); 
    d2 = pdist2(B2,S,'euclidean');
    d3 = pdist2(B3,S,'euclidean');
    d4 = pdist2(B4,S,'euclidean');
    d5 = pdist2(B5,S,'euclidean');
    d6 = pdist2(B6,S,'euclidean');
    
  % Adding Gaussian Error and calculate avg
    %1st beacon
    for i=1:10
            d1 = d1 + erf(d1);
           d11(i)=d1;
    end    
    davg1 = mean(d11);
    
    %2nd beacon
    for i=1:10
            d2 = d2 + erf(d2);
           d12(i)=d2;
    end    
    davg2 = mean(d12);
   
    %3rd beacon
    for i=1:10
            d3 = d3 + erf(d3);
            
            d13(i)=d3;
    end   
    davg3 = mean(d13);
    
    %4th beacon
    for i=1:10
            d4 = d4 + erf(d4);
            
            d14(i)=d4;
    end   
    davg4 = mean(d14);
    
    %5th beacon
    for i=1:10
            d5 = d5 + erf(d5);
            
            d15(i)=d5;
    end   
    davg5 = mean(d15);
 
    %6th beacon
    for i=1:10
            d6 = d6 + erf(d6);
            
            d16(i)=d6;
    end    
    davg6 = mean(d16);

    %davg = (davg1+davg2+davg3+davg4+davg5+davg6)/6; 
    %disp(davg);
 
    %distance from beacon nodes to generated sensor nodes
    
    D1 = pdist2(X, B1);
    D2 = pdist2(X, B2);
    D3 = pdist2(X, B3);
    D4 = pdist2(X, B4);
    D5 = pdist2(X, B5);
    D6 = pdist2(X, B6);
    
    %D = (D1+D2+D3+D4+D5+D6)/6;
    
    %difference of known avg distance and calculated distance
    
    %y1 = abs(D1-davg1); 
    %y2 = abs(D2-davg2);
    %y3 = abs(D3-davg3);
    %y4 = abs(D4-davg4);
    %y5 = abs(D5-davg5);
    %y6 = abs(D6-davg6);
    
    y = [abs(D1-davg1), abs(D2-davg2), abs(D3-davg3), abs(D4-davg4), abs(D5-davg5), abs(D6-davg6)];
    
% Calculate average distance
            %davg = mean(vecnorm(B, 2, 2) / c);
            %disp(davg);
            
%distance from beacon nodes to sensor nodes
            %D = pdist2(X, B);
            %y = abs(D-davg);
            %difference of avg distance and calculated distance
            %for m=1:50
                %for n=1:6
                    %y = abs(y(1,:)-davg);
                %end
            %end

end