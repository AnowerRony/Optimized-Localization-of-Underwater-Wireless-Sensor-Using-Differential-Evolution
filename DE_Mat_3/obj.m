
function x =  (X)

% Input and parameters.
        %c = 0.707;
        %B = [150 90; -100 -120; -80 130; 140 -70; 60 120; -90 -130];
        
 % Beacon Nodes
    B1 = [ 150 90];
    
    
  % Sensor Coordinate
    S = [80 40];
    
  %Euclidean distance from beacon nodes to known sensor node
    d1 = pdist2(B1,S,'euclidean'); 
    
    
  % Adding Gaussian Error
    %1st beacon
    for i=1:10
            d1 = d1 + erf(d1);
           d11(i)=d1;
    end
    %disp(d11);
    
    davg1 = mean(d11);
    disp(davg1);
    
    %2nd beacon
    
 
    %distance from beacon nodes to generated sensor nodes
    
    D1 = pdist2(X, B1);
    
    
    %difference of known avg distance and calculated distance
    
    x = abs(D1-86.023); 
    
    


end