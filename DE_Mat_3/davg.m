clear all
clc
% Input and parameters.
c = 0.707;
B = [150 80; -100 -120; -80 130; 140 -70; 60 120; -90 -130];

% Calculate average.
%davg1 = mean(vecnorm(B, 2, 2) / c);



clear all
	clc


%for j = 1:6
    %for i = 1:2
        %m =m+((sqrt(B(j,i)^2+B(j,i)^2))/c);
    %end
%end
%davg = m/6;

% Input and parameters.
davg = 0;
m = 0;
c = 0.707;
B = [150 90; -100 -120; -80 130; 140 -70; 60 120; -90 -130];

% Get number of data points and dimensionality.
nPoints = size(B, 1);
nDim = size(B, 2);

% Iterate every data point.
for j = 1:nPoints

  % Calculate sum of squared elements in loop (for arbitrary dimensionality).
  temp = 0;
  for i = 1:nDim
        temp = temp + B(j, i).^2;
  end

  % Apply square root afterwards.  
  m = m + sqrt(temp) / c;

end

% Calculate average.
davg = m / nPoints;



%distance from beacon nodes to sensor nodes
            D = pdist2(Y, B);
            
            disp(D(1,1));
            
            % Input and parameters
            c = 0.707;

            % Calculate average distance
            davg = mean(vecnorm(B, 2, 2) / c);
            disp(davg);
            
            %difference of avg distance and calculated distance
            for m=1:50
                for n=1:6
                    diff = abs(D(m,n)-davg);
                    
                end
                fprintf('difference');
                    disp(diff);
            end
    

    