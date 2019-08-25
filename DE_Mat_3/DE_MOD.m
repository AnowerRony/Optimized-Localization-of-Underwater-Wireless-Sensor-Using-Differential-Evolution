	clear all
	clc

	%Common Parameter Setting
	N=3; 		% Number of variables
	M=50; 		% Populations size 50
	F=0.5; 		% Mutation factor
	C=0.9;      % Crossover rate
	I_max=200; 	% Max iteration time
    Run=1;      % The number of test time
	X_max=[150,150 50];
	X_min=[-150,-150 0];

    Func=@OBJ_MOD;

	% 2.The whole test loop
    for r=1:Run
        iter=0;
        % 1.Generate MxN matrix
        for m=1:M
            for n=1:N
                X(m,n)=X_min(n)+rand()*(X_max(n)-X_min(n));
            end
            %fprintf('value of X:');
            %disp(X);
        end
		
    
        for i=1:I_max  % Stop when the iteration large than the max iteration time

            iter=iter+1;
            for m=1:M % For each individual
                % Mutation
                
				[V]=rand1(X,M,F,m);
                % Check if the element in the V matrix beyond the boundary.
                %fprintf('V matrix rand: ');
                %disp(V);
                for n=1:N
                    if V(1,n)>X_max(1,n)
                        V(1,n)=X_max(1,n);
                    end
                    if V(1,n)<X_min(1,n)
                        V(1,n)=X_min(1,n);
                    end
                    %fprintf('V matrix: ');
                    %disp(V);
                end
                

                % Crossover put the result in the U matrix
                jrand=floor(rand()*N+1);
                for n=1:N
                    R1=rand();
                    if (R1<C || n==jrand)
                        U(1,n)=V(1,n);
                    else
                        U(1,n)=X(m,n);
                    end
                    %fprintf('Offspring U: ');
                    %disp(U);
                end

                % Selection
                if Func(U(1,:)) < Func(X(m,:))
                    Tr=U(1,:);%U(1,:);
                    else
                    Tr=X(m,:);
                end
                % Use the selection result to replace the m row
                X(m,:)=Tr;
                
                % Evaluate each individual's fitness value, and put the result in the Y matrix.
                Y(m,:)=Func(X(m,:));
                Z(m,:)=X(m,:);
                
            end % Now the 1th individual generated

            % Select the lowest fitness value
            [y,ind1]=sort(Y,1);
			Y_min=y(1,:);
            
            [Ymin,ind] = min(Y);
            
            %for m=1:6
                %Z(ind,:)=X(ind,:);
            %end
            
            %fprintf('Distance difference: ');
            %disp(Y);
            
            fprintf('Least distance difference corresponding to 6 beacons: ');
            disp(Ymin);
            
            
            
            
            %fprintf('co-ordinate: ');
            %disp(Z);
            fprintf('fittest co-ordinates: ');
            disp(Z(ind,:));
            
            %w=(Z(m,:));
            W=mean(Z);
            
            %fprintf('Index: ');
            %disp(ind);
            
            fprintf('One of the most fit co-ordinate is: ');
            disp(W);
            
            z = mean(Ymin);
            fprintf('Error Rate is: ');
            disp(z);
            
            %disp(z);
          
            % plot the picture of iteration
            figure(1);
            plot(iter,z,'r+');
            xlabel('Iteration');
            ylabel('Error');
			title(sprintf('Iteration=%d, Error=%9.9f',i,z));
            grid on;
            hold on;
                        
        end % Finish I_max times iteration

    end % Run 30 times