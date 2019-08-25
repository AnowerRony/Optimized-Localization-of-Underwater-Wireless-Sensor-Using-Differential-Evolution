	clear all
	clc

	%Common Parameter Setting
	N=2; 		% Number of variables
	M=10; 		% Populations size 50
	F=0.5; 		% Mutation factor
	C=0.9;       %0.9; 		% Crossover rate
	I_max=20; 	% Max iteration time
    Run=1;      % The number of test time
	X_max=[100,100];%X_max=[5.12,5.12];
	X_min=[-100,-100];

    Func=@Rastrigin;

	% 2.The whole test loop
    for r=1:Run
        iter=0;
        % 1.Generate MxN matrix
        for m=1:M
            for n=1:N
                X(m,n)=X_min(n)+rand()*(X_max(n)-X_min(n));
            end
            fprintf('value of X:');
            disp(X);
        end
		
    
        for i=1:I_max  % Stop when the iteration large than the max iteration time

            iter=iter+1;
            for m=1:M % For each individual
                % Mutation
                
				[V]=rand1(X,M,F,m);
                % Check if the element in the V matrix beyond the boundary.
                fprintf('V matrix rand: ');
                disp(V);
                for n=1:N
                    if V(1,n)>X_max(1,n)
                        V(1,n)=X_max(1,n);
                    end
                    if V(1,n)<X_min(1,n)
                        V(1,n)=X_min(1,n);
                    end
                    fprintf('V matrix: ');
                    disp(V);
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
                    fprintf('Offspring U: ');
                    disp(U);
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
                
            end % Now the 1th individual generated

            % Select the lowest fitness value
            [y,ind1]=sort(Y,1);
			Y_min=y(1,1);
            
            [Ymin,ind] = min(Y);
            
            
            fprintf('fittest values: ');
            disp(Y);
            
            fprintf('One fittest value: ');
            disp(Ymin);
            
            
            % plot the picture of iteration
            figure(2);
            plot(iter,Ymin,'r.');
            xlabel('Iteration');
            ylabel('Fitness');
			title(sprintf('Iteration=%d, Fitness=%9.9f',i,Ymin));
            grid on;
            hold on;
			
        end % Finish I_max times iteration
	    
		hold off;
        PlotR();
        hold on;
        scatter3(X(ind,1),X(ind,1),Ymin,'fill','ro');

    end % Run 30 times