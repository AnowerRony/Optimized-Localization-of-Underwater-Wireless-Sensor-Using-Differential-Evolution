function [Err,Subset] = DEFS_chs(data_tr,data_ts,DNC,PSIZE,Ld,classif,GEN,NFPC)
%%
%    Differential Evolution based Feature Selection
%    Inputs
%    ------
%          data_tr: training dataset (with NP1 patterns x NF+1 features with last column being the training class label)
%          data_ts: testing dataset (with NP1 patterns x NF+1 features with last column being the testing class label)
%          NP:   number of patterns (NP1 and NP2 used for example)
%          NF:   number of features
%          DNC:  desired number of channels to be selected
%          PSIZE:population size
%          Ld:   either load initial population (Ld=1) or simply initialize a new population (Ld=0)
%          classif: takes text value as: 'LDA' or 'KNN' or 'NB' or 'RegTree' 
%          GEN:     number of generations or iterations
%          NFPC: Number of features per each channel
%
%    OutPuts
%    -------
%          Err: Achieved error rate across the different itreations
%          Subset: selected feature subset (feature indices)
%
%    Example:
%    -------
%          load iris.dat
%          [Err,Subset] = DEFS(iris(1:2:end,1:end),iris(2:2:end,:),3,50,0,0,100)
%
%    References:
%    [1] R. N. Khushaba, A. Al-Ani, and A. Al-Jumaily, "Feature subset selection using differential evolution and a statistical repair mechanism",
%        Expert Systems with Applications, vol. 38, no. 9, pp. 11515-11526, 2011.
% 
% DEFS by Dr. Rami Khushaba
% Research Fellow - Faculty of Engineering and IT
% University of Technology, Sydney.
% URL: www.rami-khushaba.com
%%
% CONTROL PARAMETERS %
D = DNC; % dimension of problem
NP = PSIZE; % size of population
CR = 0.5; % crossover constant
L = 1; % low boundary constraint
H = (size(data_tr,2)-1)/NFPC; % high boundary constraint
NF = H;
NE = 5;
% Ld
if nargin < 5 || isempty(Ld),
    Ld = 0;
end
if nargin < 6 || isempty(classif),
    classif = 0;
end
if nargin < 7 || isempty(GEN),
    GEN = 400; % number of generations=iterations
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
w_start = 0.95;                          %Initial inertia weight's value
w_end = 0.35;                            %Final inertia weight
w_varyfor = 1;
w_varyfor = floor(w_varyfor*GEN);       %Weight change step. Defines total number of iterations for which weight is changed.
w_now = w_start;
inertdec = (.95-.35)/w_varyfor;           %Inertia weight's change per iteration
w_start = 0.35;                          %Initial inertia weight's value
w_end = 0.95;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% *************************** %
% ** ALGORITHM VARIABLES ** %
% *************************** %
X = zeros(D,1); % trial vector
Pop = zeros(D,NP); % population
Fit = zeros(1,NP); % fitness of the population
r = zeros(3,1); % randomly selected indices
% *********************** %
% ** CREATE POPULATION ** %
% *********************** %
% initialize random number generator
rand('state',sum(100*clock));
if Ld,
    load TabK Pop
    %Pop = Tab';
else
    for j=1:NP,
        FF = randperm(H);
        Pop(:,j) = FF(1:D)'; % within b.constraints
    end
end
for j = 1:NP % initialize each individual
    val =round(Pop(:,j))';
    v2 = val*NFPC;
    v1 = v2-NFPC+1;
    val =[];
    for chidx = 1:length(v2)
        val =[val v1(chidx):v2(chidx)];
    end
    switch classif
        case 0
            %LDA classifier
            Ac1 = classify(data_ts(:,val),data_tr(:,val),data_tr(:,end),'quadratic');
            Fit(1,j)=sum(Ac1~=data_ts(:,end))/size(data_ts,1);
        case 1
            %kNN classifier
            Ac1 = knnclassify(data_ts(:,val),data_tr(:,val),data_tr(:,end),3);
            Fit(1,j)=sum(Ac1~=data_ts(:,end))/size(data_ts,1);
        otherwise
            %NaiveBase Classifier
            O1 = NaiveBayes.fit(data_tr(:,val),data_tr(:,end));
            Ac1=O1.predict(data_ts(:,val));
            Fit(1,j)=sum(Ac1~=data_ts(:,end))/size(data_ts,1);
    end
end

%load Tab.mat Pop Fit
[FFFF,x22]=sort(Fit);
Best = Pop(:,x22(1:NE));
iBest = x22;
NF1 = max([D - 5, round(D*0.65)]);

RR = zeros(1,NF);
for i=1:PSIZE
    TEMP(i,:) =zeros(1,NF);
    TEMP(i,Pop(:,i)')=1;
end
LL=[];EL=[];
LL = [LL; TEMP];                            % list of priviously tested subsets
EL = [EL;Fit'];                             % Error rates for LL subsets
PD = sum(LL((find(EL<mean(EL))),:));        % Positive feature distribution
ND = sum(LL((find(EL>=mean(EL))),:));       % Negative feature distribution
RR(find(PD+ND)==0) = 0;
RR(find(PD+ND)) = PD(find(PD+ND))./(PD(find(PD+ND))+ND(find(PD+ND)));
First = 1+ (1*RR + 1.55*(PD./max(PD)) + 0.7*(1-(PD+ND)/max(PD+ND))) ;
First = First./max(First);
Second =zeros(1,NF);
% ****************** %
% ** OPTIMIZATION ** %
% ****************** %
disp(sprintf('Iter: %d \t Acc: %.4f \t Subset Selected: %s',1,Fit(iBest(1))*100,num2str(sort(round(Pop(:,iBest(1))')))))
Err(1) = (Fit(iBest(1)))*100;
for g = 2:GEN % for each generation
    expectation = fitscalingrank(Fit,NP);
    if size(expectation,2)>size(expectation,1)
        expectation=expectation';
    end
    parents = selectionstochunif(expectation,NP);
    parents =parents(randperm(length(parents)));
    OLD =Pop;
    for j = 1:NP % for each individual
        % choose three random individuals from population,
        % mutually different and different from j
        if (g<=w_varyfor) && (g > 1)
            w_now = w_now + inertdec; %Change inertia weight
        end
        r(1) = floor(rand()* NP) + 1;
        while r(1)==j
            r(1) = floor(rand()* NP) + 1;
        end
        r(2) = floor(rand()* NP) + 1;
        while (r(2)==r(1))||(r(2)==j)
            r(2) = floor(rand()* NP) + 1;
        end
        r(3) = floor(rand()* NP) + 1;
        while (r(3)==r(2))||(r(3)==r(1))||(r(3)==j)
            r(3) = floor(rand()* NP) + 1;
        end
        %-----------------------------------------------------------%
        %                    Roulette Wheel                         %
        %-----------------------------------------------------------%
        
        % create trial individual
        % in which at least one parameter is changed
        Rnd = floor(rand()*D) + 1;
        rnd = randperm(DNC);
        zq = randperm(NE);
        for i = 1:D
            if ( rand()>CR ) || ( Rnd==i )
                X(i) = Pop(rnd(i),j) + (.5./sum([Pop(rnd(i),r(2))  Best(rnd(i),zq(1))]))  * (Best(rnd(i),zq(1)) - Pop(rnd(i),r(2)) ) ;%+ w_now*(.5*rand) * (Best(i,1) - Pop(i,parents(j)));
            else
                X(i) = Pop(rnd(i),parents(j));
            end
        end
        % verify boundary constraints
        for i = 1:D
            if (X(i)<L)||(X(i)>H)
                X(i) = L + (H-L)*rand();
            end
        end
        % select the best individual
        % between trial and current ones
        % calculate fitness of trial individual
        X=X';
        S1 = unique(round(X));
        nkkk = length(S1);
        if nkkk<D
            Second = 1+ (1*RR + 1.55*(PD./max(PD)) + ((NF-DNC)/NF).*(1-(PD+ND)/max(PD+ND))) ;
            Second = Second ./max(Second);
            Tuu = (Second - First).*Second + First;
            Tuu = Tuu - 0.5.*rand(1,NF).*(1-Tuu);
            expectation = fitscalingrank(Tuu,4*D);
            if size(expectation,2)>size(expectation,1)
                expectation=expectation';
            end
            parentss = selectionstochunif(expectation,4*D);
            parentss =parentss(randperm(length(parentss)));
            ssss = unique(parentss);ssss=setdiff(ssss,S1);
            [FFFF,ts2]=sort(Tuu(ssss),'descend');
            ssss=ssss(ts2);
            X(1,1:DNC)=[S1(randperm(length(S1))) ssss(1:D-length(S1))];
        end
        
        X=X';
        val =round(X)';
        v2 = val*NFPC;
        v1 = v2-NFPC+1;
        val =[];
        for chidx = 1:length(v2)
            val =[val v1(chidx):v2(chidx)];
        end

        switch classif
            case 0
                %LDA classifier
                Ac1 = classify(data_ts(:,val),data_tr(:,val),data_tr(:,end),'quadratic');
                f=sum(Ac1~=data_ts(:,end))/size(data_ts,1);
            case 1
                %kNN classifier
                Ac1 = knnclassify(data_ts(:,val),data_tr(:,val),data_tr(:,end),3);
                f=sum(Ac1~=data_ts(:,end))/size(data_ts,1);
            otherwise
                %NaiveBase Classifier
                O1 = NaiveBayes.fit(data_tr(:,val),data_tr(:,end));
                Ac1=O1.predict(data_ts(:,val));
                f=sum(Ac1~=data_ts(:,end))/size(data_ts,1);
        end
        % if trial is better or equal than current
        if f <= Fit(j)
            Pop(:,j) = X; % replace current by trial
            Fit(j) = f;
            % if trial is better than the best
            if f <= Fit(iBest)
                iBest = j ; % update the best�s index
                w_now = .35;
            end
        end
    end
    for i=1:PSIZE
        TEMP(i,:) =zeros(1,NF);
        TEMP(i,round(Pop(:,i)'))=1;
    end
    ds =find(Fit'<EL);
    LL(ds,:) = [TEMP(ds,:)];
    EL(ds) = [Fit(ds)'];
    
    PD =  (LL((find(EL<mean(EL))),:));
    if size(PD,1)>1
        PD = sum(PD);
    else
        PD =2*rand(1,NF);
    end
    ND =  (LL((find(EL>=mean(EL))),:));
    if size(ND,1)>1
        ND = sum(ND);
    else
        ND = 2*rand(1,NF);
    end
    
    RR(find(PD+ND)==0) = 0;
    RR(find(PD+ND)) = PD(find(PD+ND))./(PD(find(PD+ND))+ND(find(PD+ND)));
    AA = (1-(PD+ND)/max(PD+ND));
    [temp,iBest] =sort(Fit);
    Best = Pop(:,iBest(1:NE));
    disp(sprintf('Iter: %d \t Err: %.4f \t Subset Selected: %s',g,Fit(iBest(1))*100,num2str(sort(round(Pop(:,iBest(1))')))))
    First = Second;
    Err(g) = (Fit(iBest(1)))*100;
end
% ************* %
% ** RESULTS ** %
% ************* %
f = Fit(iBest(1));
Subset = sort(round(Pop(:,iBest(1))))';



function parents = ssuni(expectation,nParents)


wheel = cumsum(expectation) / nParents;

parents = zeros(1,nParents);

% we will step through the wheel in even steps.
stepSize = 1/nParents;

% we will start at a random position less that one full step
position = rand * stepSize;

% a speed optimization. Position is monotonically rising.
lowest = 1; 

for i = 1:nParents % for each parent needed,
    for j = lowest:length(wheel) % find the wheel position
        if(position < wheel(j)) % that this step falls in.
            parents(i) = j;
            lowest = j;
            break;
        end
    end
    position = position + stepSize; % take the next step.
end
   
