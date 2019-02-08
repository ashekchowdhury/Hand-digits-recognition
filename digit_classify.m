function C=digit_classify(testdata) 
% Function that performs kNN classification with given feature extraction
% method
% data entered where columns are x,y and z coordinates while the rows are
% the observational points. Only information in matrix form of one digit
% function in a loop (or adapt so that the third dimension of a matrix is
% used)

% Load the trainingdata and class

load trainprofile; % feature information of existing data (training data)
load NumdataCom;

traindata=trainprofile; % training data
trainclass=NumdataCom;
k=3; % based on crossvalidation

%% Normalization of the data
for p=1:size(testdata,2)
    testdata(:,p)=(testdata(:,p)-min(testdata(:,p)))/(max(testdata(:,p))-min(testdata(:,p)));% max-min normalization
end
%% Projections AND slope information

pos=size(testdata,1); % length of the stroke

% FOR PROJECTION ON Y and X
% define equations --> find parameter
stepsize=0.0500; % since 3 decimal numbers for the recording of the coordinates per axis
yproject=zeros(size(0:stepsize:1,2),1); % preprocessing for the amount of projections
xproject=zeros(size(0:stepsize:1,2),1); % preprocessing for the amount of projections
for y=0:stepsize:1
        n=0; % will be used for summation of the amount of projections (for projection on y)
        m=0; % will be used for summation of the amount of projections (for projection on x)
        for p=1:pos-1% '-1' since also use of p+1
             if testdata(p,1)<testdata(p+1,1) % this if-function makes sure that so that the first point (x1) is smaller the second point (x2) the solver can find a slope
                x1=testdata(p,1);
                x2=testdata(p+1,1);
                y1=testdata(p,2);
                y2=testdata(p+1,2);
            else
                x2=testdata(p,1);
                x1=testdata(p+1,1);
                y2=testdata(p,2);
                y1=testdata(p+1,2);
            end
            a=(y2-y1)/(x2-x1); %slope (y=ax +b) - for y
            b=y1-(a*x1); % intersect - for y
            xval=(y-b)/a;
            if xval>=min([x1, x2]) && xval<max([x1, x2]) % closed and open interval
                n=n+1;
            end
            yval=a*y+b;
            if yval>=min([y1, y2]) && yval<max([y1, y2]) % closed and open interval
                m=m+1;
            end
        end
        yproject(round((y+stepsize)/stepsize))=n;
        xproject(round((y+stepsize)/stepsize))=m;
end

trainprofile=[xproject; yproject]; % train vector for the profiles
%%
slopeval=zeros(4,1);
w=0.20; % 10% of initial data
slopeval(1:2,1)=[testdata(round(pos*w),2)-testdata(1,2); testdata(round(pos*w),1)-testdata(1,1)]; % first y then x
slopeval(3,1)=sum(diff(testdata(1:pos,1)));
slopeval(4,1)=sum(diff(testdata(1:pos,2)));

data=[trainprofile; slopeval]; % feature vector for the test data

%% Determination of the kNN
C=zeros(1,1); % vector of the class based on knn
% Determine the distance
for h=1:size(data,2)
    disttrain=zeros(size(traindata,2),1);
    for j=1:size(traindata,2)
    disttrain(j)=(sum((traindata(:,j)-data).^2))^0.5;
    end
    disttrainsort=sort(disttrain,'ascend');
    knearest=disttrainsort(1:1:k); % distance of k-nearest neighbors
    % Find the corresponding dataclasses
    pos=zeros(size(knearest,1),1);
    n=1;
    while n<k+1
        if size(find(disttrain==knearest(n)),1)==1 % if only one observation with the k-nearest distance
%             m=pos(n);
            pos(n)=find(disttrain==knearest(n));
            n=n+1;
            
        elseif size(find(disttrain==knearest(n)),1)<=k-n+1 % +1 since value of the n-th pos itself is also still to be assigned
%             m=pos(n:1:n+size(find(disttrain==knearest(n)),1)-1);
            pos(n:1:n+size(find(disttrain==knearest(n)),1)-1)=find(disttrain==knearest(n)); % -1 in the position since n itself is also a position
            n=n+1+size(find(disttrain==knearest(n)),1)-1;
        else % case if e.g. for k=5 5 and 6 are the same and not both values can be taken, only the first one will be included
            potv=find(disttrain==knearest(n));
%             m=pos(n:1:k);
            pos(n:1:k)=potv(1:1:k-n+1); 
            n=k+1;
        end
    end
    C=mode(trainclass(pos));% choses the mode of k nearest neighbours. If no single mode, then first class choosen
     
end
end

