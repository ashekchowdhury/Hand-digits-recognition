% %                           Practical Assignment of Pattern Recognition 
clc
clear all
close all

%% Load data
inputdata=zeros(100,3,1000); % 1000 according to input data; preprocessing and 300 since different amount of datapoints available per stroke
sizeOfData=[100 100 100 100 100 100 100 100 100 100];
dataNum=[repmat(0,sizeOfData(1),1); repmat(1,sizeOfData(2),1); repmat(2,sizeOfData(3),1); ...
    repmat(3,sizeOfData(4),1); repmat(4,sizeOfData(5),1); repmat(5,sizeOfData(6),1);...
    repmat(6,sizeOfData(7),1); repmat(7,sizeOfData(8),1); repmat(8,sizeOfData(9),1);...
    repmat(9,sizeOfData(10),1)]; % classes-vector
sizeOfDataCum=cumsum(sizeOfData);
sizeOfDataCum=[0 sizeOfDataCum]; % add zero for the loop
strokeSize=[]; % rows of the measurement size they differ in the dataset

for i = 1:10 % for numbers 0 to 9
        for j=1:sizeOfData(i)
            if j<10
            C = strcat('stroke_',num2str(i-1),'_000',num2str(j),'.mat');
            elseif (j>=10 && j<100)
                C = strcat('stroke_',num2str(i-1),'_00',num2str(j),'.mat');
                else C = strcat('stroke_',num2str(i-1),'_0',num2str(j),'.mat');
            end
            load(C);
            strokeSize=[strokeSize, size(pos,1)]; % up to 1000 strokes
            for k=1:size(pos,2)
            pos(:,k)=(pos(:,k)-min(pos(:,k)))/(max(pos(:,k))-min(pos(:,k)));% max-min normalization
            end
            inputdata(1:size(pos,1),:,sizeOfDataCum(i)+j)=pos;
        end
end

newData=inputdata;
save ('newData.mat','newData')
save('strokeSize.mat', 'strokeSize')
save('dataNum.mat','dataNum')
%% Rotation on the x-axis

degOfRot=[5,10,15,20,25,30,35,40,45]; % degress of rotation
nRotXDataSet=zeros(100,3,size(inputdata,3)*size(degOfRot,2)); % preprocessing
strokeSizeXRot=repmat(strokeSize,1,size(degOfRot,2));
for i=1:size(degOfRot,2) % size of rotation vector
for j=1:size(inputdata,3) % size of data (3rd dim)
% https://se.mathworks.com/matlabcentral/newsreader/view_thread/308859
theta = deg2rad(degOfRot(i));
R3DRot = [cos(theta) 0 sin(theta); 0 1 0; -sin(theta) 0 cos(theta)]; % rotation around x
dataSetOfRot = inputdata(1:strokeSize(j),:,j)*R3DRot;
nRotXDataSet(1:size(dataSetOfRot,1),:,(i-1)*size(inputdata,3)+j)=bsxfun(@minus,dataSetOfRot,min(dataSetOfRot,[],1))./repmat((max(dataSetOfRot,[],1)-min(dataSetOfRot,[],1)),size(dataSetOfRot,1),1); % max min
% in each step of j one new matrix is adapted from the zeros matrix in the
% third dimension of the NRotXdataset matrix
end
end

Numdataxrot=repmat(dataNum,size(degOfRot,2),1);

save ('nRotXDataSet.mat','nRotXDataSet')
save ('strokeSizeXRot.mat','strokeSizeXRot')
save ('Numdataxrot.mat','Numdataxrot')


%% Rotation on the y-axis

degOfRot=[5,10,15,20,25,30,35,40,45]; % degress of rotation
nRotYDataSet=zeros(100,3,size(inputdata,3)*size(degOfRot,2)); % preprocessing
strokeSizeYRot=repmat(strokeSize,1,size(degOfRot,2));
for i=1:size(degOfRot,2) % size of rotation vector
for j=1:size(inputdata,3) % size of data (3rd dim)
% https://se.mathworks.com/matlabcentral/newsreader/view_thread/308859
theta = deg2rad(degOfRot(i));
R3DRot =[cos(theta),0,-sin(theta);0,1,0;sin(theta),0,cos(theta)]; % rotation around y (different formula)
dataSetOfRot = inputdata(1:strokeSize(j),:,j)*R3DRot;
nRotYDataSet(1:size(dataSetOfRot,1),:,(i-1)*size(inputdata,3)+j)=bsxfun(@minus,dataSetOfRot,min(dataSetOfRot,[],1))./repmat((max(dataSetOfRot,[],1)-min(dataSetOfRot,[],1)),size(dataSetOfRot,1),1); % max min
% in each step of j one new matrix is adapted from the zeros matrix in the
% third dimension of the NRotYdataset matrix
end
end

Numdatayrot=repmat(dataNum,size(degOfRot,2),1);

save ('nRotYDataSet.mat','nRotYDataSet') 
save ('strokeSizeYRot.mat','strokeSizeYRot') 
save ('Numdatayrot.mat','Numdatayrot')



%% Rotation on the z-axis

degOfRot=[5,10,15,20,25,30,35,40,45]; % degress of rotation
nRotZDataSet=zeros(100,3,size(inputdata,3)*size(degOfRot,2)); % preprocessing
strokeSizeZRot=repmat(strokeSize,1,size(degOfRot,2));
for i=1:size(degOfRot,2) % size of rotation vector
for j=1:size(inputdata,3) % size of data (3rd dim)
% https://se.mathworks.com/matlabcentral/newsreader/view_thread/308859
theta = deg2rad(degOfRot(i));
R3DRot =[cos(theta),sin(theta),0;-sin(theta),cos(theta),0;0,0,1]; % rotation around z (different formula)
dataSetOfRot = inputdata(1:strokeSize(j),:,j)*R3DRot;
nRotZDataSet(1:size(dataSetOfRot,1),:,(i-1)*size(inputdata,3)+j)=bsxfun(@minus,dataSetOfRot,min(dataSetOfRot,[],1))./repmat((max(dataSetOfRot,[],1)-min(dataSetOfRot,[],1)),size(dataSetOfRot,1),1); % max min
% in each step of j one new matrix is adapted from the zeros matrix in the
% third dimension of the NRotZdataset matrix
end
end

Numdatazrot=repmat(dataNum,size(degOfRot,2),1);

save ('nRotZDataSet.mat','nRotZDataSet') 
save ('strokeSizeZRot.mat','strokeSizeZRot') 
save ('Numdatazrot.mat','Numdatazrot')


%% Combine the datasets to 1
sizedata=size(newData,3)+size(nRotXDataSet,3)+size(nRotYDataSet,3)+size(nRotZDataSet,3);
NdataCom=zeros(106,3,sizedata);
NdataCom(:,:,1:size(newData,3))=newData;
NdataCom(:,:,size(newData,3)+1:size(newData,3)+size(nRotXDataSet,3))=nRotXDataSet;
NdataCom(:,:,size(newData,3)+size(nRotXDataSet,3)+1:size(newData,3)+size(nRotXDataSet,3)+size(nRotYDataSet,3))=nRotYDataSet;
NdataCom(:,:,size(newData,3)+size(nRotXDataSet,3)+size(nRotYDataSet,3)+1:sizedata)=nRotZDataSet;
strokeSizeCom=[strokeSize strokeSizeXRot strokeSizeYRot strokeSizeZRot];

NumdataCom=[dataNum; Numdataxrot; Numdatayrot; Numdatazrot]; % all classes in one vector

save ('NdataCom.mat','NdataCom') 
save ('strokeSizeCom.mat','strokeSizeCom') 
save ('NumdataCom.mat','NumdataCom')
%% Adpation 10% Approach

slopeval=zeros(4,size(strokeSizeCom,2));
w=0.20;     % 10% of initial data

for m=1:size(strokeSizeCom,2)
    slopeval(1:2,m)=[NdataCom(round(strokeSizeCom(m)*w),2,m)-NdataCom(1,2,m); NdataCom(round(strokeSizeCom(m)*w),1,m)-NdataCom(1,1,m)];
    slopeval(3,m)=sum(diff(NdataCom(1:strokeSizeCom(m),1,m)));
    slopeval(4,m)=sum(diff(NdataCom(1:strokeSizeCom(m),2,m)));
end

%% Use of classifier 
load NdataCom.mat % data
load strokeSizeCom.mat
load NumdataCom.mat % the class

% FOR PROJECTION ON Y and X
% define equations --> find parameter
stepsize=0.0500; % since 3 decimal numbers for the recording of the coordinates per axis
yproject=zeros(size(0:stepsize:1,2),1); % preprocessing for the amount of projections
xproject=zeros(size(0:stepsize:1,2),1); % preprocessing for the amount of projections
profilsnumy=zeros(size(0:stepsize:1,2),size(NdataCom,3)); % preprocessing for the matrix for the amount of projections for all observations
profilsnumx=zeros(size(0:stepsize:1,2),size(NdataCom,3)); 
for k=1:size(NdataCom,3) 
    for y=0:stepsize:1
        n=0; % will be used for summation of the amount of projections (for projection on y)
        m=0; % will be used for summation of the amount of projections (for projection on x)
        for p=1:strokeSizeCom(k)-1% '-1' since also use of p+1
             if NdataCom(p,1,k)<NdataCom(p+1,1,k) % this if-function makes sure that so that the first point (x1) is smaller the second point (x2) the solver can find a slope
                x1=NdataCom(p,1,k);
                x2=NdataCom(p+1,1,k);
                y1=NdataCom(p,2,k);
                y2=NdataCom(p+1,2,k);
            else
                x2=NdataCom(p,1,k);
                x1=NdataCom(p+1,1,k);
                y2=NdataCom(p,2,k);
                y1=NdataCom(p+1,2,k);
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
    profilsnumy(:,k)=yproject;
    profilsnumx(:,k)=xproject;
end

trainprofile=[profilsnumx; profilsnumy]; % train vector for the profiles
trainprofile=[trainprofile; slopeval];
save('trainprofile.mat','trainprofile')
load dataNum.mat % the class
% % Adapt the trainclass 0-9 --> 1-10

trainratio=0.7; % size of training dataset
testratio=1-trainratio; % size of test dataset

idxTrain=randsample(size(newData,3),round(trainratio*size(newData,3))); % index of train data randomly determined
idxTest=setdiff(1:size(newData,3),idxTrain)'; % index of test data randomly determined

trainprof=trainprofile(:,idxTrain); % traindata
testprof=trainprofile(:,idxTest); % testdata
trainclassprof=dataNum(idxTrain,1); % classes traindata
testclassprof=dataNum(idxTest,1); % classes testdata

class=knn(trainclassprof,trainprof,testprof,7);
correctclass=sum(class==testclassprof)/size(testclassprof,1); %0.9483 %0.9310
avgcorrectclass=mean(correctclass);
accuracy =((correctclass)/(class))*100
%trainclass,traindata,data,k
%% Trial of the final classifier
load('stroke_8_0099.mat')
X1= pos (:,1);
Y1=pos (:,2);
plot (X1,Y1)
Class=digit_classify(pos)
