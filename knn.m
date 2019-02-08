function C = knn(trainclass, traindata, data, k)
% function C = knn(trainclass, traindata, data, k)

% Calculate Euclidan distance between (test) data and training points
for i = 1:length(data)
  for j = 1:length(traindata)
    c(j, i) = norm(traindata(:, j) - data(:, i));
  end
end

% Replicate training class for each data row
% so that when data rows are sorted the classinfo don't get lost
class = repmat(trainclass, [size(traindata), 1])';

% Sort the distances to decreasing order - each column contains
% distance from i:th datapoint to all training points
[n_nearest, index] = sort(c);

% Next sort the class information to correspond with
% the new sorted data points
train = class(index);

k_nearest = n_nearest(1:k, :); % choose only k-nearest data points
train = train(1:k, :); % choose the k-nearest classes 

% Perform the classification
for i = 1:length(data)
  for j = 1:max(trainclass) % determine the number of different classes
    % Calculate the number of samples that belong to one class
    temp(j, i) = sum(train(:, i) == j);
  end
end
% Determine the classinfo: the row that has most "hits" will be chosen
% for each data point and each row in temp represents one class
[temp, C] = max(temp);
