%% ===================== TRAIN NN TO PREDICT MAX SLANT WITH OVERFIT CHECK =====================
clear; clc; close all;

%% Load dataset
dataFile = 'satellite_training_data.csv';  % your CSV from previous step
tbl = readtable(dataFile);

% Separate features and target
X = tbl{:,1:5};          % sin_tod, cos_tod, sin_doy, cos_doy, frac_visible
y = tbl.max_slant_m;     % target

% Normalize features
Xnorm = X;
Xnorm(:,1:4) = X(:,1:4); % time features already ~[-1,1]
Xnorm(:,5) = X(:,5);     % fraction visible [0,1]

% Split data into 70% train, 15% validation, 15% test
cv = cvpartition(size(X,1),'HoldOut',0.3);
idxTrain = training(cv); idxRest = test(cv);
XTrain = Xnorm(idxTrain,:)';
yTrain = y(idxTrain)';
XRest  = Xnorm(idxRest,:)';
yRest  = y(idxRest)';

% Further split rest into validation and test (50%-50%)
Nrest = length(yRest);
valInd = 1:round(0.5*Nrest);
testInd = round(0.5*Nrest)+1:Nrest;
XVal = XRest(:,valInd); yVal = yRest(valInd);
XTest = XRest(:,testInd); yTest = yRest(testInd);

%% Define feedforward neural network
hiddenSizes = [32 16]; % two hidden layers
net = feedforwardnet(hiddenSizes);

% Configure training parameters
net.trainFcn = 'trainlm'; % Levenberg-Marquardt
net.performFcn = 'mse';
net.divideFcn = 'divideind';
net.divideParam.trainInd = 1:size(XTrain,2);
net.divideParam.valInd   = 1:size(XVal,2);
net.divideParam.testInd  = [];
net.trainParam.max_fail = 6; % early stopping based on validation

%% Train network
[net,tr] = train(net,XTrain,yTrain);

%% Predictions
yPredTrain = net(XTrain);
yPredVal   = net(XVal);
yPredTest  = net(XTest);

%% Compute errors
MAE_train = mean(abs(yPredTrain - yTrain));
MSE_train = mean((yPredTrain - yTrain).^2);
RMSE_train = sqrt(MSE_train);

MAE_val = mean(abs(yPredVal - yVal));
MSE_val = mean((yPredVal - yVal).^2);
RMSE_val = sqrt(MSE_val);

MAE_test = mean(abs(yPredTest - yTest));
MSE_test = mean((yPredTest - yTest).^2);
RMSE_test = sqrt(MSE_test);

fprintf('Training Errors:\nMAE = %.2f m | MSE = %.2f m^2 | RMSE = %.2f m\n', MAE_train,MSE_train,RMSE_train);
fprintf('Validation Errors:\nMAE = %.2f m | MSE = %.2f m^2 | RMSE = %.2f m\n', MAE_val,MSE_val,RMSE_val);
fprintf('Test Errors:\nMAE = %.2f m | MSE = %.2f m^2 | RMSE = %.2f m\n', MAE_test,MSE_test,RMSE_test);

% Overfit check
if RMSE_val > 1.5*RMSE_train
    warning('Possible overfitting detected: validation RMSE is significantly higher than training RMSE.');
else
    disp('No strong signs of overfitting detected.');
end

%% Plot results
figure;
plot(yTest/1e3, 'b.-', 'DisplayName','True Max Slant (km)'); hold on;
plot(yPredTest/1e3, 'r.-', 'DisplayName','NN Prediction (km)');
xlabel('Sample'); ylabel('Max Slant Range (km)');
grid on; legend('Location','best'); title('Neural Network Prediction vs Truth');

figure;
plot((yPredTest - yTest)/1e3, 'k.-'); grid on;
xlabel('Sample'); ylabel('Error (km)'); title('Prediction Error');

disp('Training, evaluation, and overfit check complete.');
