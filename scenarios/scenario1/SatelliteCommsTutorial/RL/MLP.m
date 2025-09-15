%% ===================== TRAIN NN TO PREDICT MAX SLANT (SAVE WEIGHTS + OVERFIT CHECK) =====================
clear; clc; close all;

%% Load dataset
dataFile = 'satellite_training_data.csv';  % your CSV from previous step
tbl = readtable(dataFile);

% Separate features and target
X = tbl{:,1:5};          % sin_tod, cos_tod, sin_doy, cos_doy, frac_visible
y = tbl.max_slant_m;     % target

% Normalize features (time features ~[-1,1], frac_visible [0,1])
Xnorm = X;
Xnorm(:,1:4) = X(:,1:4);
Xnorm(:,5) = X(:,5);

% Split data into 80% train, 20% test
cv = cvpartition(size(X,1),'HoldOut',0.2);
idxTrain = training(cv); idxTest = test(cv);

XTrain = Xnorm(idxTrain,:)';
yTrain = y(idxTrain)';
XTest  = Xnorm(idxTest,:)';
yTest  = y(idxTest)';

%% Define feedforward neural network
hiddenSizes = [32 16]; % two hidden layers
net = feedforwardnet(hiddenSizes);

% Configure training parameters
net.trainFcn = 'trainlm'; % Levenberg-Marquardt
net.performFcn = 'mse';
net.divideFcn = 'divideind';
net.divideParam.trainInd = 1:size(XTrain,2);
net.divideParam.valInd   = [];
net.divideParam.testInd  = [];

% Train network
[net,tr] = train(net,XTrain,yTrain);

%% Save trained network
save('trainedMaxSlantNN.mat','net');

%% Predictions
yPredTrain = net(XTrain);
yPredTest  = net(XTest);

%% Compute errors
MAE_train = mean(abs(yPredTrain - yTrain));
MSE_train = mean((yPredTrain - yTrain).^2);
RMSE_train = sqrt(MSE_train);

MAE_test = mean(abs(yPredTest - yTest));
MSE_test = mean((yPredTest - yTest).^2);
RMSE_test = sqrt(MSE_test);

fprintf('Training Errors:\nMAE = %.2f m\nMSE = %.2f m^2\nRMSE = %.2f m\n\n', MAE_train,MSE_train,RMSE_train);
fprintf('Test Errors:\nMAE = %.2f m\nMSE = %.2f m^2\nRMSE = %.2f m\n', MAE_test,MSE_test,RMSE_test);

%% ===================== OVERFIT CHECK =====================
if RMSE_test > 1.5 * RMSE_train
    warning('Potential overfitting detected: Test RMSE is significantly higher than training RMSE.');
end

%% Plot results
figure;
plot(yTest/1e3, 'b.-', 'DisplayName','True Max Slant (km)'); hold on;
plot(yPredTest/1e3, 'r.-', 'DisplayName','NN Prediction (km)');
xlabel('Sample'); ylabel('Max Slant Range (km)');
grid on; legend('Location','best'); title('Neural Network Prediction vs Truth');

figure;
plot((yPredTest - yTest)/1e3, 'k.-', 'DisplayName','Test Error'); hold on;
plot((yPredTrain - yTrain)/1e3, 'r.-', 'DisplayName','Train Error');
xlabel('Sample'); ylabel('Error (km)'); title('Prediction Error Comparison');
legend('Location','best'); grid on;

disp('Training, evaluation, saving, and overfit check complete.');
