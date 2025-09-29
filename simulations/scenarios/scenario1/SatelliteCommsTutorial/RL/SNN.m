%% ===================== TRAIN A SIMPLE SPIKED-LIKE NEURAL NETWORK =====================
clear; clc; close all;

%% ===================== LOAD DATA =====================
filename = 'satellite_training_data.csv';  % CSV saved from previous script
data = readmatrix(filename);   % numeric data only
% CSV header: sin_tod, cos_tod, sin_doy, cos_doy, frac_visible, max_slant_m

X = data(:,1:5)'; % features: 5xN
Y = data(:,6)';   % target: 1xN

%% ===================== NORMALIZE INPUTS =====================
Xmean = mean(X,2);
Xstd  = std(X,0,2) + 1e-6; % avoid div by 0
Xn = (X - Xmean)./Xstd;

Ymean = mean(Y);
Ystd  = std(Y) + 1e-6;
Yn = (Y - Ymean)/Ystd;

%% ===================== DEFINE NETWORK =====================
% Simple 2-layer network to emulate spiking rate coding
hiddenSize = 20;
net = feedforwardnet(hiddenSize,'trainlm'); % Levenberg-Marquardt

% Use ReLU-like activation to emulate spike rate (positive-only)
for i=1:length(net.layers)-1
    net.layers{i}.transferFcn = 'poslin'; 
end
net.layers{end}.transferFcn = 'purelin'; % linear output

net.trainParam.epochs = 200;
net.trainParam.goal   = 1e-4;

%% ===================== TRAIN NETWORK =====================
[net,tr] = train(net,Xn,Yn);

%% ===================== PREDICT =====================
Ypred_n = net(Xn);
Ypred = Ypred_n*Ystd + Ymean;

%% ===================== COMPUTE ERROR =====================
err_m = Ypred - Y;

%% ===================== PLOT RESULTS =====================
figure;
plot(Y,'b.-','DisplayName','True Max Slant'); hold on;
plot(Ypred,'r.-','DisplayName','Predicted Max Slant');
xlabel('Sample'); ylabel('Max Slant (m)');
title('Spiking-like Neural Network Prediction');
legend('Location','best'); grid on;

figure;
plot(err_m,'k.-');
xlabel('Sample'); ylabel('Error (m)');
title('Prediction Error (Spiking-like NN)');
grid on;

disp('Training and prediction complete.');
