%% ===================== LEO SNN Training =====================
clear; clc; close all;

% Load dataset
data = readtable('satellite_training_data.csv');

% Extract inputs and target
X = table2array(data(:,1:5))'; % 5 features x N samples
Y = table2array(data(:,6))';   % 1 output x N samples

% Normalize inputs
X = X ./ max(X,[],2);

% SNN parameters
Nin = size(X,1);      % number of input neurons
Nh  = 20;             % hidden layer neurons
Nout = 1;             % output neuron
Tsim = 20;            % simulation time steps per sample
dt = 1;               % time step (arbitrary units)
Vth = 1;              % firing threshold
tau = 5;              % membrane time constant

% Initialize weights
W1 = randn(Nh, Nin)*0.1;  % input -> hidden
W2 = randn(Nout, Nh)*0.1; % hidden -> output

% Training parameters
eta = 0.01;  % learning rate
epochs = 50;

% Training loop (rate-coded supervised SNN)
for ep = 1:epochs
    epochError = 0;
    
    for n = 1:size(X,2)
        % Input spike trains (rate-coded)
        spikeInput = poissrnd(X(:,n)*dt, Nin, Tsim); % binary spike trains
        
        % Hidden layer
        Vh = zeros(Nh,1); Sh = zeros(Nh, Tsim);
        for t = 1:Tsim
            Vh = Vh + W1*spikeInput(:,t) - Vh/tau;
            Sh(:,t) = Vh >= Vth;
            Vh(Sh(:,t)) = 0; % reset after spike
        end
        rateH = sum(Sh,2)/Tsim; % spike rate per hidden neuron
        
        % Output layer
        Vo = W2*rateH;  % linear combination (rate decoding)
        
        % Error
        e = Y(n) - Vo;
        epochError = epochError + e^2;
        
        % Weight update (rate-coded delta rule)
        W2 = W2 + eta * e * rateH';
        W1 = W1 + eta * (e*W2') .* sum(Sh,2)/Tsim * spikeInput(:,n)'; 
    end
    
    fprintf('Epoch %d, MSE=%.6f\n', ep, epochError/size(X,2));
end

%% ===================== INFERENCE =====================
pred = zeros(1,size(X,2));
for n = 1:size(X,2)
    spikeInput = poissrnd(X(:,n)*dt, Nin, Tsim);
    Vh = zeros(Nh,1); Sh = zeros(Nh, Tsim);
    for t = 1:Tsim
        Vh = Vh + W1*spikeInput(:,t) - Vh/tau;
        Sh(:,t) = Vh >= Vth;
        Vh(Sh(:,t)) = 0;
    end
    rateH = sum(Sh,2)/Tsim;
    pred(n) = W2*rateH;
end

% Plot results
figure;
plot(Y, 'b.-'); hold on;
plot(pred, 'r.-'); grid on;
xlabel('Sample'); ylabel('Max Slant Range (normalized)');
legend('True','SNN Prediction'); title('SNN Rate-coded Prediction');

% Calculate error
err = pred - Y;
fprintf('Overall MSE: %.6f\n', mean(err.^2));
