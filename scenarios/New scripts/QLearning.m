%% Parameters
numEpisodes = 5000;
alpha = 0.1;     % learning rate
gammaRL = 0.9;   % discount factor
epsilon = 0.1;   % exploration probability

% Discretization (for simplicity)
U_bins = linspace(0,2,10);
h_bins = linspace(0,5,10);
EbNo_bins = linspace(0,10,10);

% Initialize Q-table: states x actions
Q = zeros(length(U_bins), length(h_bins), length(EbNo_bins), 2);

%% Training loop
for ep = 1:numEpisodes
    [state, rewardOptions] = channelStep();

    % Discretize state
    [~,iU] = min(abs(U_bins - state(1)));
    [~,iH] = min(abs(h_bins - state(2)));
    [~,iE] = min(abs(EbNo_bins - state(3)));

    % ε-greedy action
    if rand < epsilon
        action = randi([0 1]);  % explore
    else
        [~,action] = max(Q(iU,iH,iE,:));  % exploit
        action = action-1;  % adjust to 0/1
    end

    % Get reward for chosen action
    r = rewardOptions(action+1);

    % Next state
    [nextState, rewardOptionsNext] = channelStep();
    [~,iU2] = min(abs(U_bins - nextState(1)));
    [~,iH2] = min(abs(h_bins - nextState(2)));
    [~,iE2] = min(abs(EbNo_bins - nextState(3)));

    % Update Q
    Q(iU,iH,iE,action+1) = Q(iU,iH,iE,action+1) + ...
        alpha * (r + gammaRL * max(Q(iU2,iH2,iE2,:)) - Q(iU,iH,iE,action+1));
end
