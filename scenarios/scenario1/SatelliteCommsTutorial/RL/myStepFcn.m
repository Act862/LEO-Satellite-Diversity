function [NextObs,Reward,IsDone,LoggedSignals] = myStepFcn(Action,LoggedSignals)

    % Example update
    LoggedSignals.State = LoggedSignals.State + Action;

    % Compute next observation
    NextObs = LoggedSignals.State;

    % Example reward
    Reward = -abs(NextObs); 

    % Example termination condition
    IsDone = abs(NextObs) > 10;

end