function [InitialObservation,LoggedSignals] = myResetFcn()

    % Reset environment state
    LoggedSignals.State = 0;

    % Initial observation
    InitialObservation = LoggedSignals.State;

end