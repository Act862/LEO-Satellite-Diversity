function [Observation, Reward, IsDone, LoggedSignals] = localStep(Action, LoggedSignals, obsValid, yValid, Nvalid, maxStepsPerEp)
% Step function with exactly two inputs (Action, LoggedSignals).
% Outputs Observation as 5x1 column, scalar Reward, logical IsDone.

    % Extract current pointer
    ptr = LoggedSignals.ptr;

    % Ensure numeric scalar action in [0,1]
    a = double(Action);
    if iscell(a), a = a{1}; end
    a = min(max(a,0),1);

    % Ground-truth normalized target at current step
    y = yValid(ptr);

    % Reward: negative squared error
    Reward = - (y - a).^2;

    % Advance pointer and check termination
    ptrNext = ptr + 1;
    horizonEnd = min(Nvalid, LoggedSignals.ptr + maxStepsPerEp - 1);
    IsDone = ptrNext > horizonEnd;

    if IsDone
        ptrNext = min(ptrNext, Nvalid);
    end

    % Next observation as 5x1 column
    Observation        = obsValid(ptrNext, :)';
    LoggedSignals.ptr  = ptrNext;
end