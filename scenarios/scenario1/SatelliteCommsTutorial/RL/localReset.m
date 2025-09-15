function [InitialObservation, LoggedSignals] = localReset(obsValid, Nvalid, maxStepsPerEp)
% Reset function with zero input args as required by rlFunctionEnv.
% Returns a 5x1 column InitialObservation to match obsInfo.

    % Pick a random starting pointer, ensuring room for an episode window
    startMax = max(1, Nvalid - maxStepsPerEp + 1);
    currPtr  = randi([1, startMax], 1, 1);

    % Column vector (5x1) as required by obsInfo
    InitialObservation   = obsValid(currPtr, :)';
    LoggedSignals.ptr    = currPtr;
end