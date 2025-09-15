% ==========================================================
% Script: max_slant_range.m
% Purpose: Predict maximum slant range between 3 satellites
% Author: ChatGPT
% ==========================================================

clear; clc;

%% --- INPUT PARAMETERS ---
% Time setup
startTime = datetime(2025,8,15,0,0,0);
endTime   = datetime(2025,8,15,1,0,0);  % 1-hour window
timeStep  = seconds(60);                % 1-minute step
timeVec   = startTime:timeStep:endTime;

% Example: satellite positions in ECEF [m] for each time step
% Replace this with your ephemeris propagation results
% satPos{i}(t,:) = [x, y, z] in meters

% Example dummy data (replace with actual)
nSteps = numel(timeVec);
sat1Pos = repmat([6771e3, 0, 0], nSteps, 1); % ~LEO altitude
sat2Pos = repmat([0, 6771e3, 0], nSteps, 1);
sat3Pos = repmat([0, 0, 6771e3], nSteps, 1);

%% --- COMPUTATION ---
maxSlantRange = 0;
maxPair = "";
maxTime = NaT;

for k = 1:nSteps
    % Positions at this time
    p1 = sat1Pos(k,:);
    p2 = sat2Pos(k,:);
    p3 = sat3Pos(k,:);

    % Compute distances between pairs
    d12 = norm(p1 - p2);
    d13 = norm(p1 - p3);
    d23 = norm(p2 - p3);

    % Find maximum for this time step
    [dMax, idx] = max([d12, d13, d23]);

    % Update global maximum if needed
    if dMax > maxSlantRange
        maxSlantRange = dMax;
        maxTime = timeVec(k);
        switch idx
            case 1, maxPair = "Sat1–Sat2";
            case 2, maxPair = "Sat1–Sat3";
            case 3, maxPair = "Sat2–Sat3";
        end
    end
end

%% --- OUTPUT ---
fprintf('Maximum slant range: %.2f km\n', maxSlantRange/1e3);
fprintf('Between: %s\n', maxPair);
fprintf('At time: %s\n', datestr(maxTime));

