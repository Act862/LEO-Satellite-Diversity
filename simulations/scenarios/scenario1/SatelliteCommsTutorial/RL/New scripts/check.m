%   Satellite "first to go out of sight" simulation
clear; clc;
tleFile = "leoSatelliteConstellation.tle";

startTime = datetime("now");
endTime   = startTime + hours(12);   % longer window to capture passes
sampleTime = 60;                     % seconds

%   Initialize the satellite scenario
sc = satelliteScenario(startTime,endTime,sampleTime);

%   Add satellites
sat = satellite(sc,tleFile);

%   Add ground station (Athens)
gsLat = 37.9838; gsLon = 23.7275; gsAlt = 0.1;  % km
gs = groundStation(sc,gsLat,gsLon,"Altitude",gsAlt,...
    "MinElevationAngle",5,"Name","Athens");

% -------------------------------------------------------
% Find the satellite that will FIRST go out of sight
% -------------------------------------------------------
numSat = numel(sat);
endTimes = NaT(1,numSat);

for k = 1:numSat
    ac = access(sat(k), gs);
    intervals = accessIntervals(ac);

    if ~isempty(intervals)
        % The first access interval (row 1, col 2 = end time)
        endTimes(k) = intervals(1,2);
    end
end

% Ignore satellites that never had access
validIdx = ~isnat(endTimes);
if any(validIdx)
    [firstEndTime, satIdx] = min(endTimes(validIdx));
    satNum = find(validIdx);
    satNum = satNum(satIdx);

    fprintf("Satellite %d will be the FIRST to go out of sight at %s\n", ...
        satNum, string(firstEndTime));
else
    disp("No satellites are visible in this window.");
end
