startTime = datetime(2025,6,28,0,0,0);
stopTime = startTime + hours(5);
sampleTime = 1;

sc = satelliteScenario(startTime,stopTime,sampleTime);

sat = satellite(sc,"threeSatelliteConstellation.tle");
show(sat)
groundTrack(sat,"LeadTime",1200);

ele1 = orbitalElements(sat(1))
ele2 = orbitalElements(sat(2))
ele3 = orbitalElements(sat(3))


