startTime = datetime(2023,5,5,0,0,0);
stopTime = startTime + hours(2);
sampleTime = 30;

sc = satelliteScenario(startTime,stopTime,sampleTime);
sat1 = satellite(sc,"Zarya_TLE.txt");
play(sc);