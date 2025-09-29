function angle = getUang(velocity,time,radius,altitude)
    angle = velocity*time./(radius+altitude);
end