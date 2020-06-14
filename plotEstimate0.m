function plotEstimate0(data,color,PRF,vel,coord)
%Plots surface scattering pattern of input level 0 PSD data
%Input: 1)Data, 2)color fx. 'b', 3)PRF, 4)Velocity vector, 5)min/max
%coordinates of scene of interest in a vector vel = [x1,x2]

%Color choice
colorCHAR = color;
%Mean velocity in scene
vel_mean = mean(vel(coord(1):coord(2)));
%Wavelength
c_light=physconst('LightSpeed');
lambda=c_light/435000000;
%Frequency resolution 
df=PRF/length(data);
%Doppler frequency vector
fd=-PRF/2:df:PRF/2-df;
%Inversion of (fd=2vsin(theta)/lambda)
theta = (asin(fd.*lambda/(2*vel_mean)))*180/pi;

plot(theta,10*log10(fftshift(data/max(data))),colorCHAR,'LineWidth',1.5);
xlabel('Look angle \theta');ylabel('Power dB')
grid on
end

