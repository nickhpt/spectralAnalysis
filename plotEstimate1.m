function plotEstimate1(data,color)
%Plots surface scattering pattern of input level 1 PSD data
%Input
%Data, and color fx. 'b'

colorCHAR = color;
c_light=physconst('LightSpeed');
lambda=c_light/435000000;
max_thet = lambda/(4*1.5);

data=data/max(data);
N=length(data);
theta = linspace(-max_thet,max_thet,N);

plot(theta*180/pi,10*log10(fftshift(data)),colorCHAR,'LineWidth',1.5);
xlabel('Look angle \theta');ylabel('Power dB')
grid on
end

