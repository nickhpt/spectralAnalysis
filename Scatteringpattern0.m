function Scatteringpattern0(data,pstruct,nav,coord,Range,window_length,Normalized)
%Plots scattering pattern of input level 1 data
%Input: 1)data, 2)p-struct, 3)navigation file, 4)Sceene coordinates, (vector), 5)Range
%in ice (vector), 6)PSD window length, 7)Normalization fx. 'yes';

dat0=data;
p=pstruct;
c1=coord(1);c2=coord(2);
d0=Range(1);d=Range(2);
c_light=physconst('LightSpeed');
lambda=c_light/p.Fc;
vel=sqrt(nav.ve.^2+nav.vn.^2);
%Mean velocity in scene
vel_mean = mean(vel(c1:c2));
%Frequency resolution 
df=(p.PRF)/length(data);
%Doppler frequency vector
fd=-p.PRF/2:df:p.PRF/2-df;
%Inversion of (fd=2vsin(theta)/lambda)
theta = (asin(fd.*lambda/(2*vel_mean)))*180/pi;

incr=1;
RaSpacing = 1/2*c_light*1/p.Fs;
delay = 1/2*c_light*p.RxDelay;

    for depth=d0:-1:d-1
        Range_cells=round((nav.dif-delay-depth)/RaSpacing);
        Ra_cells(:,incr) = Range_cells;

        for i=1:p.Naz
            cell = Range_cells(i);
            surface_dat(i) = dat0(cell,i);         
        end

    data_volume(:,incr) = fftshift(PSD_Welch2(surface_dat(c1:c2),50,window_length));
 
    if contains(Normalized,"yes")
            max_val = max(data_volume(:,incr));
        else
            max_val = 1;
    end

    data_volume(:,incr) = data_volume(:,incr)/max_val;
    incr=incr+1;
    
    end
    
if contains(Normalized,'yes'); string = 'Normalized Power [dB]';else; string ='Power [dB]'; end

Range=0:1/p.Fs*0.5*c_light:size(data_volume,2)*1/p.Fs*0.5*c_light-1/p.Fs*0.5*c_light;
imagesc(theta,Range,10*log10(data_volume.'));
h = colorbar;
ylabel(h, string)
xlabel('Look angle \theta [deg]');ylabel('Range [m]');
set(gca,'FontSize',12)
end

