function Scatteringpattern1(data,pstruct,nav,coord,Range,window_length,Normalized)
%Plots scattering pattern of input level 1 data
%Input: 1)data, 2)p-struct, 3)navigation file, 4)Sceene coordinates, (vector), 5)Range
%in ice (vector), 6)PSD window length, 7)Normalization fx. 'yes';

dat1=data;
p=pstruct;  
nav=nav;
c1=coord(1);c2=coord(2);
d0=Range(1);d=Range(2);
c_light=physconst('LightSpeed');
lambda=c_light/p.Fc;
max_thet = lambda/(4*p.AzSpacing);
incr=1;

    for depth=d0:-1:d-1
        Range_cells=round((p.ZOff-depth-(nav.alt-nav.dif))/abs(p.ZSpacing)).';
        Ra_cells(:,incr) = Range_cells;

        for i=1:p.Naz
            cell = Range_cells(i);
            surface_dat(i) = dat1(cell,i);         
        end

    surface_dat = surface_dat(c1:c2);
    data_volume(:,incr) = fftshift(PSD_Welch2(surface_dat,50,window_length));

    if contains(Normalized,"yes")
            max_val = max(data_volume(:,incr));
        else
            max_val = 1;
    end

    data_volume(:,incr) = data_volume(:,incr)/max_val;
    incr=incr+1;
    
    end
    
if contains(Normalized,'yes'); string = 'Normalized Power [dB]';else; string ='Power [dB]'; end

Range=0:abs(p.ZSpacing):size(data_volume,2)*abs(p.ZSpacing)-1;
N=size(data_volume,1);
theta = linspace(-max_thet,max_thet,N)*180/pi;
imagesc(theta,Range,10*log10(data_volume.'));
set(gca,'FontSize',12)
xlabel('Look angle \theta [deg]');ylabel('Range [m]')
h = colorbar;
ylabel(h, string)
  
end

