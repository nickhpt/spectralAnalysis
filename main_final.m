close all; clear all; clc
%Read Polaris data
[dat1_KM,p1_KM]=read_polaris('p110219_m125706_KM5KM6_1c_fvv0');
[dat1_EP,p1_EP]=read_polaris('p110218_m100213_epica_1c_svv0');
[dat0_KM,p0_KM]=read_polaris('p110219_m125706_KM5KM6_0c_fvv0');
[dat0_EP,p0_EP]=read_polaris('p110218_m100213_epica_0c_svv0');
nav1_KM=read_nav('p110219_m125706_KM5KM6_1c_fvv0_nav');
nav1_EP=read_nav('p110218_m100213_epica_1c_svv0_nav');
nav0_KM=read_nav('p110219_m125706_KM5KM6_0c_fvv0_nav');
nav0_EP=read_nav('p110218_m100213_epica_0c_svv0_nav');
sbi1_KM = load('p110219_m125706_KM5KM6_1c_fvv0_sbi');
sbi1_EP = load('p110218_m100213_epica_1c_svv0_sbi');
sbi0_EP = load('p110218_m100213_epica_0c_svv0_sbi');

%Variable definitions
c_light=physconst('LightSpeed');
lambda=c_light/p1_KM.Fc;
max_thet = lambda/(4*p1_KM.AzSpacing); %Same for EP/KM 
k=0:p1_KM.Nra;
el=p1_KM.ZOff+(k-1)*(p1_KM.ZSpacing);
dist = (0:p1_KM.AzSpacing:p1_KM.Naz*p1_KM.AzSpacing);
dist_EP = (0:p1_EP.AzSpacing:p1_EP.Naz*p1_EP.AzSpacing);
PRF_KM = p1_KM.PRF; 
PRF_EP = p1_EP.PRF; 
vel_KM0 = sqrt(nav0_KM.ve.^2+nav0_KM.vn.^2);
vel_EP0 = sqrt(nav0_EP.ve.^2+nav0_EP.vn.^2);

%Nomenclature
%EP: EPICA data
%KM: KM5KM6 data
%0,1: Level 0,1 
%%
%Plot the data
subplot(2,2,1)
generate_Dataplot(dat1_KM,p1_KM,1,'KM');
subplot(2,2,2)
generate_Dataplot(dat1_EP,p1_EP,1,'EP');
subplot(2,2,3)
generate_Dataplot(dat0_KM,p0_KM,0,'KM');
subplot(2,2,4)
generate_Dataplot(dat0_EP,p0_EP,0,'EP');

%Replace corrupt values in altitude and dif parameters for LVL1 KM data
nav1_KM.dif(1381)=nav1_KM.dif(1380);nav1_KM.dif(1382)=nav1_KM.dif(1380);
nav1_KM.alt(1381)=nav1_KM.alt(1380);nav1_KM.alt(1382)=nav1_KM.alt(1380);
 
%Get surface range cells LEVEL 1 
offset=0;
%Range cells KM1
Ra_cells_KM1=round((p1_KM.ZOff-offset-(nav1_KM.alt-nav1_KM.dif))/abs(p1_KM.ZSpacing)).';
%Range cells KM1 from SBI
Ra_cells_KM1_SBI = sbi1_KM(:,2);
%Replace zero data in Range cells KM1 (SBI)
Ra_cells_KM1_SBI(1:88)=60;Ra_cells_KM1_SBI(20438:p1_KM.Naz)=137; 
%Range cells EP
Ra_cells_EP1=round((p1_EP.ZOff-offset-(nav1_EP.alt-nav1_EP.dif))/abs(p1_EP.ZSpacing)).';
%Range cells EP from SBI
Ra_cells_EP1_SBI = sbi1_EP(:,2); Ra_cells_EP1_SBI(Ra_cells_EP1_SBI==0)=NaN;
%Bottom range cells KM1
Ra_cells_KM1_Bottom = sbi1_KM(:,3);
%Put zero data to NaN
Ra_cells_KM1_Bottom(Ra_cells_KM1_Bottom==0)=NaN; 

%Get surface range cells (level 0c)
%Range spacings
RaSpacing_KM = 1/2*c_light*1/p0_KM.Fs;
RaSpacing_EP = 1/2*c_light*(1/p0_EP.Fs);

%First sample delay (level 0c)
Rx_Del_KM=1/2*c_light*p0_KM.RxDelay;
Rx_Del_EP=1/2*c_light*p0_EP.RxDelay;

%Surface range cells KM0/EP0
Ra_cells_KM0 = round((nav0_KM.dif-Rx_Del_KM)/RaSpacing_KM)-2;
Ra_cells_EP0 = round((nav0_EP.dif-Rx_Del_EP)/RaSpacing_EP); %4 sample offset from sbi
Ra_cells_SBI_EP0 = sbi0_EP(:,2); 

%Bottom range cells EP0 from SBI
Ra_cells_EP0_Bottom_SBI = sbi0_EP(:,3);

%Obtain data from known range cells (Surface/Basal data)
%Data is obtained from dif,alt parameters and SBI file
%For KM data level 1
for i = 1:p1_KM.Naz
    
    %Case for parametric extraction
    cell_EX_KM1 = Ra_cells_KM1(i);
    %Case for SBI file
    cell_SBI_KM1 = Ra_cells_KM1_SBI(i);
    cell_SBI_Bottom_KM1 = Ra_cells_KM1_Bottom(i);
    
    %Get surface data
    surf_KM1_EX(i)  = dat1_KM(cell_EX_KM1,i); 
    surf_KM1_SBI(i) = dat1_KM(cell_SBI_KM1,i);
    
    %Get bottom data
    if isnan(cell_SBI_Bottom_KM1)
        bottom_KM1_SBI(i) = NaN;
        else
        bottom_KM1_SBI(i) = dat1_KM(cell_SBI_Bottom_KM1,i);  
    end
    
end

%FOR KM Level 0 
for i=1:p0_KM.Naz
    
    %Extracted surface range cells 
    cell_EX_KM0 =  Ra_cells_KM0(i);
    
    %Get surface data
    surf_KM0_EX(i)  = dat0_KM(cell_EX_KM0,i); 
end

%For EPICA Level 1
for i =1:p1_EP.Naz
     %Extracted surface range cells 
     cell_EX_EP1 =  Ra_cells_EP1(i);
     %SBI case
     cell_SBI_EP1= Ra_cells_EP1_SBI(i);
    
    %Get surface data 
    surf_EP1_EX(i)  = dat1_EP(cell_EX_EP1,i);  
    
    %NaN treatment in SBI
    if isnan(cell_SBI_EP1)
        surf_EP1_SBI(i) = NaN;
    else
        surf_EP1_SBI(i) = dat1_EP(cell_SBI_EP1,i);  
    end 
end

%For EPICA Level 0
for i = 1:p0_EP.Naz
    %Extracted surface range cells 
    cell_EX_EP0 = Ra_cells_EP0(i);
    
    %Get surface data 
    surf_EP0_EX(i)  = dat0_EP(cell_EX_EP0,i);  

end

%%
%Visualize scene (KM1) and get coordinates
factor = 1;
coord=visualizeSceneKM1(dat1_KM,Ra_cells_KM1,factor,1);
coord_KM0 = [3165 17230; 32040 46105; 79100 93165  ]; %Along track sample numbers for parts on shelf
coord_EP1 = [1 1000; 3500 4500; 6400 7400 ]; %Along track sample numbers of parts in EPICA SCENE

%CASE 1: Following the surface
ice_rise1=surf_KM1_EX(coord(1,1):coord(1,2));
ice_rise2=surf_KM1_EX(coord(2,1):coord(2,2));
ice_shelf=surf_KM1_EX(coord(3,1):coord(3,2));
ice_bottom_KM1 = bottom_KM1_SBI(coord(3,1):coord(3,2));
ice_bottom_part2_KM1 =  bottom_KM1_SBI(9728:13300);
ice_rise1_KM0 = surf_KM0_EX(coord_KM0(1,1):coord_KM0(1,2));
ice_rise2_KM0 = surf_KM0_EX(coord_KM0(2,1):coord_KM0(2,2));
ice_shelf_KM0 = surf_KM0_EX(coord_KM0(3,1):coord_KM0(3,2));
ice_complete=surf_KM1_EX(1:end);
ice_complete_KM0 = surf_KM0_EX(1:end-3);
ice_complete_EP1 = surf_EP1_EX(1:end);
ice_EP1 = surf_EP1_EX(coord_EP1(1,1):coord_EP1(1,2));
ice_EP2 = surf_EP1_EX(coord_EP1(2,1):coord_EP1(2,2));
ice_EP3 = surf_EP1_EX(coord_EP1(3,1):coord_EP1(3,2));

%Case 2: Constant range extraction
[const_EX,mismatch] = constantExtraction(128,Ra_cells_KM1,dat1_KM,p1_KM);
[const_EX_KM0,mismatch_KM0] = constantExtraction(128,Ra_cells_KM0,dat0_KM,p0_KM);
ice_rise1_const = const_EX(round(mismatch*coord(1,1):mismatch*coord(1,2)));
ice_rise2_const = const_EX(round(mismatch*coord(2,1):mismatch*coord(2,2)));
ice_shelf_const = const_EX(round(mismatch*coord(3,1):mismatch*coord(3,2)));
ice_rise1_const_KM0 = const_EX_KM0(round(mismatch_KM0*coord_KM0(1,1):mismatch_KM0*coord_KM0(1,2)));
ice_rise2_const_KM0 = const_EX_KM0(round(mismatch_KM0*coord_KM0(2,1):mismatch_KM0*coord_KM0(2,2)));
ice_shelf_const_KM0 = const_EX_KM0(round(mismatch_KM0*coord_KM0(3,1):mismatch_KM0*coord_KM0(3,2)));

%Welch estimates
winlen = 128;overlap=50;
%Estimates for exact extraction
PSD_rise1 = PSD_Welch2(ice_rise1,overlap,winlen);
PSD_rise11= PSD_Welch2(ice_rise1,overlap,512); %Window 512 
PSD_rise2 = PSD_Welch2(ice_rise2,overlap,winlen);
PSD_shelf = PSD_Welch2(ice_shelf,overlap,winlen);
PSD_bottom= PSD_Welch2(ice_bottom_KM1,50,512);
PSD_bottom2= PSD_Welch2(ice_bottom_part2_KM1,50,512);
PSD_rise1_KM0 = PSD_Welch2(ice_rise1_KM0,overlap,winlen);
PSD_rise2_KM0 = PSD_Welch2(ice_rise2_KM0,overlap,winlen);
PSD_shelf_KM0 = PSD_Welch2(ice_shelf_KM0,overlap,winlen);
PSD_EP1 = PSD_Welch2(ice_EP1,overlap,winlen);
PSD_EP2 = PSD_Welch2(ice_EP2,overlap,winlen);
PSD_EP3 = PSD_Welch2(ice_EP3,overlap,winlen);
PSD_EP_comlete = PSD_Welch2(ice_complete_EP1,overlap,winlen);

%Estimates for constant range extractions
PSD_rise1_const = PSD_Welch2(ice_rise1_const,overlap,winlen);
PSD_rise2_const = PSD_Welch2(ice_rise2_const,overlap,winlen);
PSD_shelf_const = PSD_Welch2(ice_shelf_const,overlap,winlen);
PSD_bot = PSD_Welch2(bottom_KM1_SBI(14670:19970),overlap,winlen); %Basal return
PSD_rise1_const_KM0 = PSD_Welch2(ice_rise1_const_KM0,overlap,winlen);
PSD_rise2_const_KM0 = PSD_Welch2(ice_rise2_const_KM0,overlap,winlen);
PSD_shelf_const_KM0 = PSD_Welch2(ice_shelf_const_KM0,overlap,winlen);

%%
close all
%PLOTS OF RESULTS:
figure(1)
plotEstimate1(PSD_rise1,'-.b');
title('PSD Estimates for Fimbul Ice shelf')
hold on
plotEstimate1(PSD_rise2,'-.g');
hold on
plotEstimate1(PSD_shelf,'r');
legend('Positive Doppler','Negative Doppler','Zero Doppler')

figure(2)
plotEstimate0(PSD_rise1_KM0,'-.b',p0_KM.PRF,vel_KM0,[coord_KM0(1,1) coord_KM0(1,2)]);
hold on
plotEstimate0(PSD_rise2_KM0,'-.g',p0_KM.PRF,vel_KM0,[coord_KM0(2,1) coord_KM0(2,2)]);
hold on
plotEstimate0(PSD_shelf_KM0,'r',p0_KM.PRF,vel_KM0,[coord_KM0(3,1) coord_KM0(3,2)]);
title('Level 0c PSD Estimates for Fimbul Ice shelf');
legend('Positive inclination','Negative inclination','Flat part');
set(gca,'FontSize',20)
%%
close all
figure(1)
plotEstimate1(PSD_rise11,'b');
title('Ice rise (positive), window length 512')
set(gca,'FontSize',20)

figure(2)
plotEstimate1(PSD_rise1,'b');
hold on
plotEstimate1(PSD_rise1_const,'g');
title('PSD Fimbul ice shelf, area 1');legend('Exact extraction','Constant range extraction')
figure(3)
plotEstimate1(PSD_rise2,'b');
hold on
plotEstimate1(PSD_rise2_const,'g');
title('PSD Fimbul ice shelf, area 2');legend('Exact extraction','Constant range extraction')
figure(4)
plotEstimate1(PSD_shelf,'b');
hold on
plotEstimate1(PSD_shelf_const,'g');
title('PSD Fimbul ice shelf, area 3');legend('Exact extraction','Constant range extraction')

figure(5)
plotEstimate1(PSD_bottom,'b');
title('PSD Estimate of Fimbul bottom, area 3');

figure(6) 
plotEstimate0(PSD_rise1_const_KM0,'-.b',p0_KM.PRF,vel_KM0,[coord_KM0(1,1),coord_KM0(1,2)]);
hold on
plotEstimate0(PSD_rise2_const_KM0,'-.g',p0_KM.PRF,vel_KM0,[coord_KM0(2,1),coord_KM0(2,2)]);
hold on
plotEstimate0(PSD_shelf_const_KM0,'r',p0_KM.PRF,vel_KM0,[coord_KM0(3,1),coord_KM0(3,2)]);
title('0C Constant Extraction Estimates for Fimbul Ice shelf');
legend('Positive inclination','Negative inclination','Flat part');

figure(7)
hold on
plotEstimate0(PSD_rise2_KM0,'b',p0_KM.PRF,vel_KM0,[coord_KM0(3,1) coord_KM0(3,2)]);
title('Level 0 PSD')

%%
close all
plotEstimate1(PSD_EP1,'r')
hold on
plotEstimate1(PSD_EP2,'g')
hold on
plotEstimate1(PSD_EP3,'b')
title('Level 1 EPICA PSD Estimates');
d1 = [dist_EP(coord_EP1(1,1)) dist_EP(coord_EP1(1,2))];
d2 = [dist_EP(coord_EP1(2,1)) dist_EP(coord_EP1(2,2))];
d3 = [dist_EP(coord_EP1(3,1)) dist_EP(coord_EP1(3,2))];
legend('Along track: 0 to 1.5 km','Along track: 5.2 to 6.7 km','Along track: 9.6 to 11.1 km');
figure
plotEstimate1(PSD_EP_comlete,'b')
title('PSD Estimate of complete EPICA Scene');

%%
close all
%Scattering patterns (takes a little while to run with the defined ranges, in line 285 and 286)
%Azimuth and range areas for analysis (samples)
area_KM1=[coord(3,1) coord(3,2)];range_KM1=[10 -30];
area_KM0 = [coord_KM0(3,1) coord_KM0(3,2)];range_KM0 = [25 -350]; 
area_EP1=[coord_EP1(3,1) coord_EP1(3,2)];range_EP1=[10 -250];range_EP12=[5 -18];
area_EP0 = [1 1000];range_EP0 = [5 -18];

%Plot the scattering patterns
Scatteringpattern1(dat1_KM,p1_KM,nav1_KM,area_KM1,range_KM1,128,'yes')
figure
Scatteringpattern0(dat0_KM,p0_KM,nav0_KM,area_KM0,range_KM0,128,'yes')
figure
Scatteringpattern1(dat1_EP,p1_EP,nav1_EP,area_EP1,range_EP12,128,'yes')
figure
Scatteringpattern0(dat0_EP,p0_EP,nav0_EP,area_EP0,range_EP0,128,'yes')

%%
%SPECTROGRAM GENERATION
%Spectrogram generation 

rise1mat = reshape(ice_complete(89:20088),20,1000).*rectwin(20); %end:20436 start:89, len=20356, (28,727)
spec_rise1mat = fftshift(fft(rise1mat,2^16)/2^16);
theta=linspace(-max_thet,max_thet,size(spec_rise1mat,1))*180/pi;
Dx = (0:p1_KM.AzSpacing:p1_KM.Naz*p1_KM.AzSpacing-1);
Dxmat = reshape(Dx(89:20088),20,1000);
dist_vec = Dxmat(1,:);

%Convert image
figure(1)
scale = 1.2;pow=0.8;
Realimage = abs(spec_rise1mat);
Maximum = 20*mean(mean(Realimage))/scale;
ImageByte = 128 * ( ((Realimage*255)/(Maximum*128)).^pow );
ImageByte=(circshift(ImageByte,-500,2));
image(dist_vec/1000,[theta],(ImageByte))
title('Spectrogram Complete Surface')
ylabel('Look angle [\theta]');xlabel('Along track distance [km]')

figure(2)
generate_Dataplot(dat1_KM,p1_KM,1,'KM');

figure(3) %Bottom KM1
Dx = (0:p1_KM.AzSpacing:p1_KM.Naz*p1_KM.AzSpacing-1);
Dxmat = reshape(Dx(coord(3,1):coord(3,2)),40,100);
dist_vec = Dxmat(1,:);
bottom_KM1 = reshape(ice_bottom_KM1,40,100);
spec_KM1 = fftshift(fft(bottom_KM1,2^16)/2^16);
theta=linspace(-max_thet,max_thet,size(spec_KM1,1))*180/pi;
scale = 2;pow=0.8;
Realimage = abs(spec_KM1);
Maximum = 20*mean(mean(Realimage))/scale;
ImageByte = 128 * ( ((Realimage*255)/(Maximum*128)).^pow );
ImageByte=(circshift(ImageByte,-50,2));
image(dist_vec/1000,theta,ImageByte);
title('Spectrogram Fimbul: Bottom of Ice');
xlabel('Along track distance [km]');ylabel('Look angle \theta [deg]');

figure(4) %Shelf 0C 
complete_KM0 = reshape(ice_complete_KM0,61,1669);
spec_KM0_comp = fftshift(fft(complete_KM0,2^16)/2^16);
scale = 1;pow=0.8;
Realimage = abs(spec_KM0_comp);
Maximum = 20*mean(mean(Realimage))/scale;
ImageByte = 128 * ( ((Realimage*255)/(Maximum*128)).^pow );
ImageByte=(circshift(ImageByte,-835,2));
vel_KM0 = vel_KM0(1:end-3);
vel_mean_KM0=mean(vel_KM0);
df=PRF_KM/length(ice_complete_KM0);
fd=-PRF_KM/2:df:PRF_KM/2-df;
theta = (asin(fd.*lambda/(2*vel_mean_KM0)))*180/pi;
time  = (1/p0_KM.PRF).*(1:p0_KM.Naz-3);
time_mat = reshape(time,61,1669);
time_vec = time_mat(1,:); 
%theta=linspace(-max_thet,max_thet,size(spec_KM0_comp,1))*180/pi;
image(time_vec,[theta],(ImageByte));
ylabel('Look angle \theta [deg]');xlabel('Time along track [s]')
title('Level 0c Fimbul Ice Shelf Spectrogram');

figure(5) %Complete 1C Epica
Dx = (0:p1_EP.AzSpacing:p1_EP.Naz*p1_EP.AzSpacing-2);
Dxmat = reshape(Dx(1:end),39,190);
dist_vec = Dxmat(1,:);
bottom_KM1 = reshape(ice_complete_EP1(1:end-1),39,190);
spec_KM1 = fftshift(fft(bottom_KM1,2^16)/2^16);
theta=linspace(-max_thet,max_thet,size(spec_KM1,1))*180/pi;
scale = 1;pow=0.8;
Realimage = abs(spec_KM1);
Maximum = 20*mean(mean(Realimage))/scale;
ImageByte = 128 * ( ((Realimage*255)/(Maximum*128)).^pow );
%ImageByte=(circshift(ImageByte,-95,2));
image(dist_vec/1000,theta,ImageByte);
title('Spectrogram Complete EPICA Scene');
xlabel('Along track distance [km]');ylabel('Look angle \theta [deg]');


