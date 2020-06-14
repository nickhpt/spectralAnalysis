function generate_Dataplot(data,pstruct,level,scene)
%Input
%1: Data
%2: Navigation file structure
%3: Data level: 1 or 0
%4: Scene: 'EP' or 'KM'
%Example Level 1: generate_Dataplot(dat1_EP,p1_EP,1,'EP');
%Output
%Void function generating desired plot

%Convert image
scale = 1;pow=0.4;
Realimage = abs(data);
Maximum = 20*mean(mean(Realimage))/scale;
ImageByte = 128 * ( ((Realimage*255)/(Maximum*128)).^pow );

if contains(scene,"KM"); type = 1; else; type = 0; end
    
switch level   
    case 1 %Case for level 1 data
        if type==1
            p1_KM = pstruct;
            title_string = 'Level 1: KM5KM6 Data Scene';
            k=0:p1_KM.Nra;
            el=p1_KM.ZOff+(k-1)*(p1_KM.ZSpacing);
            dist = (0:p1_KM.AzSpacing:p1_KM.Naz*p1_KM.AzSpacing);
            
            image(dist/1000,el,ImageByte)
            ax = gca;ax.YDir = 'normal';
            xlabel('Along-track distance [km]');ylabel('Elevation [m]');
            title(title_string);  
            
        else %Else generate EPICA lvl 1       
            p1_EP = pstruct;
            title_string = 'Level 1: EPICA Data Scene';
            k=0:p1_EP.Nra;
            el=p1_EP.ZOff+(k-1)*(p1_EP.ZSpacing);
            dist = (0:p1_EP.AzSpacing:p1_EP.Naz*p1_EP.AzSpacing);
            
            image(dist/1000,el/1000,ImageByte)
            ax = gca;ax.YDir = 'normal';
            xlabel('Along-track distance [km]');ylabel('Elevation [km]');
            title(title_string);  
        end   
    
    case 0 %Case for level 0 data
        if type == 1
            p0_KM = pstruct;
            title_string = 'Level 0: KM5KM6 Data Scene';
            delay = [p0_KM.RxDelay p0_KM.RxDelay+(2:p0_KM.Nra)*1/p0_KM.Fs];
            time  = (1/p0_KM.PRF).*(1:p0_KM.Naz);
            image(time,delay*10^6,ImageByte);
            xlabel('Time along track [s]');ylabel('Delay [us]');
            title(title_string);  
                    
        else %Else generate EPICA lvl 0 
            p0_EP = pstruct;
            title_string = 'Level 0: EPICA Data Scene';
            delay = [p0_EP.RxDelay p0_EP.RxDelay+(2:p0_EP.Nra)*1/p0_EP.Fs];
            time  = (1/p0_EP.PRF).*(1:p0_EP.Naz);
            image(time,delay*(10^6),ImageByte);
            xlabel('Time along track [s]');ylabel('Delay [us]');
            title(title_string);  
            
        end      
end
end

