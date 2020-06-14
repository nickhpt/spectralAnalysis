function [coordinates] = visualizeSceneKM1(dat1_KM,cells,factor,p)
%Output: Plot of KM1 scene + coordinates of 3 parts (2 on rise, 1 on shelf)
%Input
%1) Data
%2) Range cells for surface return
%3) Factor to decrease length of scenes
%4) Plot wished = 1 else 0

%Save range cells
Ra_cells_surf = cells;
%Original fixed scene size
start_Len = 4000;
%Resulting length after decreasing length
tot_Len = start_Len/factor;
%Correction term for applying factor 
c = start_Len - tot_Len;

%Azimuth coordinates (samples) for parts in ice rise and shelf
co = [713+c/2 4712-c/2;...
               8450+c/2 12449-c/2;...
               15360+c/2 19359-c/2];

%dist = (0:1.5:20542*1.5); 
%co=dist(co);     
coordinates = co;

% k=0:560;
% el=410+(k-1)*(-2.5);
% dist = (0:1.5:20542*1.5);

if p==1
%generate_Dataplot(dat1_KM,p1_KM,1,'KM');
imagesc(abs(log(dat1_KM)));
hold on
plot(Ra_cells_surf,'r','LineWidth',1.5)
hold on
plot([co(1,1) co(1,1)],[0.5 80],'k-.*','LineWidth',1.5)
hold on
plot([co(1,2) co(1,2)],[0.5 65],'k-.*','LineWidth',1.5)
hold on
plot([co(2,1) co(2,1)],[25 100],'k-.*','LineWidth',1.5)
hold on
plot([co(2,2) co(2,2)],[25+40 100+40],'k-.*','LineWidth',1.5)
hold on
plot([co(3,1) co(3,1)],[115 190],'k-.*','LineWidth',1.5)
hold on
plot([co(3,2) co(3,2)],[115 190],'k-.*','LineWidth',1.5)
hold on
set(gcf,'position',[320 342 560 420])
%xlabel('Azimuth samples');ylabel('Range samples');
end           
end

