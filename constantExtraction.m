function [surface_const,mismatch] = constantExtraction(segLen,cells,data,pstruct)
%Output: Constant extraction data, mismatch to compare plots
%Mismatch value: Due to the constant extraction, these plots are double
%length compared to the exact extraction case, therefore mismatch value is
%used to correct for dimensions when comparing constant/exact
%Inputs:
%1) Length of segments
%2) Range cells for surface return
%3) Data file
%4) p structure

%Constant range surface data extraction
Ra_cells_surf = cells;
start = 1;
stop = segLen;
%Create matrix (Az) that contains the azimuth sample numbers for each of
%the constant range segments

 for i = 1:2*floor(pstruct.Naz/segLen)-1
     temp = start:stop;
     Az(i,:) = temp;
     midval = floor(median(temp)); 
     start = midval+1;
     stop=midval+segLen;   
 end

%Reshape the matrix into a vector 
Az1 = reshape(Az.',1,numel(Az));

%Get the surface range indicies for all the azimuth samples
idx_EX = Ra_cells_surf(Az1);

%Number of azimuth samples for constant extraction is not identical to the case
%where the surface is traced. Therefore mismatch is used for proper
%comparison
mismatch = length(idx_EX)/pstruct.Naz;

%Get the surface data for the case of constant extraction: The mean range
%cell (avg_range) for a segment is used as the representative range value.

for i=1:2*round(pstruct.Naz/segLen)-1 %round changed from floor
    
    avg_range(i) = round(mean(idx_EX(round(mismatch*min(Az(i,:))):round(mismatch*max(Az(i,:))))));
  
    %Pick out surface from the data, given the average range cell in a
    %given azimuth segment
    surface_const(:,i) = data(avg_range(i),min(Az(i,:)):max(Az(i,:)));
  
end

%Reshaping data matrix to a vector
surface_const=reshape(surface_const,1,numel(surface_const));
end

