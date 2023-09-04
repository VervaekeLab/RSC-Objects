function cmap = colorPicker(name)
% Based on the name, returns the correct colormap.

switch name
    case 'placeCells'
        cmap = cbrewer('seq','Purples',3);
        cmap = cmap(3,:);
        
    case 'landmarkCellsNonLasting'
   
        cmap = cbrewer('seq','Greens',3);
        cmap = cmap(3,:);
         
    case 'landmarkCellsLasting'
        
        cmap = cbrewer('seq','Oranges',3);
        cmap = cmap(3,:);
        
    otherwise
        warning('Colormap doesnt exist');
        cmap = [0,0,0];
    
end