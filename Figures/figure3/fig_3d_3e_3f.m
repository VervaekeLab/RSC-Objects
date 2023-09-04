%% Create the classification plot for all spatially tuned cells
% Find all unique place cells
cell_counter = 1;

response_map = zeros(3,5000);    
for s = 1:length(mData)
   
    % Find all unique place cell inds in the current session
    pcs = unique([mData(s).rmaps.level1.pcs,mData(s).rmaps.level2.pcs,mData(s).rmaps.level3.pcs,mData(s).rmaps.level1.lcs,mData(s).rmaps.level2.lcs,mData(s).rmaps.level3.lcs]);
    levels = {'level1','level2','level3'};
    
    for c = 1:length(pcs)
        
        for l = 1:length(levels)
           level = levels{l};
           
           % If the place cell is in the list of pcs at the current light level, add a 1 to the given cell at the given light level.
           if sum(ismember(mData(s).rmaps.(level).pcs,pcs(c)))
              response_map(l,cell_counter) = 1;
           elseif sum(ismember(mData(s).rmaps.(level).lcs,pcs(c)))
              response_map(l,cell_counter) = 2;               
           end
        end
        cell_counter = cell_counter + 1;
    end
end

response_map = response_map(:,1:cell_counter-1);

% Let us sort the responses
response_map = response_map';

sorted_response_map = sortrows(response_map,'descend');
figure(5); clf;
subplot(1,3,1)
imagesc(sorted_response_map,[0,2])
yticks([1,size(sorted_response_map,1)])

ylabel('Cell #')
xticks([1,2,3])
xtickangle(0)
xticklabels({'No light','Dim light','Bright light'})
set(gca,'FontSize',16);
title("Evolution of cell classification")

colormap([1,1,1;place_cells_color;landmark_cells_color]);

%% Probability of changing classification 
% This is the probability that a cell will either loose, maintain or change
% its class (place cell or landmark cell) when the condition changes from
% no light -> dim light or dim light -> bright light.
% 
% subplot(9,1,6:7);
% 
% % Landmark cells
% lcs_flip = [];
% for c = 1:size(response_map,1)
%     
%    % --- From no light -> dim light
%    % Let us look at landmark cells 
%    if response_map(c,1) == 2
% 
%        % If it looses tuning
%        if response_map(c,2) == 0
%             the_flip = [1,0,0];
%        end
% 
%        % If tuning stays
%        if response_map(c,2) == 2
%             the_flip = [0,1,0];
%        end
% 
%        % If tuning changes category
%        if response_map(c,2) == 1
%             the_flip = [0,0,1];
%        end
% 
%        lcs_flip = [lcs_flip; the_flip];
%    end
%    
% %    % --- From dim light -> bright light
% %    % Let us look at landmark cells 
% %    if response_map(c,2) == 2
% % 
% %        % If it looses tuning
% %        if response_map(c,3) == 0
% %             the_flip = [1,0,0];
% %        end
% % 
% %        % If tuning stays
% %        if response_map(c,3) == 2
% %             the_flip = [0,1,0];
% %        end
% % 
% %        % If tuning changes category
% %        if response_map(c,3) == 1
% %             the_flip = [0,0,1];
% %        end
% % 
% %        lcs_flip = [lcs_flip; the_flip];
% %    end
% 
%   % --- From no light -> bright light
%    % Let us look at landmark cells 
%    if response_map(c,1) == 2
% 
%        % If it looses tuning
%        if response_map(c,3) == 0
%             the_flip = [1,0,0];
%        end
% 
%        % If tuning stays
%        if response_map(c,3) == 2
%             the_flip = [0,1,0];
%        end
% 
%        % If tuning changes category
%        if response_map(c,3) == 1
%             the_flip = [0,0,1];
%        end
% 
%        lcs_flip = [lcs_flip; the_flip];
%    end
%     
% end
% lcs_results = sum(lcs_flip) / sum(sum(lcs_flip)) *100;
% 
% % Place cells
% pcs_flip = [];
% for c = 1:size(response_map,1)
%     
%    % --- From no light -> dim light
%    % Let us look at landmark cells 
%    if response_map(c,1) == 1
% 
%        % If it looses tuning
%        if response_map(c,2) == 0
%             the_flip = [1,0,0];
%        end
% 
%        % If tuning stays
%        if response_map(c,2) == 1
%             the_flip = [0,1,0];
%        end
% 
%        % If tuning changes category
%        if response_map(c,2) == 2
%             the_flip = [0,0,1];
%        end
% 
%        pcs_flip = [pcs_flip; the_flip];
%    end
%    
%    % --- From dim light -> bright light
%    % Let us look at landmark cells 
%    if response_map(c,2) == 2
% 
%        % If it looses tuning
%        if response_map(c,3) == 0
%             the_flip = [1,0,0];
%        end
% 
%        % If tuning stays
%        if response_map(c,3) == 1
%             the_flip = [0,1,0];
%        end
% 
%        % If tuning changes category
%        if response_map(c,3) == 2
%             the_flip = [0,0,1];
%        end
% 
%        pcs_flip = [pcs_flip; the_flip];
%    end
% 
% %   % --- From no light -> bright light
% %    % Let us look at landmark cells 
% %    if response_map(c,1) == 2
% % 
% %        % If it looses tuning
% %        if response_map(c,3) == 0
% %             the_flip = [1,0,0];
% %        end
% % 
% %        % If tuning stays
% %        if response_map(c,3) == 1
% %             the_flip = [0,1,0];
% %        end
% % 
% %        % If tuning changes category
% %        if response_map(c,3) == 2
% %             the_flip = [0,0,1];
% %        end
% % 
% %        pcs_flip = [pcs_flip; the_flip];
% %    end
%     
% end
% 
% pcs_results = sum(pcs_flip) / sum(sum(pcs_flip)) * 100;
% b = bar([lcs_results;pcs_results]');
% b(1).FaceColor = landmark_cells_color;
% b(1).EdgeColor = 'none';
% b(2).FaceColor = place_cells_color;
% b(2).EdgeColor = 'none';
% ylim([0,60])
% xticklabels({'Loose tuning','Maintain tuning','Change tuning'})
% xtickangle(25)
% ylabel('Probability (%)')
% legend({'Landmark cells','Place cells'})
% title('Probability of changing classification');
% box off
% set(gca,'FontSize',16)


%% Probability of keeping the tuning through each context

% Tactile cells
tactile_cells_maintain = 0;
tactile_cells_total = 0;
position_cells_maintain = 0;
position_cells_total = 0;

for c = 1:size(response_map,1)
    
    % Tactile cells
    if response_map(c,1) == 2
        tactile_cells_total = tactile_cells_total + 1;
        
        % Maintain tuning
       if sum(response_map(c,2:3) == [2,2]) == 2
            tactile_cells_maintain = tactile_cells_maintain + 1;
       end
       
    end
    
    % Position tuned cells
    if response_map(c,1) == 1
        position_cells_total = position_cells_total + 1;
        
        % Maintain tuning
        if sum(response_map(c,2:3) == [1,1]) == 2
            position_cells_maintain = position_cells_maintain + 1;
        end
        
    end
    
end


tactile_cells_maintain = (tactile_cells_maintain / tactile_cells_total) * 100;
position_cells_maintain = (position_cells_maintain / position_cells_total) * 100;



subplot(1,3,2);
b = bar([tactile_cells_maintain,position_cells_maintain]);

ylim([0,100])
xticklabels({'Tactile cells','Position-tuned cells'})
xtickangle(25)
ylabel('Probability (%)')
title('Probability of changing classification');
box off
set(gca,'FontSize',16)

set(gcf,'renderer', 'painters', 'Position', [100,100,500,600])

%% Venn diagram
% Create a venn diagram
% Venn diagram of cells that are classified in each of the light conditions
% Init vars
l1 = 0;
l2 = 0;
l3 = 0;
l1l2 = 0;
l1l3 = 0;
l2l3 = 0;
l1l2l3 = 0;

% For each session find the landmark cells that match
for s = 1:length(mData)
    
    l1_cells = mData(s).rmaps.level1.lcs;
    l2_cells = mData(s).rmaps.level2.lcs;
    l3_cells = mData(s).rmaps.level3.lcs;
   
    l1 = l1 + length(l1_cells);
    l2 = l2 + length(l2_cells);
    l3 = l3 + length(l3_cells);
    
    l1l2 = l1l2 + sum(ismember(l1_cells,l2_cells));
    l1l3 = l1l3 + sum(ismember(l1_cells,l3_cells));
    l2l3 = l2l3 + sum(ismember(l2_cells,l3_cells));
    
    l1l2l3 = l1l2l3 + sum(ismember(l1_cells,l2_cells(ismember(l2_cells,l3_cells))));

end

no_all_cells = 0;
for s = 1:length(mData)
    no_all_cells =no_all_cells+size(mData(s).rmaps.dff,3);
end
subplot(1,3,3);
cmap = cbrewer('seq','Greys',8);
venn([l1,l2,l3],[l1l2,l1l3,l2l3,l1l2l3],'FaceColor',{[cmap(8,:)],[cmap(6,:)],[cmap(2,:)]});
xticks([]);
yticks([]);
legend({'Dark','Dim light','Bright light'})
axis('square')
box off

% Set correct size of the figure
set(gcf,'renderer', 'painters', 'Position', [100,100,1500,500])
