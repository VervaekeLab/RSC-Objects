
% Plot one example in each category
figure(3); clf;

shifted_color = [0.5,0.5,0.7];
shifted_lineStyle = '-';

% Biased to first landmark 1
plot(normalize(avg_tuning_map_on(sorted_amplitude_change_index(8),:),'range',[0,1])-1,'Color','k'); 
hold on;
%plot(circshift(normalize(avg_tuning_map_on(sorted_amplitude_change_index(8),:),'range',[0,1])-1,40),'Color',shifted_color,'LineStyle','-'); 

% Biased to first landmark 2
plot(normalize(avg_tuning_map_on(sorted_amplitude_change_index(20),:),'range',[0,1])-2,'Color','k'); 
hold on;
%plot(circshift(normalize(avg_tuning_map_on(sorted_amplitude_change_index(20),:),'range',[0,1])-2,40),'Color',shifted_color,'LineStyle','-'); 

% Equally biased
plot(normalize(avg_tuning_map_on(sorted_amplitude_change_index(203),:),'range',[0,1])-3,'Color','k'); 
%plot(circshift(normalize(avg_tuning_map_on(sorted_amplitude_change_index(203),:),'range',[0,1])-3,100),'Color',repmat(0.65,1,3),'LineStyle','--'); 

% Biased to second landmark 1
plot(normalize(avg_tuning_map_on(sorted_amplitude_change_index(344),:),'range',[0,1])-4,'Color','k'); 
%plot(circshift(normalize(avg_tuning_map_on(sorted_amplitude_change_index(344),:),'range',[0,1])-4,-40),'Color',shifted_color,'LineStyle','-'); 

% Biased to second landmark 2
% plot(normalize(avg_tuning_map_on(sorted_amplitude_change_index(371),:),'range',[0,1])-5,'Color','k'); 
%plot(circshift(normalize(avg_tuning_map_on(sorted_amplitude_change_index(364),:),'range',[0,1])-5,-40),'Color',shifted_color,'LineStyle','-'); 

% plot(normalize(avg_tuning_map_on(sorted_amplitude_change_index(417),:),'range',[0,1])-6,'Color','k'); 

plot(normalize(avg_tuning_map_on(sorted_amplitude_change_index(474),:),'range',[0,1])-5,'Color','k'); 

% plot(normalize(avg_tuning_map_on(sorted_amplitude_change_index(486),:),'range',[0,1])-8,'Color','k'); 

plot(normalize(avg_tuning_map_on(sorted_amplitude_change_index(528),:),'range',[0,1])-6,'Color','k'); 
% plot(normalize(avg_tuning_map_on(sorted_amplitude_change_index(538),:),'range',[0,1])-6,'Color','k'); 
% plot(normalize(avg_tuning_map_on(sorted_amplitude_change_index(542),:),'range',[0,1])-6,'Color','k'); 
% plot(normalize(avg_tuning_map_on(sorted_amplitude_change_index(551),:),'range',[0,1])-10,'Color','k'); 
% % plot(normalize(avg_tuning_map_on(sorted_amplitude_change_index(560),:),'range',[0,1])-6,'Color','k'); 
% plot(normalize(avg_tuning_map_on(sorted_amplitude_change_index(561),:),'range',[0,1])-11,'Color','k'); 
% plot(normalize(avg_tuning_map_on(sorted_amplitude_change_index(573),:),'range',[0,1])-12,'Color','k'); 
% % plot(normalize(avg_tuning_map_on(sorted_amplitude_change_index(582),:),'range',[0,1])-6,'Color','k'); 
% plot(normalize(avg_tuning_map_on(sorted_amplitude_change_index(590),:),'range',[0,1])-13,'Color','k'); 
plot(normalize(avg_tuning_map_on(sorted_amplitude_change_index(591),:),'range',[0,1])-7,'Color','k'); 
% plot(normalize(avg_tuning_map_on(sorted_amplitude_change_index(596),:),'range',[0,1])-15,'Color','k'); 

% Draw lines to indicate the landmark onsets
line([32,32],[-10,10],'Color','r')
line([72,72],[-10,10],'Color','r')
line([92,92],[-10,10],'Color','k')
yticks({})
ylim([-7.1,0.1])
xlim([0,105])
xticks([0,105]);
xticklabels({0,157});
xlabel('Position (cm)');
ylabel('Cell #');
title('Example cells');
set(gca,'FontSize',16)
box off

% Set correct size of the figure
set(gcf,'renderer', 'painters', 'Position', [-1500,100,600,500])



% % Plot most biased in both directions
% n_examples = 5;
% figure(3); clf;
% subplot(1,2,1);
% for i = 1:n_examples
%    plot(normalize(avg_tuning_map_on(sorted_amplitude_change_index(i),:),'range',[0,1])-i,'Color','k'); 
%    hold on;
%    plot(circshift(normalize(avg_tuning_map_on(sorted_amplitude_change_index(i),:),'range',[0,1])-i,40),'Color',repmat(0.65,1,3),'LineStyle','--'); 
% 
% end
% 
% % Draw lines to indicate landmarks
% line([32,32],[-100,100],'Color','r')
% line([72,72],[-100,100],'Color','r')
% line([92,92],[-100,100],'Color','k')
% 
% ylim([-n_examples,0]);
% xlim([0,105])
% 
% 
% subplot(1,2,2);
% for i = 1:n_examples
%    plot(normalize(avg_tuning_map_on(sorted_amplitude_change_index(end-i),:),'range',[0,1])-i,'Color','k'); 
%     hold on;
%    plot(circshift(normalize(avg_tuning_map_on(sorted_amplitude_change_index(end-i),:),'range',[0,1])-i,-40),'Color',repmat(0.65,1,3),'LineStyle','--'); 
% end
% 
% ylim([-n_examples,0]);
% line([32,32],[-100,100],'Color','r')
% line([72,72],[-100,100],'Color','r')
% line([92,92],[-100,100],'Color','k')
