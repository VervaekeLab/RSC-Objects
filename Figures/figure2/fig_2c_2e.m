%% Find the number of responses in each category and the intersections between them
E = 0; % Enganged responsive
S = 0; % Passive Still responsive
R = 0; % Passive Running responsive

EiS = 0; % intersection between E and S
EiR = 0; % ...
SiR = 0;
EiSiR = 0; % ...

% Getting the numbers per session
S_session = zeros(1,7);
R_session = zeros(1,7);
E_session = zeros(1,7);

for i = 1:length(analysis_venn)
    
    % Engaged cells
    if analysis_venn(i).isLcs
        
        E = E+1;
        E_session(analysis_venn(i).session) = E_session(analysis_venn(i).session) + 1;

        if analysis_venn(i).still.response
            EiS = EiS+1;
            S_session(analysis_venn(i).session) = S_session(analysis_venn(i).session) + 1;
            
            if analysis_venn(i).moving.response
                EiSiR = EiSiR+1;
            end
            
        end
        
        if analysis_venn(i).moving.response
            EiR = EiR+1;
            R_session(analysis_venn(i).session) = R_session(analysis_venn(i).session) + 1;

        end
    end
    
    % Still
    if analysis_venn(i).still.response
        S = S+1;
        
        if analysis_venn(i).moving.response
            SiR = SiR+1;
        end
        
    end
    
    % Running
    if analysis_venn(i).moving.response
       R = R+1; 
    end

    
end

% Every other session is 0 because it is not a passive stimulus session,
% discard these
S_session = S_session([1,3,5,7]);
R_session = R_session([1,3,5,7]);
E_session = E_session([1,3,5,7]);


%% Create venn diagram
landmark_color = [239,35,60]/255;
    
figure(5); clf;
subplot(1,3,1);

cmap = cbrewer('div','RdYlBu',7)';

e_color = cmap(:,1);
s_color = cmap(:,3);
r_color = cmap(:,7);

venn([E,S,R],[EiS,EiR,SiR,EiSiR],'FaceColor',{e_color,s_color,r_color});
xticks([])
yticks([])
title('Venn diagram of responses')
legend({'Engaged','Passive (immobile)', 'Passive (running)'})

hold on
text(-4,-1.4,sprintf('%i\n(%.1f%%)',E,([E]/length(analysis_venn))*100),'FontSize',20,'Color','w')
text(1.4,8.3,sprintf('%i\n(%.1f%%)',R,([R]/length(analysis_venn))*100),'FontSize',20,'Color','w')
text(7.4,-2,sprintf('%i\n(%.1f%%)',S,([S]/length(analysis_venn))*100),'FontSize',20,'Color','k')
set(gca,'FontSize',20)

%% Plot the response distribtuion in pasive stimulation for landmark cells

session_lcs_any = zeros(7,1);
session_lcs_all = zeros(7,1);
session_never = zeros(7,1);
session_cells= zeros(7,1);
for i = 1:length(analysis_venn)
     if analysis_venn(i).isLcs || analysis_venn(i).moving.response || analysis_venn(i).still.response
         session_lcs_any(analysis_venn(i).session)= session_lcs_any(analysis_venn(i).session)+1;
     end
     if ~analysis_venn(i).isLcs && ~analysis_venn(i).moving.response && ~analysis_venn(i).still.response
         session_never(analysis_venn(i).session)= session_never(analysis_venn(i).session)+1;
     end
     
     session_cells(analysis_venn(i).session) = session_cells(analysis_venn(i).session)+1;
     
end

session_cells(2:2:end) = [];
session_lcs_any(2:2:end) = [];
session_never(2:2:end) = [];
session_lcs_all(2:2:end) = [];



subplot(1,3,2);
bar(1,nanmean(100*session_lcs_any./session_cells))
hold on
errorbar(1,nanmean(100*session_lcs_any./session_cells),nanstd(100*session_lcs_any./session_cells)/sqrt(length(session_cells)),'k','LineWidth',1.5)
scatter(ones(4,1),100*session_lcs_any./session_cells,70,'k')
ylabel('Cells (%)')

bar(2,nanmean(100*session_never./session_cells))
hold on
errorbar(2,nanmean(100*session_never./session_cells),nanstd(100*session_never./session_cells)/sqrt(length(session_cells)),'k','LineWidth',1.5)
scatter(2*ones(4,1),100*session_never./session_cells,70,'k')
ylabel('Cells (%)')
xticks([1:2])
xticklabels({'Any','Never'})

subplot(1,3,3);
pie([sum(session_lcs_any)/sum(sum(session_cells)),
sum(session_never)/sum(sum(session_cells))])
legend({'Any','Never'})



figure()
% cla;

b = bar(1,100);
hold on
b(1).FaceColor = e_color;
xticks([])

b = bar(2,(EiR/E)*100);
b(1).FaceColor = r_color;

b = bar(3,(EiS/E)*100);
b(1).FaceColor = s_color;


% Add errorbars
R_low_error = -(std(R_session./E_session)/sqrt(length(R_session./E_session)))*100;
R_high_error = (std(R_session./E_session)/sqrt(length(R_session./E_session)))*100;
er = errorbar(2,(EiR/E)*100,R_low_error, R_high_error);
er.Color = [0,0,0];
er.LineStyle = 'none';
er.LineWidth = 2;

S_low_error = -(std(S_session./E_session)/sqrt(length(S_session./E_session)))*100;
S_high_error = (std(S_session./E_session)/sqrt(length(S_session./E_session)))*100;
er = errorbar(3,(EiS/E)*100,S_low_error, S_high_error);
er.Color = [0,0,0];
er.LineStyle = 'none';
er.LineWidth = 2;

ylim([0,110])
yticks([0,20,40,60,80,100])
ylabel('% of landmark cells')

title('Response distribution of landmark cells')
set(gca,'FontSize',20)

% Set correct size of the figure
% set(gcf,'renderer', 'painters', 'Position', [2200,100,1000,400])

% 
% 
