
data.sessionIDs = {'m5115-20201203-02',...
    'm5115-20201203-01',...
    'm5115-20201202-01',...
    'm5114-20201203-01',...
    'm5114-20201202-02',...
    'm5114-20201202-01',...
    'm5108-20201120-03',...
    'm5108-20201120-02',...
    'm5108-20201120-01'};

dirct = '';
binGapsPos = 0:1.5:157;

for s = 1:length(data.sessionIDs)

    mData(s).sData = load(fullfile(dirct,strcat(data.sessionIDs{s},'.mat')));
    mData(s).sData.behavior.wheelLapDsBinned=mData(s).sData.trials.trialLength-2;
    mData(s).sData.behavior.wheelLapDs=mData(s).sData.trials.trialLength-2;
    for t = 1:length(mData(s).sData.behavior.wheelPosDs)
        [~,ind] = min(abs(binGapsPos-mData(s).sData.behavior.wheelPosDs(t)));
        mData(s).sData.behavior.wheelPosDsBinned(t) = ind;
    end
    
    
    mData(s).rmaps = sb.generate.rasterMapsForSession(mData(s).sData,"dff_signal", 'dff','manipulated_motor',0);
      
end
