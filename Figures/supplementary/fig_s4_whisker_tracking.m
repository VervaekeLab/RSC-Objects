
M2=open("D:\GoodTrack\m1412-20191010-01_11-18-19.000DLC_mobnet_100_WhiskerDetectFeb7shuffle1_50000.csv");
M3=open("D:\GoodTrack\m1414-20191023-01_13-15-19.000DLC_mobnet_100_WhiskerDetectFeb7shuffle1_50000.csv");
M4=open("D:\GoodTrack\m1415-20191023-01_13-48-25.000DLC_mobnet_100_WhiskerDetectFeb7shuffle1_50000.csv");
M6=open("D:\GoodTrack\m1408-20190818-01_14-28-47.000DLC_mobnet_100_WhiskerDetectFeb7shuffle1_50000.csv");


M2_data = M2.data;
M3_data = M3.data;
M4_data = M4.data;
M6_data = M6.data;

% Data flipped from video (up is down)

%% error estimate

[Lc2,Le2]=hist(M2_data(:,7),20);
disp(['Removed % :',num2str(sum(Lc2(1:18))/sum(Lc2)),' at  Likelyhood: ',num2str(Le2(18))])
[Rc2,Re2]=hist(M2_data(:,13),20);
disp(['Removed % :',num2str(sum(Rc2(1:17))/sum(Rc2)),' at  Likelyhood: ',num2str(Re2(17))])
disp([' '])

[Lc3,Le3]=hist(M3_data(:,7),20);
disp(['Removed % :',num2str(sum(Lc3(1:8))/sum(Lc3)),' at  Likelyhood: ',num2str(Le3(8))])
[Rc3,Re3]=hist(M3_data(:,13),20);
disp(['Removed % :',num2str(sum(Rc3(1:13))/sum(Rc3)),' at  Likelyhood: ',num2str(Re3(13))])
disp([' '])

[Lc4,Le4]=hist(M4_data(:,7),20);
disp(['Removed % :',num2str(sum(Lc4(1:10))/sum(Lc4)),' at  Likelyhood: ',num2str(Le4(10))])
[Rc4,Re4]=hist(M4_data(:,13),20);
disp(['Removed % :',num2str(sum(Rc4(1:15))/sum(Rc4)),' at  Likelyhood: ',num2str(Re4(15))])


disp([' '])
[Lc4,Le4]=hist(M6_data(:,7),20);
disp(['Removed % :',num2str(sum(Lc4(1:11))/sum(Lc4)),' at  Likelyhood: ',num2str(Le4(10))])
[Rc4,Re4]=hist(M4_data(:,13),20);
disp(['Removed % :',num2str(sum(Rc4(1:15))/sum(Rc4)),' at  Likelyhood: ',num2str(Re4(15))])
%% remove outliers then interpolate missing 

% Remove outliers based on wrong side og mouse

mid=mean([M2_data(:,2);M2_data(:,8)]);   %removes base at wrong side
RemoveIndx=find(M2_data(:,8) < mid);
M2_data(RemoveIndx,:)=nan;
RemoveIndx=find(M2_data(:,11) < mid);    %removes perifery at wrong side
M2_data(RemoveIndx,:)=nan;
disp(['Estimated error M2: ',num2str(numel(RemoveIndx)/numel(M2_data(:,1)))])

mid=mean([M3_data(:,2);M3_data(:,8)]);   %removes base at wrong side
RemoveIndx=find(M3_data(:,8) < mid);
M3_data(RemoveIndx,:)=nan;
RemoveIndx=find(M3_data(:,11) < mid);    %removes perifery at wrong side
M3_data(RemoveIndx,:)=nan;
disp(['Estimated error M3: ',num2str(numel(RemoveIndx)/numel(M3_data(:,1)))])

mid=mean([M4_data(:,2);M4_data(:,8)]);   %removes base  at wrong side
RemoveIndx=find(M4_data(:,8) < mid);
M4_data(RemoveIndx,:)=nan;
RemoveIndx=find(M4_data(:,11) < mid);    %removes perifery at wrong side
M4_data(RemoveIndx,:)=nan;
disp(['Estimated error M4: ',num2str(numel(RemoveIndx)/numel(M4_data(:,1)))])


RemoveIndx=find(M4_data(:,11) < 510);
M4_data(RemoveIndx,:)=nan;
RemoveIndx=find(M3_data(:,11) < 450);
M3_data(RemoveIndx,:)=nan;
RemoveIndx=find(M2_data(:,11) < 450);
M2_data(RemoveIndx,:)=nan;
RemoveIndx=find(M6_data(:,2) < 490)
M6_data(RemoveIndx,:)=nan;




'done'
%% downsample to 60Hz
fs=31;        %sampling frequancy
dM2 = downsample(M2_data,round(299/fs));
dM3 = downsample(M3_data,round(299/fs));
dM4 = downsample(M4_data,round(299/fs));
dM6 = downsample(M6_data,round(299/fs));

disp([num2str(mean(sum(isnan(dM2))/numel(dM2(:,1))))])
disp([num2str(mean(sum(isnan(dM3))/numel(dM3(:,1))))])
disp([num2str(mean(sum(isnan(dM4))/numel(dM4(:,1))))])
disp([num2str(mean(sum(isnan(dM6))/numel(dM6(:,1))))])
%% interpolate missing data

IndFixList=find(isnan(dM2(:,8)));
for I =1:numel(IndFixList)
    dM2(IndFixList(I),8)  = mean([dM2(IndFixList(I)-1,8),dM2(IndFixList(I)+1,8)]);
    dM2(IndFixList(I),11) = mean([dM2(IndFixList(I)-1,11),dM2(IndFixList(I)+1,11)]);
    dM2(IndFixList(I),9) = mean([dM2(IndFixList(I)-1,9),dM2(IndFixList(I)+1,9)]);
    dM2(IndFixList(I),12) = mean([dM2(IndFixList(I)-1,12),dM2(IndFixList(I)+1,12)]);
end
[sum(isnan(dM2(:,8))),sum(isnan(dM2(:,9))),sum(isnan(dM2(:,11))),sum(isnan(dM2(:,12)))] 


IndFixList=find(isnan(dM3(:,8)));
for I =1:numel(IndFixList)
    dM3(IndFixList(I),8)  = mean([dM3(IndFixList(I)-1,8),dM3(IndFixList(I)+1,8)]);
    dM3(IndFixList(I),11) = mean([dM3(IndFixList(I)-1,11),dM3(IndFixList(I)+1,11)]);
    dM3(IndFixList(I),9) = mean([dM3(IndFixList(I)-1,9),dM3(IndFixList(I)+1,9)]);
    dM3(IndFixList(I),12) = mean([dM3(IndFixList(I)-1,12),dM3(IndFixList(I)+1,12)]);
end
[sum(isnan(dM3(:,8))),sum(isnan(dM3(:,9))),sum(isnan(dM3(:,11))),sum(isnan(dM3(:,12)))] 


IndFixList=find(isnan(dM4(:,8)));
for I =1:numel(IndFixList)
    dM4(IndFixList(I),8)  = mean([dM4(IndFixList(I)-1,8),dM4(IndFixList(I)+1,8)]);
    dM4(IndFixList(I),11) = mean([dM4(IndFixList(I)-1,11),dM4(IndFixList(I)+1,11)]);
    dM4(IndFixList(I),9) = mean([dM4(IndFixList(I)-1,9),dM4(IndFixList(I)+1,9)]);
    dM4(IndFixList(I),12) = mean([dM4(IndFixList(I)-1,12),dM4(IndFixList(I)+1,12)]);
end
[sum(isnan(dM4(:,8))),sum(isnan(dM4(:,9))),sum(isnan(dM4(:,11))),sum(isnan(dM4(:,12)))] 

IndFixList=find(isnan(dM6(:,11)));
for I =1:numel(IndFixList)
    dM6(IndFixList(I),8)  = mean([dM6(IndFixList(I)-1,8),dM6(IndFixList(I)+1,8)]);
    dM6(IndFixList(I),11) = mean([dM6(IndFixList(I)-1,11),dM6(IndFixList(I)+1,11)]);
    dM6(IndFixList(I),9) = mean([dM6(IndFixList(I)-1,9),dM6(IndFixList(I)+1,9)]);
    dM6(IndFixList(I),12) = mean([dM6(IndFixList(I)-1,12),dM6(IndFixList(I)+1,12)]);
end
[sum(isnan(dM6(:,8))),sum(isnan(dM6(:,9))),sum(isnan(dM6(:,11))),sum(isnan(dM6(:,12)))] 



'done'
%% Plot


figure()
hold on
plot(dM2(:,2),dM2(:,3),'ok')
plot(dM2(:,5),dM2(:,6),'or')
plot(dM2(:,8),dM2(:,9),'ok')
plot(dM2(:,11),dM2(:,12),'or')
title('Mouse 2')

figure()
hold on
plot(dM3(:,2),dM3(:,3),'ok')
%plot(dM3(:,5),dM3(:,6),'or')
plot(dM3(:,8),dM3(:,9),'ok')
plot(dM3(:,11),dM3(:,12),'or')
title('Mouse 3')

figure()
hold on
plot(dM4(:,2),dM4(:,3),'ok')
%plot(dM4(:,5),dM4(:,6),'or')
plot(dM4(:,8),dM4(:,9),'ok')
plot(dM4(:,11),dM4(:,12),'or')
title('Mouse 4')

figure()
hold on
plot(dM6(:,5),dM6(:,6),'or')
plot(dM6(:,2),dM6(:,3),'ok')
plot(dM6(:,8),dM6(:,9),'ok')
%plot(dM6(:,11),dM6(:,12),'or')
title('Mouse 6')

%% Calculating angles with atan2 from mouse left whisker


%P = atan2(Y,X)
%2 = x, 3 = y

W1M2=(180/pi)*(atan2(dM2(:,3)-dM2(:,6),abs(dM2(:,5)-dM2(:,2))));  % convert rad to deg 180/pi
W2M2=(180/pi)*(atan2(dM2(:,9)-dM2(:,12),(dM2(:,11)-dM2(:,8))));
FixW2M2=(180/pi)*(atan2(dM2(1,9)-dM2(:,12),(dM2(:,11)-dM2(1,8))));
W1M2=W1M2-min(W1M2);
t1M2=linspace(0,numel(W1M2)/fs,numel(W1M2));
t2M2=linspace(0,numel(W2M2)/fs,numel(W2M2));

W1M3=(180/pi)*(atan2(dM3(:,3)-dM3(:,6),abs(dM3(:,5)-dM3(:,2))));  % convert rad to deg 180/pi
W2M3=(180/pi)*(atan2(dM3(:,9)-dM3(:,12),(dM3(:,11)-dM3(:,8))));
FixW2M3=(180/pi)*(atan2(dM3(1,9)-dM3(:,12),(dM3(:,11)-dM3(1,8))));
W1M3=W1M3-min(W1M3);
t1M3=linspace(0,numel(W1M3)/fs,numel(W1M3));
t2M3=linspace(0,numel(W2M3)/fs,numel(W2M3));

W1M4=(180/pi)*(atan2(dM4(:,3)-dM4(:,6),abs(dM4(:,5)-dM4(:,2))));  % convert rad to deg 180/pi
W2M4=(180/pi)*(atan2(dM4(:,9)-dM4(:,12),(dM4(:,11)-dM4(:,8))));
FixW2M4=(180/pi)*(atan2(dM4(1,9)-dM4(:,12),(dM4(:,11)-dM4(1,8))));
W1M4=W1M4-min(W1M4);
t1M4=linspace(0,numel(W1M4)/fs,numel(W1M4));
t2M4=linspace(0,numel(W2M4)/fs,numel(W2M4));


W1M6=(180/pi)*(atan2(dM6(:,3)-dM6(:,6),abs(dM6(:,5)-dM6(:,2))));  % convert rad to deg 180/pi
W2M6=(180/pi)*(atan2(dM6(:,9)-dM6(:,12),(dM6(:,11)-dM6(:,8))));
FixW2M6=(180/pi)*(atan2(dM6(1,9)-dM6(:,12),(dM6(:,11)-dM6(1,8))));
W1M6=W1M6-min(W1M6);
t1M6=linspace(0,numel(W1M6)/fs,numel(W1M6));
t2M6=linspace(0,numel(W2M6)/fs,numel(W2M6));

'done Calc angel'
%% plot angels

figure();
subplot(2,1,1)
plot(t1M2,smooth(W1M2,4)-mean(W1M2),'k')  %-min to ground
title('Mouse 2 Right whisker')
ylim([-pi*(180/pi)*(0.5) pi*(180/pi)*0.5])
ylabel('Degrees [Deg]')
xlim([0 60])
xlabel('Time [sec]')
subplot(2,1,2)
plot(t2M2,smooth(W2M2,4),'k')
title('Mouse 2 Left whisker')
xlim([0 60])
xlabel('Time [sec]')
ylim([-pi*(180/pi)*(0.5) pi*(180/pi)*0.5])
ylabel('Degrees [Deg]')

figure();
subplot(2,1,1)
plot(t1M3,smooth(W1M3,4)-mean(W1M3),'k')
title('Mouse 3 Right whisker')
ylim([-pi*(180/pi)*(0.5) pi*(180/pi)*0.5])
ylabel('Degrees [Deg]')
xlim([0 60])
xlabel('Time [sec]')
subplot(2,1,2)
plot(t2M3,smooth(W2M3,4),'k')
title('Mouse 3 Left whisker')
xlim([0 60])
xlabel('Time [sec]')
ylim([-pi*(180/pi)*(0.5) pi*(180/pi)*0.5])
ylabel('Degrees [Deg]')

figure();
subplot(2,1,1)
plot(t1M4,smooth(W1M4,4)-mean(W1M4),'k')
title('Mouse 4 Right whisker')
ylim([-pi*(180/pi)*(0.5) pi*(180/pi)*0.5])
ylabel('Degrees [Deg]')
xlim([0 60])
xlabel('Time [sec]')
subplot(2,1,2)
plot(t2M4,smooth(W2M4,4),'k')
title('Mouse 4 Left whisker')
xlim([0 60])
xlabel('Time [sec]')
ylim([-pi*(180/pi)*(0.5) pi*(180/pi)*0.5])
ylabel('Degrees [Deg]')


figure();
subplot(2,1,1)
plot(t1M6,W1M6,'k')
title('Mouse 6 Right whisker')
ylim([-pi*(180/pi)*(0.5) pi*(180/pi)*0.5])
ylabel('Degrees [Deg]')
xlim([0 60])
xlabel('Time [sec]')
subplot(2,1,2)
plot(t2M6,smooth(W2M6,4),'k')
title('Mouse 6 Left whisker')
xlim([0 60])
xlabel('Time [sec]')
ylim([-pi*(180/pi)*(0.5) pi*(180/pi)*0.5])
ylabel('Degrees [Deg]')



%% Creates summary plots of whisker left data 

% codes assumes that behaviour is in workspace as eg m2data and whisker
% angles in the for om of W2M2 (Whiskerside 2 mouse2)

Mouse = 2


grayColor = [.7 .7 .7]; scarlet=[0.8706 0.1490 0.2431];Skybule=[0.502 1 1]; YellowDark = [0.9294,0.6941,0.1255];
    
if Mouse > 1     % sets current data to data of mouse of choice

    onset=max(sData.daqdata.wheelDiode);

    X_CM = 0:1.5:156.5;

    if Mouse == 2
        Wdata=-W2M2;
        sData=m2data;
        VideoTime = 13 + (11/60);
    elseif Mouse ==3
        Wdata=-W2M3;
        sData=m3data;
        VideoTime = 15.3;
    elseif Mouse ==4
        Wdata=-W2M4;
        sData=m4data;
        VideoTime = 11 +(5/60);
    elseif Mouse ==5
        Wdata=-W2M5;
        sData=m5data;
        VideoTime = 16 +(35/60);
    elseif Mouse ==6
        
        Wdata=-W2M6;
        if Right ==1
            Wdata=W1M6; %W2
            Wdata=Wdata-100;
            Wdata(Wdata< -40) =  nanmean(Wdata);
            Wdata(isnan(Wdata)) = nanmean(Wdata);

            Wdata=-Wdata;
        end
        sData=m6data;
        VideoTime = 14 +(45/60);

    end
    
    data=Wdata;
    
    Frame = (sData.daqdata.frameIndex);
    
    Spd = sData.behavior.runSpeedDs;   % removes outlier points if pressent
    Spd(Spd < - 100) = 0;
    Spd(Spd > 100) = mean(Spd);
    

    ref=sData.daqdata.frameOnsetReferenceFrame;
    if Mouse ~= 6
        Pos = sData.behavior.wheelPosDsBinned; 
         if min(Pos)== 0
            Pos=Pos+1;
        end
        Lap = sData.behavior.wheelLapDsBinned;
        
        lick=sData.daqdata.lickSignal;
        lick=lick(Frame);
    else % mouse6
         Pos = sData.vr.positionBinnedDs;
         Lap = sData.vr.lapDs;
         lick=sData.daqdata.lickSignal;
         lick=lick(Frame);
         X_CM = 0:1:149;

    end
   
    tsm  = str2num(sData.sessionInfo.sessionStartTime(4:5)); % time start min
    tss  = str2num(sData.sessionInfo.sessionStartTime(7:8)); % time start sec
    tem = str2num(sData.sessionInfo.sessionStopTime(4:5));   % time end min
    tes = str2num(sData.sessionInfo.sessionStopTime(7:8));   % time end sec
    if Mouse == 4;
        tem =tem + 60;
    end
    Time = tem-tsm + (tes - tss)/60;
    Differ=abs(numel(Wdata)-numel(Spd));
end

    Wdata=Wdata(217:end);                        % removes the 7 second time delay between behaviour recoding and whisker recording   
    xq=0:(numel(Wdata)/numel(Spd)):numel(Wdata);
    Wdatat=interp1(1:numel(Wdata),Wdata,xq);
    Wdatat=Wdatat(2:end);
    Wdata=Wdatat;
    


%% whisking stats
if isnan(Wdata(1)) == 1
    Wdata(1) =0;
end

y = hilbert(Wdata);
sigphase = atan2(imag(y),real(y)); % Phase
hilb = (abs(y));                   % Hilbert 
if Mouse == 6 
    Amp = atan2(imag(y),real(y)); 
else
    Amp =(envelope(abs(y),300));       % Amplitude  
end

%% trial binner

LapIndex =diff(Lap);
LapIndex =find(LapIndex ~= 0);


% Zero-phase implementation - delay compensation

d1 = designfilt("lowpassiir",FilterOrder=12, ...
    HalfPowerFrequency=0.15,DesignMethod="butter");
ZeroPhase = filtfilt(d1,Wdata);
SetOff = ZeroPhase;


LickTrials ={};
SpeedTrials={};
PosTrials  ={};
WhiskerData={};
PhaseData  ={};
HilbData   ={};
AmpData    ={};
SFData     ={};
for i =1:numel(LapIndex)-1
    if i == 1
        LickTrials{i}  = lick(1:LapIndex(i)-1);
        SpeedTrials{i} = Spd(1:LapIndex(i)-1);
        PosTrials{i}   = Pos(1:LapIndex(i)-2);
        WhiskerData{i} = Wdata(1:LapIndex(i)-1);
        PhaseData{i}  =sigphase(1:LapIndex(i)-1);
        HilbData{i}   =hilb(1:LapIndex(i)-1);
        AmpData{i}    =Amp(1:LapIndex(i)-1);
        SFData{i}     =SetOff(1:LapIndex(i)-1);
    elseif i ==numel(LapIndex)-1
        LickTrials{i}  = lick(LapIndex(i):end);
        SpeedTrials{i} = Spd(LapIndex(i):end);
        PosTrials{i}   = Pos(LapIndex(i):end);
        WhiskerData{i} = Wdata(LapIndex(i):end);
        PhaseData{i}  =sigphase(LapIndex(i):end);
        HilbData{i}   =hilb(LapIndex(i):end);
        AmpData{i}    =Amp(LapIndex(i):end);
        SFData{i}     =SetOff(LapIndex(i):end);
    else
        LickTrials{i}  = lick(LapIndex(i-1)-1:LapIndex(i)-1);
        SpeedTrials{i} = Spd(LapIndex(i-1)-1:LapIndex(i)-1);
        PosTrials{i}   = Pos(LapIndex(i-1):LapIndex(i)-1);
        WhiskerData{i} = Wdata(LapIndex(i-1):LapIndex(i)-1);
        PhaseData{i}  =sigphase(LapIndex(i-1):LapIndex(i)-1);
        HilbData{i}   =hilb(LapIndex(i-1):LapIndex(i)-1);
        AmpData{i}    =Amp(LapIndex(i-1):LapIndex(i)-1);
        SFData{i}     =SetOff(LapIndex(i-1):LapIndex(i)-1);
    end

end


LickT=zeros(numel(LapIndex)-3,max(Pos)); 
SpeedT=zeros(numel(LapIndex)-3,max(Pos)); 
WhiskT=zeros(numel(LapIndex)-3,max(Pos));

PhaseT=zeros(numel(LapIndex)-3,max(Pos));
HilbT=zeros(numel(LapIndex)-3,max(Pos));
AmpT=zeros(numel(LapIndex)-3,max(Pos));

SetOffT=zeros(numel(LapIndex)-3,max(Pos));

for i = 2:numel(LapIndex)-2  % skips first trial
    L = LickTrials{i};
    S = SpeedTrials{i};
    P = PosTrials{i};
    W = WhiskerData{i};

    H = HilbData{i};
    A = AmpData{i};
    Ph= PhaseData{i};
    SF= SFData{i};

    for j= 1:max(Pos)
        PosIndex=find(P == j);

        LickT(i-1,j)  = nanmean(L(PosIndex));
        SpeedT(i-1,j) = nanmean(S(PosIndex));
        WhiskT(i-1,j) = nanmean(W(PosIndex));

        PhaseT(i-1,j)  = nanmean(Ph(PosIndex));
        AmpT(i-1,j) = nanmean(A(PosIndex));
        HilbT(i-1,j) = nanmean(H(PosIndex));
        SetOffT(i-1,j) = nanmean(SF(PosIndex));

    end


end

IndList=[33*1.5 73*1.5];

if Mouse ==6
    IndList=[50 110];
end

%% plots

X_time = linspace(0,Time,numel(Wdata(1:numel(Spd))));

Wn=78;
if Wn > numel(WhiskT(:,1))
    Wn = numel(WhiskT(:,1));
end

figure()
set(gcf,'color','w');
box off
subplot(4,1,1)
hold on
for i = 1:Wn
    plot(X_CM,WhiskT(i,:), color=grayColor)
end
plot(X_CM,nanmean(WhiskT(1:Wn,:)), 'k')
xline(IndList(1),'--',color=scarlet)
xline(IndList(2),'--',color=scarlet)
xlim([0 158])
ylabel({'Whisking', '(Deg)'})
%xlabel(['Distance [Cm]'])
%title('Whisking')

subplot(4,1,2)
hold on
for i = 1:Wn
    plot(X_CM,SetOffT(i,:), color=grayColor)
end
plot(X_CM,nanmean(SetOffT(1:Wn,:)), 'k')
xline(IndList(1),'--',color=scarlet)
xline(IndList(2),'--',color=scarlet)
xlim([0 158])
ylabel({'Offset' ,'(Deg)'})
%xlabel(['Distance [Cm]'])
%title('Offset')

subplot(4,1,3)
hold on
for i = 1:Wn
    plot(X_CM,SpeedT(i,:), color=grayColor)
end
plot(X_CM,nanmean(SpeedT(1:Wn,:)), 'k')
xline(IndList(1),'--',color=scarlet)
xline(IndList(2),'--',color=scarlet)
xlim([0 158])
ylabel({'Running speed', '(cm/s)'})

%title('Speed')

subplot(4,1,4)
Figpos = get(gcf, 'Position');
plot([0 50 110 157],[0 0 0 0],'|k')
set(gca, 'ytick',[])
ylim([0 0.1])
xlim([0 157.1])
set(gca,'color',[0.8 0.8 0.8])
xticks([0 50 110 157])
xlabel(['Distance (cm)'])



%% Summarypolots, Different xlims for nice segment for different mouse
motorPos1=find(Pos == 34);
motorPos2=find(Pos == 74);

if Mouse == 6    % binned differently then the others, motor position adjusted acordingly
    motorPos1=find(Pos == 50);
    motorPos2=find(Pos == 110);
end

t = linspace(0,numel(W2M3)/30,numel(Wdata)); 


if  Mouse == 2
    mouse1408=imread("C:\Pictures\Mouse 1412 Real.png");
elseif  Mouse == 3
    mouse1408=imread("C:\Pictures\Mouse 1414 Real.png");
elseif  Mouse == 4
    mouse1408=imread("C:\Pictures\5215V2.png");
else 
    mouse1408=imread("C:\mouse1408_3.png");
end

figure()

set(gcf,'color','w');
subplot(3,4,[1 5 ])
image([50 50], [40 40],mouse1408);
set(gca,'xtick',[],'ytick',[])

subplot(3,4,[2 3 4])
plot(Wdata,'k')
hold on
plot(ZeroPhase,'color',YellowDark,'linewidth',1.5)
%xlim([250/30 1052/30])
xlim([4230 4200+30*30])
%title('Whisking')
ylabel('Whisking (Deg)')
xline(motorPos1(20:40),'--',color=scarlet)
xline(motorPos2(20:40),'--',color=scarlet)
set(gca,'box', 'off','xtick',[])

subplot(3,4,[6 7 8])
plot(Pos*1.5,'k')
if Mouse == 6
    plot(Pos,'k')
end
%xlim([250 1052])
%xlim([4230 4200+30*30])
xlim([1713 2557])
xline(motorPos1(20:40),'--',color=scarlet)
xline(motorPos2(20:40),'--',color=scarlet)
%title('Position')
ylabel('Position (cm)')
set(gca,'box', 'off','xtick',[])

subplot(3,4,[10 11 12])
plot(Spd,'k')
xline(motorPos1(20:40),'--',color=scarlet)
xline(motorPos2(20:40),'--',color=scarlet)
%xlim([250 1052])
%xlim([4230 4200+30*30])
xlim([1713 2557])
%title('Velocity')
ylabel('Running speed (cm/s)')
xlabel('Time (s)')
set(gca,'box', 'off','xtick',[])

