MyParticipants = ls('../DATA/*_Y*');
for i = 1:size(MyParticipants,1)
    
[Num,Text,All] = xlsread(sprintf('../DATA/%s/%s',MyParticipants(i,:),ls(sprintf('../DATA/%s/ROOMQ*.csv',MyParticipants(i,:)))));
[NumBSC,TextBSC,AllBSC] = xlsread(sprintf('../DATA/%s/%s',MyParticipants(i,:),ls(sprintf('../DATA/%s/BSCQ*.csv',MyParticipants(i,:)))));
TextBSC = TextBSC(2:end,:)

P1 = [2 2 3 2 2 1 2 2 1 3 1 3 1 3 3 1 1 1 1 1 3 3 3 2 2 2]; %1PP synch environement 1 = Loft/living room, 2 = chnaging room, 3 = cabin
A1 = [1 3 1 1 3 3 3 3 3 2 3 1 2 1 1 3 2 2 2 3 2 2 2 1 1 3]; % 1PP asynch environement
A3 = [3 1 2 3 1 2 1 1 2 1 2 2 3 2 2 2 3 3 3 2 1 1 1 3 3 1]; %3PP asynch environement

Text_1 = contains(Text(:,1),'Q1]');nnz(Text_1)
TextQ1 = Text(Text_1,:);
Text_2 = contains(Text(:,1),'Q2]');nnz(Text_2)
TextQ2 = Text(Text_2,:);
 NumQ1 = Num(find(Text_1)-1,:);
NumQ2 = Num(find(Text_2)-1,:);
%to test resuzlts only of the first 10 trials
%NumQ1 = NumQ1(1:30,:);TextQ1 = TextQ1(1:30,:);

IndexEnv1 = contains(TextQ1(:,1),'Env1]');nnz(IndexEnv1)
IndexEnv2 = contains(TextQ1(:,1),'Env2]');nnz(IndexEnv2)
IndexEnv3 = contains(TextQ1(:,1),'Env4]');nnz(IndexEnv3)
TextQ1Env1 = TextQ1(IndexEnv1,1);
TextQ1Env2 = TextQ1(IndexEnv2,1);
TextQ1Env3 = TextQ1(IndexEnv3,1);
NumQ1Env1 = NumQ1(IndexEnv1,:);
NumQ1Env2 = NumQ1(find(IndexEnv2),:);
NumQ1Env3 = NumQ1(find(IndexEnv3),:);

if i == 4 || i == 12 ||  i ==22 || i ==26 %Mock Y04, Y12,Y23,Y30 the first 4 participants. The excel file obtained was slightly different from the rest. (we chaneg the way to obtain the data but notjing in the experiment)
NumQ1Env1(:,1) = NumQ1Env1(:,[1]);
NumQ1Env2(:,1) = NumQ1Env2(:,[1]);
NumQ1Env3(:,1) = NumQ1Env3(:,[1]);
%NumQ2(:,1) = NumQ2(:,2);
end
if i == 4 ||  i ==22
TimeQ1Env1 = TextQ1(IndexEnv1,3);
TimeQ2Env1 = TextQ2(IndexEnv1,3);
TimeQ1Env2 = TextQ1(IndexEnv2,3);
TimeQ2Env2 = TextQ2(IndexEnv2,3);
TimeQ1Env3 = TextQ1(IndexEnv3,3);
TimeQ2Env3 = TextQ2(IndexEnv3,3);
for j = 1:length(TimeQ1Env1)
TimeQ1Env1_(j) = str2num(TimeQ1Env1{j}(6:end-2))
TimeQ2Env1_(j) = str2num(TimeQ2Env1{j}(6:end-2))
TimeQ1Env2_(j) = str2num(TimeQ1Env2{j}(6:end-2))
TimeQ2Env2_(j) = str2num(TimeQ2Env2{j}(6:end-2))
TimeQ1Env3_(j) = str2num(TimeQ1Env3{j}(6:end-2))
TimeQ2Env3_(j) = str2num(TimeQ2Env3{j}(6:end-2))
end
else
TimeQ1Env1 = TextQ1(IndexEnv1,1);
TimeQ2Env1 = TextQ2(IndexEnv1,1);
TimeQ1Env2 = TextQ1(IndexEnv2,1);
TimeQ2Env2 = TextQ2(IndexEnv2,1);
TimeQ1Env3 = TextQ1(IndexEnv3,1);
TimeQ2Env3 = TextQ2(IndexEnv3,1);

TimeQ1Env1_ =zeros(45,1);TimeQ1Env2_ =zeros(45,1) ;TimeQ1Env3_ =zeros(45,1);
TimeQ2Env1_ =zeros(44,1);TimeQ2Env2_ =zeros(44,1) ;TimeQ2Env3_ =zeros(44,1);
c = 1;
for j = 1:length(TimeQ1Env1)
        if size(TimeQ1Env1{j},2) <65 %to avoid catch trial
TimeQ1Env1_(c,:) = str2num(TimeQ1Env1{j}(50:end-2))
TimeQ2Env1_(c,:) = str2num(TimeQ2Env1{j}(50:end-2))
        elseif size(TimeQ1Env1{j},2) <68
TimeQ1Env1_(c,:) = str2num(TimeQ1Env1{j}(54:end-2))
TimeQ2Env1_(c,:) = str2num(TimeQ2Env1{j}(54:end-2))
        else
TimeQ1Env1_(c,:) = str2num(TimeQ1Env1{j}(58:end-2))
TimeQ2Env1_(c,:) = str2num(TimeQ2Env1{j}(58:end-2))
        end
c = c+1;
        end
c = 1;
for j = 1:length(TimeQ1Env2)
        if size(TimeQ1Env2{j},2) <65 %to avoid catch trial
TimeQ1Env2_(c,:) = str2num(TimeQ1Env2{j}(50:end-2))
TimeQ2Env2_(c,:) = str2num(TimeQ2Env2{j}(50:end-2))
        elseif size(TimeQ1Env2{j},2) <68
TimeQ1Env2_(c,:) = str2num(TimeQ1Env2{j}(54:end-2))
TimeQ2Env2_(c,:) = str2num(TimeQ2Env2{j}(54:end-2))
        else
TimeQ1Env2_(c,:) = str2num(TimeQ1Env2{j}(58:end-2))
TimeQ2Env2_(c,:) = str2num(TimeQ2Env2{j}(58:end-2))
        end
c = c+1;
end
        

c = 1;
for j = 1:length(TimeQ1Env3)
        if size(TimeQ1Env3{j},2) <65 %to avoid catch trial
TimeQ1Env3_(c,:) = str2num(TimeQ1Env3{j}(50:end-2))
TimeQ2Env3_(c,:) = str2num(TimeQ2Env3{j}(50:end-2))
        elseif size(TimeQ1Env3{j},2) <68
TimeQ1Env3_(c,:) = str2num(TimeQ1Env3{j}(54:end-2))
TimeQ2Env3_(c,:) = str2num(TimeQ2Env3{j}(54:end-2))
        else
TimeQ1Env3_(c,:) = str2num(TimeQ1Env3{j}(58:end-2))
TimeQ2Env3_(c,:) = str2num(TimeQ2Env3{j}(58:end-2))
        end
c = c+1;
        end
end
% compute the time duration of the routine for the first question
%we just do the difference between the time ate which Q2 start and the time
%at which Q1 started
for j = 1:length(TimeQ1Env1)-1
    DurationQ1Env1_(j) = (TimeQ2Env1_(j)-TimeQ1Env1_(j))/1000
    DurationQ1Env2_(j) = (TimeQ2Env2_(j)-TimeQ1Env2_(j))/1000
    DurationQ1Env3_(j) = (TimeQ2Env3_(j)-TimeQ1Env3_(j))/1000
end

%for the second question, we do the difference between the time at which Q1
%from the NEXT routine start minus the time at which (q2 started + 15000 ms
%which corresponds to the new environnement observation and the black
%cross)
for j = 1:length(TimeQ1Env1)-1
  DurationQ2Env1_(j) = (TimeQ1Env2_(j)-(TimeQ2Env1_(j)+15000))/1000
  DurationQ2Env2_(j) = (TimeQ1Env3_(j)-(TimeQ2Env2_(j)+15000))/1000
  DurationQ2Env3_(j) = (TimeQ1Env1_(j+1)-(TimeQ2Env3_(j)+15000))/1000
end

MeanDEnv1_Q1 = mean(DurationQ1Env1_);
MeanDEnv1_Q2 = mean(DurationQ2Env1_);
MeanDEnv2_Q1 = mean(DurationQ1Env2_);
MeanDEnv2_Q2 = mean(DurationQ2Env2_);
MeanDEnv3_Q1 = mean(DurationQ1Env3_);
MeanDEnv3_Q2 = mean(DurationQ2Env3_);
IndexF_Env1 = contains(TextQ1Env1,'-F');
TotFalse_Env1 = nnz(IndexF_Env1);
IndexT_Env1 = contains(TextQ1Env1,'-T');
TotTrue_Env1 = nnz(IndexT_Env1)
T_Env1 = find(IndexT_Env1);
F_Env1 = find(IndexF_Env1);

IndexF_Env2 = contains(TextQ1Env2,'-F');nnz(IndexF_Env2)
IndexT_Env2= contains(TextQ1Env2,'-T');nnz(IndexT_Env2)
TotTrue_Env2 = nnz(IndexT_Env2)
TotFalse_Env2 = nnz(IndexF_Env2);

T_Env2 = find(IndexT_Env2);
F_Env2 = find(IndexF_Env2);


IndexF_Env3 = contains(TextQ1Env3,'-F');nnz(IndexF_Env3)
IndexT_Env3= contains(TextQ1Env3,'-T');nnz(IndexT_Env3)
TotTrue_Env3 = nnz(IndexT_Env3)
TotFalse_Env3 = nnz(IndexF_Env3);

T_Env3 = find(IndexT_Env3);
F_Env3 = find(IndexF_Env3);

% nature of the trials
TrialsEnv1 = zeros(length(NumQ1Env1),1); TrialsEnv2 = zeros(length(NumQ1Env2),1); TrialsEnv3 = zeros(length(NumQ1Env3),1); 
TrialsEnv1(IndexT_Env1) = 1; TrialsEnv1(IndexF_Env1) = 2; % T = 1, F = 2, Catch = 0 TialsEnv1 = 1 --> no change in room , so if answer correct CR. Trials ENV1 = 2 --> change in room so if answer correct = hit
TrialsEnv2(IndexT_Env2) = 1; TrialsEnv2(IndexF_Env2) = 2; % T = 1, F = 2, Catch = 0
TrialsEnv3(IndexT_Env3) = 1; TrialsEnv3(IndexF_Env3) = 2; % T = 1, F = 2, Catch = 0

PerfEnv1 = zeros(length(NumQ1Env1),1);PerfEnv2 = zeros(length(NumQ1Env2),1);PerfEnv3 = zeros(length(NumQ1Env3),1);
% PerfbisEnv1 = [];PerfbisEnv2 = zeros(45,1);PerfbisEnv3 = zeros(45,1);

for f = 1:length(NumQ1Env1)
    if NumQ1Env1(f) == 0 && TrialsEnv1(f) == 1
    PerfbisEnv1{f} = {'cr'};
    PerfEnv1(f) = 1;
    elseif NumQ1Env1(f) == 1 && TrialsEnv1(f) == 2
    PerfbisEnv1{f} = {'hit'};
    PerfEnv1(f) = 1;
    elseif NumQ1Env1(f) == 1 && TrialsEnv1(f) == 1
    PerfbisEnv1{f} = {'fa'};
    PerfEnv1(f) = 0;
    elseif NumQ1Env1(f) == 0 && TrialsEnv1(f) == 2
    PerfbisEnv1{f} = {'miss'};
    PerfEnv1(f) = 0;
    elseif NumQ1Env1(f) == 1 && TrialsEnv1(f) == 0
    PerfbisEnv1{f} = {'hit'};
    PerfEnv1(f) = 1;
    
    elseif NumQ1Env1(f) == 0 && TrialsEnv1(f) == 0
    PerfbisEnv1{f} = {'miss'};
    PerfEnv1(f) = 1;
    end
   
end

for f = 1:length(NumQ1Env2)
    if NumQ1Env2(f) == 0 && TrialsEnv2(f) == 1
    PerfbisEnv2{f} = {'cr'};
    PerfEnv2(f) = 1;
    elseif NumQ1Env2(f) == 1 && TrialsEnv2(f) == 2
    PerfbisEnv2{f} = {'hit'};
    PerfEnv2(f) = 1;
    elseif NumQ1Env2(f) == 1 && TrialsEnv2(f) == 1
    PerfbisEnv2{f} = {'fa'};
    PerfEnv2(f) = 0;
    elseif NumQ1Env2(f) == 0 && TrialsEnv2(f) == 2
    PerfbisEnv2{f} = {'miss'};
    PerfEnv2(f) = 0;
    elseif NumQ1Env2(f) == 1 && TrialsEnv2(f) == 0
    PerfbisEnv2{f} = {'hit'};
    PerfEnv2(f) = 1;
    elseif NumQ1Env2(f) == 0 && TrialsEnv2(f) == 0
    PerfbisEnv2{f} = {'miss'};
    PerfEnv2(f) = 1;
    end
end

for f = 1:length(NumQ1Env3)
    if NumQ1Env3(f) == 0 && TrialsEnv3(f) == 1
    PerfbisEnv3{f} = {'cr'};
    PerfEnv3(f) = 1;
    elseif NumQ1Env3(f) == 1 && TrialsEnv3(f) == 2
    PerfbisEnv3{f} = {'hit'};
    PerfEnv3(f) = 1;
    elseif NumQ1Env3(f) == 1 && TrialsEnv3(f) == 1
    PerfbisEnv3{f} = {'fa'};
    PerfEnv3(f) = 0;
    elseif NumQ1Env3(f) == 0 && TrialsEnv3(f) == 2
    PerfbisEnv3{f} = {'miss'};
    PerfEnv3(f) = 0;
    elseif NumQ1Env3(f) == 1 && TrialsEnv3(f) == 0
    PerfbisEnv3{f} = {'hit'};
    PerfEnv3(f) = 1;
    elseif NumQ1Env3(f) == 0 && TrialsEnv3(f) == 0
    PerfbisEnv3{f} = {'miss'};
    PerfEnv3(f) = 1;
    end
end

% % performance

[correctrej_Env1,falsealarm_Env1,correct] =Mistake(NumQ1Env1,T_Env1);
[correctrej_Env2,falsealarm_Env2] =Mistake(NumQ1Env2,T_Env2);
[correctrej_Env3,falsealarm_Env3] =Mistake(NumQ1Env3,T_Env3);
[misscount_Env1,hitcount_Env1] =Mistake(NumQ1Env1,F_Env1);
[misscount_Env2,hitcount_Env2] =Mistake(NumQ1Env2,F_Env2);
[misscount_Env3,hitcount_Env3] =Mistake(NumQ1Env3,F_Env3);

Performance_ENV1 = [hitcount_Env1 misscount_Env1 falsealarm_Env1 correctrej_Env1 length(T_Env1) length(F_Env1)];Performance_ENV2 = [hitcount_Env2 misscount_Env2 falsealarm_Env2 correctrej_Env2 length(T_Env2) length(F_Env2)];Performance_ENV3 = [hitcount_Env3 misscount_Env3 falsealarm_Env3 correctrej_Env3 length(T_Env3) length(F_Env3)];


[hitcount_Env1 hitcount_Env2; misscount_Env1 misscount_Env2;falsealarm_Env1 falsealarm_Env2; correctrej_Env1 correctrej_Env2]
% xlswrite('Analysis_Mock_YoungHealthy_020221_FullSample_Sheet.xlsx',[{'ParticipantsName'},{'HitEnv1'},{'MissEnv1'},{'FalseAlarmEnv1'},{'CorrectRejectionEnv1'},{'N_TrueEnv1'},{'N_FalseEnv1'},{'HitEnv2'},{'MissEnv2'},{'FalseAlarmEnv2'},{'CorrectRejectionEnv2'},{'N_TrueEnv2'},{'N_FalseEnv2'},{'HitEnv3'},{'MissEnv3'},{'FalseAlarmEnv3'},{'CorrectRejectionEnv3'},{'N_TrueEnv3'},{'N_FalseEnv3'}],'Performance_Env','A1:S1')
% xlswrite('Analysis_Mock_YoungHealthy_020221_FullSample_Sheet.xlsx',{MyParticipants(i,:)},'Performance_Env',sprintf('A%s',num2str(i+1)))
% xlswrite('Analysis_Mock_YoungHealthy_020221_FullSample_Sheet.xlsx',[hitcount_Env1 misscount_Env1 falsealarm_Env1 correctrej_Env1 length(T_Env1) length(F_Env1) hitcount_Env2 misscount_Env2 falsealarm_Env2 correctrej_Env2 length(T_Env2) length(F_Env2) hitcount_Env3 misscount_Env3 falsealarm_Env3 correctrej_Env3 length(T_Env3) length(F_Env3)],'Performance_Env',sprintf('B%s:S%s',num2str(i+1),num2str(i+1)))
% 
% xlswrite('Analysis_Mock_YoungHealthy_020221_FullSample_Sheet.xlsx',[{'ParticipantsName'},{'TimeQ1Env1'},{'TimeQ2Env1'},{'TimeQ1Env2'},{'TimeQ2Env2'},{'TimeQ1Env3'},{'TimeQ2Env3'}],'Time_Env','A1:G1')
% xlswrite('Analysis_Mock_YoungHealthy_020221_FullSample_Sheet.xlsx',{MyParticipants(i,:)},'Time_Env',sprintf('A%s',num2str(i+1)))
% xlswrite('Analysis_Mock_YoungHealthy_020221_FullSample_Sheet.xlsx',[MeanDEnv1_Q1 MeanDEnv1_Q2 MeanDEnv2_Q1 MeanDEnv2_Q2 MeanDEnv3_Q1 MeanDEnv3_Q2],'Time_Env',sprintf('B%s:G%s',num2str(i+1),num2str(i+1)))
% 
% 
% xlswrite('Analysis_Mock_YoungHealthy_020221_FullSample_Sheet.xlsx',[{'ParticipantsName'},{'Hit_Synch'},{'Miss_Synch'},{'FalseAlarm_Synch'},{'CorrectRejection_Synch'},{'N_True_Synch'},{'N_False_Synch'},{'Hit_Asynch'},{'Miss_Asynch'},{'FalseAlarm_Asynch'},{'CorrectRejection_Asynch'},{'N_True_Asynch'},{'N_False_Asynch'},{'Hit_Asynch3PP'},{'Miss_Asynch3PP'},{'FalseAlarm_Asynch3PP'},{'CorrectRejection_Asynch3PP'},{'N_True_Asynch3PP'},{'N_False_Asynch3PP'}],'Performance_Cond','A1:S1')
% xlswrite('Analysis_Mock_YoungHealthy_020221_FullSample_Sheet.xlsx',{MyParticipants(i,:)},'Performance_Cond',sprintf('A%s',num2str(i+1)))
% xlswrite('Analysis_Mock_YoungHealthy_020221_FullSample_Sheet.xlsx',[{'ParticipantsName'},{'TimeQ1Env1'},{'TimeQ2Env1'},{'TimeQ1Env2'},{'TimeQ2Env2'},{'TimeQ1Env3'},{'TimeQ2Env3'}],'Time_Cond','A1:G1')
% xlswrite('Analysis_Mock_YoungHealthy_020221_FullSample_Sheet.xlsx',{MyParticipants(i,:)},'Time_Cond',sprintf('A%s',num2str(i+1)))
 
if P1(i) == 1 && A1(i) == 2 && A3(i) == 3
% xlswrite('Analysis_Mock_YoungHealthy_020221_FullSample_Sheet.xlsx',[Performance_ENV1 Performance_ENV2 Performance_ENV3],'Performance_Cond',sprintf('B%s:S%s',num2str(i+1),num2str(i+1)))
% xlswrite('Analysis_Mock_YoungHealthy_020221_FullSample_Sheet.xlsx',[MeanDEnv1_Q1 MeanDEnv1_Q2 MeanDEnv2_Q1 MeanDEnv2_Q2 MeanDEnv3_Q1 MeanDEnv3_Q2],'Time_Cond',sprintf('B%s:G%s',num2str(i+1),num2str(i+1)))

elseif  P1(i) == 2 && A1(i) == 3 && A3(i) == 1
% xlswrite('Analysis_Mock_YoungHealthy_020221_FullSample_Sheet.xlsx',[Performance_ENV2 Performance_ENV3 Performance_ENV1],'Performance_Cond',sprintf('B%s:S%s',num2str(i+1),num2str(i+1)))
% xlswrite('Analysis_Mock_YoungHealthy_020221_FullSample_Sheet.xlsx',[MeanDEnv2_Q1 MeanDEnv2_Q2 MeanDEnv3_Q1 MeanDEnv3_Q2 MeanDEnv1_Q1 MeanDEnv1_Q2],'Time_Cond',sprintf('B%s:G%s',num2str(i+1),num2str(i+1)))

elseif  P1(i) == 2 && A1(i) == 1 && A3(i) == 3
% xlswrite('Analysis_Mock_YoungHealthy_020221_FullSample_Sheet.xlsx',[Performance_ENV2 Performance_ENV1 Performance_ENV3],'Performance_Cond',sprintf('B%s:S%s',num2str(i+1),num2str(i+1)))
% xlswrite('Analysis_Mock_YoungHealthy_020221_FullSample_Sheet.xlsx',[MeanDEnv2_Q1 MeanDEnv2_Q2 MeanDEnv1_Q1 MeanDEnv1_Q2 MeanDEnv3_Q1 MeanDEnv3_Q2],'Time_Cond',sprintf('B%s:G%s',num2str(i+1),num2str(i+1)))

elseif P1(i) == 1 && A1(i) == 3 && A3(i) == 2
% xlswrite('Analysis_Mock_YoungHealthy_020221_FullSample_Sheet.xlsx',[Performance_ENV1 Performance_ENV3 Performance_ENV2],'Performance_Cond',sprintf('B%s:S%s',num2str(i+1),num2str(i+1)))
% xlswrite('Analysis_Mock_YoungHealthy_020221_FullSample_Sheet.xlsx',[MeanDEnv1_Q1 MeanDEnv1_Q2 MeanDEnv3_Q1 MeanDEnv3_Q2 MeanDEnv2_Q1 MeanDEnv2_Q2],'Time_Cond',sprintf('B%s:G%s',num2str(i+1),num2str(i+1)))

elseif P1(i) == 3 && A1(i) == 1 && A3(i) == 2
% xlswrite('Analysis_Mock_YoungHealthy_020221_FullSample_Sheet.xlsx',[Performance_ENV3 Performance_ENV1 Performance_ENV2],'Performance_Cond',sprintf('B%s:S%s',num2str(i+1),num2str(i+1)))
% xlswrite('Analysis_Mock_YoungHealthy_020221_FullSample_Sheet.xlsx',[MeanDEnv3_Q1 MeanDEnv3_Q2 MeanDEnv1_Q1 MeanDEnv1_Q2 MeanDEnv2_Q1 MeanDEnv2_Q2],'Time_Cond',sprintf('B%s:G%s',num2str(i+1),num2str(i+1)))

elseif P1(i) == 3 && A1(i) == 2 && A3(i) == 1
% xlswrite('Analysis_Mock_YoungHealthy_020221_FullSample_Sheet.xlsx',[Performance_ENV3 Performance_ENV2 Performance_ENV1],'Performance_Cond',sprintf('B%s:S%s',num2str(i+1),num2str(i+1)))
% xlswrite('Analysis_Mock_YoungHealthy_020221_FullSample_Sheet.xlsx',[MeanDEnv3_Q1 MeanDEnv3_Q2 MeanDEnv2_Q1 MeanDEnv2_Q2 MeanDEnv1_Q1 MeanDEnv1_Q2],'Time_Cond',sprintf('B%s:G%s',num2str(i+1),num2str(i+1)))

end

% Confidence
MeanConfidenceT_Env1 = nanmean(NumQ2(T_Env1,1));
MeanConfidenceT_Env2 = nanmean(NumQ2(T_Env2,1));
MeanConfidenceT_Env3 = nanmean(NumQ2(T_Env3,1));
MeanConfidenceF_Env1 = nanmean(NumQ2(F_Env1,1));
MeanConfidenceF_Env2 = nanmean(NumQ2(F_Env2,1));
MeanConfidenceF_Env3 = nanmean(NumQ2(F_Env3,1));
[MeanConfidenceT_Env1 MeanConfidenceT_Env2; MeanConfidenceF_Env1 MeanConfidenceF_Env2]
CONF_ENV1 = [MeanConfidenceT_Env1 MeanConfidenceF_Env1]; CONF_ENV2 =  [MeanConfidenceT_Env2 MeanConfidenceF_Env2];CONF_ENV3 =  [MeanConfidenceT_Env3 MeanConfidenceF_Env3]
% xlswrite('Analysis_Mock_YoungHealthy_020221_FullSample_Sheet.xlsx',[{'ParticipantsName'},{'Confidence_T_Env1'},{'Confidence_F_Env1'},{'Confidence_T_Env2'},{'Confidence_F_Env2'},{'Confidence_T_Env3'},{'Confidence_F_Env3'}],'Confidence','A1:G1')
% xlswrite('Analysis_Mock_YoungHealthy_020221_FullSample_Sheet.xlsx',{MyParticipants(i,:)},'Confidence',sprintf('A%s',num2str(i+1)))
% xlswrite('Analysis_Mock_YoungHealthy_020221_FullSample_Sheet.xlsx',[MeanConfidenceT_Env1 MeanConfidenceF_Env1 MeanConfidenceT_Env2 MeanConfidenceF_Env2 MeanConfidenceT_Env3 MeanConfidenceF_Env3],'Confidence',sprintf('B%s:G%s',num2str(i+1),num2str(i+1)))
% 
% xlswrite('Analysis_Mock_YoungHealthy_020221_FullSample_Sheet.xlsx',[{'ParticipantsName'},{'Confidence_T_Env1'},{'Confidence_F_Env1'},{'Confidence_T_Env2'},{'Confidence_F_Env2'},{'Confidence_T_Env3'},{'Confidence_F_Env3'}],'Confidence_Cond','A1:G1')
% xlswrite('Analysis_Mock_YoungHealthy_020221_FullSample_Sheet.xlsx',{MyParticipants(i,:)},'Confidence_Cond',sprintf('A%s',num2str(i+1)))
% if P1(i) == 1 && A1(i) == 2 && A3(i) == 3
% xlswrite('Analysis_Mock_YoungHealthy_020221_FullSample_Sheet.xlsx',[CONF_ENV1 CONF_ENV2 CONF_ENV3],'Confidence_Cond',sprintf('B%s:G%s',num2str(i+1),num2str(i+1)))
% elseif  P1(i) == 2 && A1(i) == 3 && A3(i) == 1
% xlswrite('Analysis_Mock_YoungHealthy_020221_FullSample_Sheet.xlsx',[CONF_ENV2 CONF_ENV3 CONF_ENV1],'Confidence_Cond',sprintf('B%s:G%s',num2str(i+1),num2str(i+1)))
% elseif  P1(i) == 2 && A1(i) == 1 && A3(i) == 3
% xlswrite('Analysis_Mock_YoungHealthy_020221_FullSample_Sheet.xlsx',[CONF_ENV2 CONF_ENV1 CONF_ENV3],'Confidence_Cond',sprintf('B%s:G%s',num2str(i+1),num2str(i+1)))
% elseif P1(i) == 1 && A1(i) == 3 && A3(i) == 2
% xlswrite('Analysis_Mock_YoungHealthy_020221_FullSample_Sheet.xlsx',[CONF_ENV1 CONF_ENV3 CONF_ENV2],'Confidence_Cond',sprintf('B%s:G%s',num2str(i+1),num2str(i+1)))
% elseif P1(i) == 3 && A1(i) == 1 && A3(i) == 2
% xlswrite('Analysis_Mock_YoungHealthy_020221_FullSample_Sheet.xlsx',[CONF_ENV3 CONF_ENV1 CONF_ENV2],'Confidence_Cond',sprintf('B%s:G%s',num2str(i+1),num2str(i+1)))
% elseif P1(i) == 3 && A1(i) == 2 && A3(i) == 1
% xlswrite('Analysis_Mock_YoungHealthy_020221_FullSample_Sheet.xlsx',[CONF_ENV3 CONF_ENV2 CONF_ENV1],'Confidence_Cond',sprintf('B%s:G%s',num2str(i+1),num2str(i+1)))
% end


%%BSCQ
TextBSC_1 = contains(TextBSC(:,1),':SYNCH_1PP');nnz(TextBSC_1)
TextBSC_Env1 = TextBSC(TextBSC_1,:);
TextBSC_2 = contains(TextBSC(:,1),':ASYNCH_1PP');nnz(TextBSC_2)
TextBSC_Env2 = TextBSC(TextBSC_2,:);
TextBSC_3 = contains(TextBSC(:,1),':ASYNCH_3PP');nnz(TextBSC_3)
TextBSC_Env3 = TextBSC(TextBSC_3,:);
NumBSC_1 = NumBSC(find(TextBSC_1),:)
NumBSC_2 = NumBSC(find(TextBSC_2),:)
NumBSC_3 = NumBSC(find(TextBSC_3),:)

Ownership_1 = NumBSC_1(contains(TextBSC_Env1(:,1),'OWNERSHIP]'));nnz(Ownership_1)
Agency_1 = NumBSC_1(contains(TextBSC_Env1(:,1),'AGENCY]'));nnz(Agency_1)
Control_1 = NumBSC_1(contains(TextBSC_Env1(:,1),'CONTROL]'));nnz(Control_1)
Control2_1 = NumBSC_1(contains(TextBSC_Env1(:,1),'CONTROL2]'));nnz(Control2_1)
Threat_1 = NumBSC_1(contains(TextBSC_Env1(:,1),'THREAT]'));nnz(Threat_1)

Ownership_2 = NumBSC_2(contains(TextBSC_Env2(:,1),'OWNERSHIP]'));nnz(Ownership_2)
Agency_2 = NumBSC_2(contains(TextBSC_Env2(:,1),'AGENCY]'));nnz(Agency_2)
Control_2 = NumBSC_2(contains(TextBSC_Env2(:,1),'CONTROL]'));nnz(Control_2)
Control2_2 = NumBSC_2(contains(TextBSC_Env2(:,1),'CONTROL]'));nnz(Control2_2)
Threat_2 = NumBSC_2(contains(TextBSC_Env2(:,1),'THREAT]'));nnz(Threat_2)


Ownership_3 = NumBSC_3(contains(TextBSC_Env3(:,1),'OWNERSHIP]'));nnz(Ownership_3)
Agency_3 = NumBSC_3(contains(TextBSC_Env3(:,1),'AGENCY]'));nnz(Agency_3)
Control_3 = NumBSC_3(contains(TextBSC_Env3(:,1),'CONTROL]'));nnz(Control_3)
Control2_3 = NumBSC_3(contains(TextBSC_Env3(:,1),'CONTROL]'));nnz(Control2_3)
Threat_3 = NumBSC_3(contains(TextBSC_Env3(:,1),'THREAT]'));nnz(Threat_3)

Ownership_Env1 = nanmean(Ownership_1); Agency_Env1 = nanmean(Agency_1); Control_Env1 = nanmean(Control_1);Control2_Env1 = nanmean(Control2_1); Threat_Env1 = nanmean(Threat_1);
Ownership_Env2 = nanmean(Ownership_2); Agency_Env2 = nanmean(Agency_2); Control_Env2 = nanmean(Control_2); Control2_Env2 = nanmean(Control2_2); Threat_Env2 = nanmean(Threat_2);
Ownership_Env3 = nanmean(Ownership_3); Agency_Env3 = nanmean(Agency_3); Control_Env3 = nanmean(Control_3); Control2_Env3 = nanmean(Control2_3); Threat_Env3 = nanmean(Threat_3);

% Here, for BSC, BSC... ENV1 corresponds to synch, Env2 correspond to asynch and env3 to
% asynch 3pp
BSC_SYNCH = [Ownership_Env1 Agency_Env1 Control_Env1 Control2_Env1 Threat_Env1]; BSC_ASYNCH = [Ownership_Env2 Agency_Env2 Control_Env2 Control2_Env2 Threat_Env2]; BSC_ASYNCH3PP = [Ownership_Env3 Agency_Env3 Control_Env3 Control2_Env3 Threat_Env3];


Conditions = [P1' A1' A3']
    if find(Conditions(i,:) ==1) ==1
        Env1_Cond = 'SYNCH1PP';
        BSC1 = BSC_SYNCH;
       
        PERFSYNCH(i,:) = PerfEnv1;
    elseif find(Conditions(i,:) ==1) ==2
        Env1_Cond = 'ASYNCH1PP';
       BSC1 = BSC_ASYNCH;
         PERFASYNCH(i,:) = PerfEnv1;
    elseif find(Conditions(i,:) ==1) ==3
        Env1_Cond = 'ASYNCH3PP';
       BSC1 = BSC_ASYNCH3PP;
         PERFASYNCH3PP(i,:) = PerfEnv1;
    end
    
    if find(Conditions(i,:) ==2) ==1
        Env2_Cond = 'SYNCH1PP';
        BSC2 = BSC_SYNCH
         PERFSYNCH(i,:) = PerfEnv2;
    elseif find(Conditions(i,:) ==2) ==2
        Env2_Cond = 'ASYNCH1PP';
    BSC2 = BSC_ASYNCH;
         PERFASYNCH(i,:) = PerfEnv2;
    elseif find(Conditions(i,:) ==2) ==3
        Env2_Cond = 'ASYNCH3PP';
        BSC2 = BSC_ASYNCH3PP;
         PERFASYNCH3PP(i,:) = PerfEnv2;
    end
    
    if find(Conditions(i,:) ==3) ==1
        Env3_Cond = 'SYNCH1PP';
       BSC3 = BSC_SYNCH;
         PERFSYNCH(i,:) = PerfEnv3;
    elseif find(Conditions(i,:) ==3) ==2
        Env3_Cond = 'ASYNCH1PP';
    BSC3 = BSC_ASYNCH;
          PERFASYNCH(i,:) = PerfEnv3;
    elseif find(Conditions(i,:) ==3) ==3
        Env3_Cond = 'ASYNCH3PP';
        BSC3 = BSC_ASYNCH3PP;
         PERFASYNCH3PP(i,:) = PerfEnv3;
    end



% xlswrite('Analysis_Mock_YoungHealthy_020221_FullSample_Sheet.xlsx',[{'ParticipantsName'},{'OWNERSHIP_1PPsynch'},{'AGENCY_1PPsynch'},{'CONTROL_1PPsynch'},{'CONTROL2_1PPsynch'},{'THREAT_1PPsynch'},{'OWNERSHIP_1PPasynch'},{'AGENCY_1PPasynch'},{'CONTROL_1PPasynch'},{'CONTROL2_1PPasynch'},{'THREAT_1PPasynch'},{'OWNERSHIP_3PPasynch'},{'AGENCY_3PPasynch'},{'CONTROL_3PPasynch'},{'CONTROL2_3PPasynch'},{'THREAT_3PPasynch'}],'BSC_Cond','A1:P1')
% xlswrite('Analysis_Mock_YoungHealthy_020221_FullSample_Sheet.xlsx',{MyParticipants(i,:)},'BSC_Cond',sprintf('A%s',num2str(i+1)))
% xlswrite('Analysis_Mock_YoungHealthy_020221_FullSample_Sheet.xlsx',[Ownership_Env1 Agency_Env1 Control_Env1 Control2_Env1 Threat_Env1 Ownership_Env2 Agency_Env2 Control_Env2 Control2_Env2 Threat_Env2 Ownership_Env3 Agency_Env3 Control_Env3 Control2_Env3 Threat_Env3],'BSC_Cond',sprintf('B%s:P%s',num2str(i+1),num2str(i+1)))
%for 10 trials
xlswrite('Analysis_Mock_YoungHealthy_020221_FullSample.xlsx',[{'ParticipantsName'},{'Hit'},{'Miss'},{'FalseAlarm'},{'CorrectRejection'},{'N_True'},{'N_False'},{'Environment'},{'Conditions'},{'Answer'},{'Performance'},{'Performance_HitCR'},{'Trials'},{'Ownership'},{'Agency'},{'Control1'},{'Control2'},{'Threat'},{'MeanQ1'},{'MeanQ2'}],'Summary','A1:T1')
xlswrite('Analysis_Mock_YoungHealthy_020221_FullSample.xlsx',[repmat(string(MyParticipants(i,:)),length(NumQ1Env1),1) repmat(hitcount_Env1,length(NumQ1Env1),1) repmat(misscount_Env1,length(NumQ1Env1),1) repmat(falsealarm_Env1,length(NumQ1Env1),1)  repmat(correctrej_Env1,length(NumQ1Env1),1) repmat(length(T_Env1),length(NumQ1Env1),1)  repmat(length(F_Env1),length(NumQ1Env1),1) repmat('ENV1',length(NumQ1Env1),1) repmat(Env1_Cond,length(NumQ1Env1),1) NumQ1Env1(:,1) PerfEnv1 string(PerfbisEnv1') TrialsEnv1 repmat(BSC1,length(NumQ1Env1),1) repmat(MeanDEnv1_Q1,length(NumQ1Env1),1) repmat(MeanDEnv1_Q2,length(NumQ1Env1),1) ],'Summary',sprintf('A%s:T%s',num2str(2 + (i-1)*3*length(NumQ1Env1)),num2str(2+(i-1)*3*length(NumQ1Env1)+(length(NumQ1Env1)-1))))
xlswrite('Analysis_Mock_YoungHealthy_020221_FullSample.xlsx',[repmat(string(MyParticipants(i,:)),length(NumQ1Env1),1) repmat(hitcount_Env2,length(NumQ1Env1),1) repmat(misscount_Env2,length(NumQ1Env1),1) repmat(falsealarm_Env2,length(NumQ1Env1),1)  repmat(correctrej_Env2,length(NumQ1Env1),1) repmat(length(T_Env2),length(NumQ1Env1),1)  repmat(length(F_Env2),length(NumQ1Env1),1) repmat('ENV2',length(NumQ1Env1),1) repmat(Env2_Cond,length(NumQ1Env1),1) NumQ1Env2(:,1) PerfEnv2 string(PerfbisEnv2') TrialsEnv2 repmat(BSC2,length(NumQ1Env1),1) repmat(MeanDEnv2_Q1,length(NumQ1Env2),1) repmat(MeanDEnv2_Q2,length(NumQ1Env2),1)],'Summary',sprintf('A%s:T%s',num2str(2+length(NumQ1Env1) + (i-1)*3*length(NumQ1Env1)),num2str(2+(i-1)*3*length(NumQ1Env1)+(2*length(NumQ1Env1)-1))))
xlswrite('Analysis_Mock_YoungHealthy_020221_FullSample.xlsx',[repmat(string(MyParticipants(i,:)),length(NumQ1Env1),1) repmat(hitcount_Env3,length(NumQ1Env1),1) repmat(misscount_Env3,length(NumQ1Env1),1) repmat(falsealarm_Env3,length(NumQ1Env1),1)  repmat(correctrej_Env3,length(NumQ1Env1),1) repmat(length(T_Env3),length(NumQ1Env1),1)  repmat(length(F_Env3),length(NumQ1Env1),1) repmat('ENV3',length(NumQ1Env1),1) repmat(Env3_Cond,length(NumQ1Env1),1) NumQ1Env3(:,1) PerfEnv3 string(PerfbisEnv3') TrialsEnv3 repmat(BSC3,length(NumQ1Env1),1) repmat(MeanDEnv3_Q1,length(NumQ1Env3),1) repmat(MeanDEnv3_Q2,length(NumQ1Env3),1)],'Summary',sprintf('A%s:T%s',num2str(2+(2*length(NumQ1Env1)) + (i-1)*3*length(NumQ1Env1)),num2str(2+(i-1)*3*length(NumQ1Env1)+(3*length(NumQ1Env1)-1))))

%for 45 trials
% xlswrite('Analysis_Mock_YoungHealthy_020221_FullSample_Sheet.xlsx',[{'ParticipantsName'},{'Hit'},{'Miss'},{'FalseAlarm'},{'CorrectRejection'},{'N_True'},{'N_False'},{'Environment'},{'Conditions'},{'Answer'},{'Performance'},{'Performance_HitCR'},{'Trials'},{'Ownership'},{'Agency'},{'Control1'},{'Control2'},{'Threat'}],'Summary','A1:R1')
% xlswrite('Analysis_Mock_YoungHealthy_020221_FullSample_Sheet.xlsx',[repmat(string(MyParticipants(i,:)),45,1) repmat(hitcount_Env1,45,1) repmat(misscount_Env1,45,1) repmat(falsealarm_Env1,45,1)  repmat(correctrej_Env1,45,1) repmat(length(T_Env1),45,1)  repmat(length(F_Env1),45,1) repmat('ENV1',45,1) repmat(Env1_Cond,45,1) NumQ1Env1(:,1) PerfEnv1 string(PerfbisEnv1') TrialsEnv1 repmat(BSC1,45,1) ],'Summary',sprintf('A%s:R%s',num2str(2 + (i-1)*135),num2str(2+(i-1)*135+44)))
% xlswrite('Analysis_Mock_YoungHealthy_020221_FullSample_Sheet.xlsx',[repmat(string(MyParticipants(i,:)),45,1) repmat(hitcount_Env2,45,1) repmat(misscount_Env2,45,1) repmat(falsealarm_Env2,45,1)  repmat(correctrej_Env2,45,1) repmat(length(T_Env2),45,1)  repmat(length(F_Env2),45,1) repmat('ENV2',45,1) repmat(Env2_Cond,45,1) NumQ1Env2(:,1) PerfEnv2 string(PerfbisEnv2') TrialsEnv2 repmat(BSC2,45,1) ],'Summary',sprintf('A%s:R%s',num2str(47 + (i-1)*135),num2str(2+(i-1)*135+89)))
% xlswrite('Analysis_Mock_YoungHealthy_020221_FullSample_Sheet.xlsx',[repmat(string(MyParticipants(i,:)),45,1) repmat(hitcount_Env3,45,1) repmat(misscount_Env3,45,1) repmat(falsealarm_Env3,45,1)  repmat(correctrej_Env3,45,1) repmat(length(T_Env3),45,1)  repmat(length(F_Env3),45,1) repmat('ENV3',45,1) repmat(Env3_Cond,45,1) NumQ1Env3(:,1) PerfEnv3 string(PerfbisEnv3') TrialsEnv3 repmat(BSC3,45,1) ],'Summary',sprintf('A%s:R%s',num2str(92 + (i-1)*135),num2str(2+(i-1)*135+134)))



% figure(1)
% subplot(5,2,i)
% plot(1:length(NumQ2(IndexEnv1)),NumQ2(IndexEnv1),'ob')
% hold on
sprintf('../DATA/%s/%s',MyParticipants(i,:),ls(sprintf('../DATA/%s/ROOMQ*.csv',MyParticipants(i,:))))
[F1_r1,F2_r1,F3_r1,F4_r1,F5_r1,F6_r1,F7_r1,F8_r1] = MistakeDetailed(sprintf('../DATA/%s/%s',MyParticipants(i,:),ls(sprintf('../DATA/%s/ROOMQ*.csv',MyParticipants(i,:)))),MyParticipants,i,'Env1') % check object difficulty for kitchen (env1)
[F1_r2,F2_r2,F3_r2,F4_r2,F5_r2,F6_r2,F7_r2,F8_r2] = MistakeDetailed(sprintf('../DATA/%s/%s',MyParticipants(i,:),ls(sprintf('../DATA/%s/ROOMQ*.csv',MyParticipants(i,:)))),MyParticipants,i,'Env2')% check object difficulty for labroom (env2)
[F1_r3,F2_r3,F3_r3,F4_r3,F5_r3,F6_r3,F7_r3,F8_r3] = MistakeDetailed(sprintf('../DATA/%s/%s',MyParticipants(i,:),ls(sprintf('../DATA/%s/ROOMQ*.csv',MyParticipants(i,:)))),MyParticipants,i,'Env4')% check object difficulty for cabin (env4)

Object_ENV1 = [F1_r1,F2_r1,F3_r1,F4_r1,F5_r1,F6_r1,F7_r1,F8_r1];
Object_ENV2 = [F1_r2,F2_r2,F3_r2,F4_r2,F5_r2,F6_r2,F7_r2,F8_r2];
Object_ENV3 = [F1_r3,F2_r3,F3_r3,F4_r3,F5_r3,F6_r3,F7_r3,F8_r3];
% xlswrite('Analysis_Mock_YoungHealthy_020221_FullSample_Sheet.xlsx',[{MyParticipants(i,:)}],'Difficulty',sprintf('A%s',num2str(i+1)))
% xlswrite('Analysis_Mock_YoungHealthy_020221_FullSample_Sheet.xlsx',[{'ParticipantsName'},{'F1_Env1'},{'F2_Env1'},{'F3_Env1'},{'F4_Env1'},{'F5_Env1'},{'F6_Env1'},{'F7_Env1'},{'F8_Env1'},{'F1_Env2'},{'F2_Env2'},{'F3_Env2'},{'F4_Env2'},{'F5_Env2'},{'F6_Env2'},{'F7_Env2'},{'F8_Env2'},{'F1_Env4'},{'F2_Env4'},{'F3_Env4'},{'F4_Env4'},{'F5_Env4'},{'F6_Env4'},{'F7_Env4'},{'F8_Env4'}],'Difficulty','A1:Z1')
% xlswrite('Analysis_Mock_YoungHealthy_020221_FullSample_Sheet.xlsx',[F1_r1 F2_r1 F3_r1 F4_r1 F5_r1 F6_r1 F7_r1 F8_r1 F1_r2 F2_r2 F3_r2 F4_r2 F5_r2 F6_r2 F7_r2 F8_r2 F1_r3,F2_r3,F3_r3,F4_r3,F5_r3,F6_r3,F7_r3,F8_r3],'Difficulty',sprintf('B%s:Z%s',num2str(i+1),num2str(i+1)))
% 


 Num = []; text = []; All = [];NumQ1 = []; NumQ2 = []; NumQ1Env1 = []; NumQ2Env1 = [];NumQ1Env2 = []; NumQ2Env2 = [];NumBSC_1 = []; NumBSC_2 = [];NumBSC_3 = []; TextBSC_1 = []; TextBSC_2 = []; TextBSC_3 = [];TimeQ1Env1_ = [];TimeQ2Env1_ = [];TimeQ1Env2_ = [];TimeQ2Env2_ = [];TimeQ1Env3_ = [];TimeQ2Env3_ = [];
 TimeQ1Env1 = []; TimeQ2Env1 = [];TimeQ1Env2 = []; TimeQ2Env2 = [];TimeQ1Env3 = []; TimeQ2Env3 = [];
end


%Mock Y30 and Mock Y20 are outliers, so we don't want to plot them
% Plot performance per trial of 5
c = 1;k = 1;
for i = 1:5:45
    for j = 1:size(MyParticipants,1)
        if j ~= 20  && j ~= 26 
        PerfSynch_5(j,c) = nnz(PERFSYNCH(j,i:i+4))/5;
        PerfAsynch_5(j,c) = nnz(PERFASYNCH(j,i:i+4))/5;
        PerfAsynch3pp_5(j,c) = nnz(PERFASYNCH3PP(j,i:i+4))/5;
       
        end
    end
    c = c+1; 
end
 PerfSynch_5(20,:) = [];  PerfAsynch_5(20,:) = []; PerfAsynch3pp_5(20,:) = [];

figure(1)
subplot(1,3,1)

for i = 1:size(MyParticipants,1)-2
    plot(1:9,PerfSynch_5(i,:),'-o')
    hold on
end
 plot(1:9,mean(PerfSynch_5,1),'-o','LineWidth',2,'color','k')
 title('Synch')
 subplot(1,3,2)

for i = 1:size(MyParticipants,1)-2
    plot(1:9,PerfAsynch_5(i,:),'-o')
    hold on
end
 plot(1:9,mean(PerfAsynch_5,1),'-o','LineWidth',2,'color','k')
  title('Asynch')
 subplot(1,3,3)
 
 for i = 1:size(MyParticipants,1)-2
    plot(1:9,PerfAsynch3pp_5(i,:),'-o')
    hold on
end
 plot(1:9,mean(PerfAsynch3pp_5,1),'-o','LineWidth',2,'color','k')
 title('Asynch3pp')
 
 
 
 c = 1; k = 1;
 for i = 1:45-4
    for j = 1:size(MyParticipants,1)
          if j ~= 20  && j ~= 26 
        PerfSynch_slide(j,c) = nnz(PERFSYNCH(j,i:i+4))/5;
        PerfAsynch_slide(j,c) = nnz(PERFASYNCH(j,i:i+4))/5;
        PerfAsynch3pp_slide(j,c) = nnz(PERFASYNCH3PP(j,i:i+4))/5;
      
          end
          
    end
    c = c+1;  
 end
 PerfSynch_slide(20,:) = [];  PerfAsynch_slide(20,:) = []; PerfAsynch3pp_slide(20,:) = [];

 figure(2)
subplot(1,3,1)

for i = 1:size(MyParticipants,1)-2
    plot(1:41,PerfSynch_slide(i,:),'-o')
    hold on
end
 plot(1:41,mean(PerfSynch_slide,1),'-o','LineWidth',2,'color','k')
 title('Synch')
 subplot(1,3,2)

for i = 1:size(MyParticipants,1)-2
    plot(1:41,PerfAsynch_slide(i,:),'-o')
    hold on
end
 plot(1:41,mean(PerfAsynch_slide,1),'-o','LineWidth',2,'color','k')
  title('Asynch')
 subplot(1,3,3)
 
 for i = 1:size(MyParticipants,1)-2
    plot(1:41,PerfAsynch3pp_slide(i,:),'-o')
    hold on
end
 plot(1:41,mean(PerfAsynch3pp_slide,1),'-o','LineWidth',2,'color','k')
 title('Asynch3pp')
 