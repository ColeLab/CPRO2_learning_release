function quantify_rule_type_subject(subj_id)

%Author: Ravi Mill, rdm146@rutgers.edu
%Last Update: June 7th 2023

%Accompanies Mill & Cole (2023): "Neural representation dynamics reveal computational 
%principles of cognitive task learning"

%DESCRIPTION: Function quantifies the strength of compositional
%(single rules: logic, sensory, motor) and conjunctive (rule interaction term)
%representations in block-to-block regional activation patterns evoked by the 
%C-PRO2 task paradigm.
%Rule templates are created using held-out task data (i.e. averaging over relevant novel tasks 
%presented in the C-PRO2 Test session) and then fit via regression to the block-to-block 
%Practice session data. This yields beta coefficients capturing rule strength 
%for the 4 rule templates, which is repeated separately for each practiced task and block
%presentation, to quantify strength of each rule type over time.
%The method of quantifying the rule representations (estimating rule beta coefficients)
%is robust multiple linear regression. 
%The function also outputs behavioral associations (Spearman correlation)
%between the rule type detectability timecourses and behavioral
%accuracy/RT.

%INPUT:
%subj_id: numeric, id for subject being run; this sets paths to load in required
%   inputs: behavioral data (Excel sheet) and preprocessed vertex-level fMRI data for this subject

%OUTPUT:
%prac_rep: matlab structure containing timecourses of rule detectability for the practiced tasks,
%   and behavioral associations of those timecourses with accuracy and RT

%DEPENDENCIES:
%1. Connectome workbench must be installed and accessible to your system path, 
%   including the command line version (wb_command), 
%   https://www.humanconnectome.org/software/get-connectome-workbench
%2. (Included) Relevant files from Matlab Gifti toolbox
%   https://github.com/gllmflndn/gifti, allowing for opening (ciftiopen) and
%   saving (ciftisave and ciftisavereset) HCP CIFTI files
%3. (Included) CAB-NP network partition .dlabel file from https://github.com/ColeLab/ColeAnticevicNetPartition
%   that links vertex/voxel data to affiliated Glasser atlas regions:
%   CortexSubcortex_ColeAnticevic_NetPartition_wSubcorGSR_parcels_LR.dlabel.nii

%See section below for paths and parameters, *=paths to the data set by the user
%See lines 145-173 for manual treatment of subjects with missing/partial runs

%% Set up paths, directories etc

%*set path to base directory (where input data are stored and where output data will be saved)
basedir='/projectsn/f_mc1689_1/CPRO2_learning/';

%list of all subjects (for information's sake)
subjects=[1,2,3,5,6,7,9,10,11,13,14,16,17,18,21,22,23,24,26,29,30,31,32,34,35,36,37,38,39,41,42,43,44,45,46,48,49,50,51,53,54,55,56,57];
num_subs=length(subj_id);

%*set path to behavioral data (excel file containing Eprime output for all subjects)
%used here to identify tasks and rules for each block across Practiced and Test sessions
behav_xls=[basedir,'/data/results/CPRO2_fMRI_EEG_main.xlsx'];

%set number of regions (based on CAB-NP network partition)
num_regions=718; %360 surface, 358 subcortex

%set path to directory with preprocessed fMRI data (vertex/voxel-level, in hdf5
%format) for each subject
%num total vertices/voxels=91282
vertex_dir=[basedir,'/data/preprocessed/MRI/vertexWise/'];
vertex_suff='_glmOutput_vertex_data.h5'; %suffix for each subject file

%input file keys used to reference h5 files with post-GLM output (activity betas)
%note that constant is always the first row
%120 task regs (one for each Test task block; entire miniblock)
model1_test_key='/taskRegression/model1_Task_Test_24pXaCompCorXVolterra_taskReg_betas_canonical';
%144 task regs (one for each Prac task block; entire miniblock)
model2_prac_key='/taskRegression/model2_Task_Practice_24pXaCompCorXVolterra_taskReg_betas_canonical';

%set path to CIFTI file containing labels that affiliate each vertex/voxel to a Glasser atlas region 
glasser_dlabel='CortexSubcortex_ColeAnticevic_NetPartition_wSubcorGSR_parcels_LR.dlabel.nii';
cii_label=ciftiopen(glasser_dlabel,'wb_command');
CIFTIlabel=cii_label.cdata;
%find unique region labels
region_labels=unique(CIFTIlabel);

%task info
num_prac_blocks=144;
num_test_blocks=120;
blocks_per_run=12;
num_prac_tasks=4;

%set col numbers that will reference relevant conditions, behavior in the excel sheet
%note all columns have same ordering across Practice/Test session data
sub_col=2; %sub ID
logic_col=117; %rule info
motor_col=118;
sensory_col=143;
cond_col=102; %determines practiced (Switch1, Switch2, NoSwitch1, NoSwitch2) versus novel tasks

acc_col=120; %behavioral acc
rt_col=132;

%output directory
outputdir=[basedir,'/data/results/Prac_Conjunctive_robust_publicRelease/'];
if ~exist(outputdir);mkdir(outputdir);end
%where subjects' data will be saved
subdir=[outputdir,'/Regionwise/'];
if ~exist(subdir);mkdir(subdir);end

%% load in regionwise behav/fMRI data for prac/test (as in other scripts)

tic;

%load in behavior for Practice session
[~,~,group_edat_prac]=xlsread(behav_xls,'fMRI_prac');
group_edat_prac=group_edat_prac(2:end,1:161);%exclude first row header

%load in behavior for Test session 
[~,~,group_edat_test]=xlsread(behav_xls,'fMRI_test');
group_edat_test=group_edat_test(2:end,1:161);%exclude first row header

%initialize output variables
%rule strength
prac_rep=struct;
prac_rep.sensory_r=zeros(num_prac_tasks,num_prac_blocks/num_prac_tasks,num_regions,num_subs);
prac_rep.logic_r=zeros(num_prac_tasks,num_prac_blocks/num_prac_tasks,num_regions,num_subs);
prac_rep.motor_r=zeros(num_prac_tasks,num_prac_blocks/num_prac_tasks,num_regions,num_subs);
prac_rep.conj_r=zeros(num_prac_tasks,num_prac_blocks/num_prac_tasks,num_regions,num_subs);
%spearman correlation of rule strength x behavioral accuracy
prac_rep.sensory_acc=zeros(num_prac_tasks,num_regions,num_subs);
prac_rep.logic_acc=zeros(num_prac_tasks,num_regions,num_subs);
prac_rep.motor_acc=zeros(num_prac_tasks,num_regions,num_subs);
prac_rep.conj_acc=zeros(num_prac_tasks,num_regions,num_subs);
%spearman correlation of rule strength x behavioral reaction time (RT)
prac_rep.sensory_rt=zeros(num_prac_tasks,num_regions,num_subs);
prac_rep.logic_rt=zeros(num_prac_tasks,num_regions,num_subs);
prac_rep.motor_rt=zeros(num_prac_tasks,num_regions,num_subs);
prac_rep.conj_rt=zeros(num_prac_tasks,num_regions,num_subs);

%store rt on its own
prac_rep.rt=zeros(num_prac_tasks,num_prac_blocks/num_prac_tasks,num_subs);

for i=1:length(subj_id)
    %*prac data
    sub_h5=[vertex_dir,'sub-',num2str(subj_id),vertex_suff];
    sub_data=h5read(sub_h5,model2_prac_key);
    sub_data(1,:)=[]; %remove constant
    %check for subs with partial data
    if subj_id==9
        nan_block=nan(4*blocks_per_run,size(sub_data,2)); %sub9=missing prac runs=9-12
        sub_data=[sub_data;nan_block];
    elseif subj_id==44
        nan_block=nan(1*blocks_per_run,size(sub_data,2)); %sub44=missing prac runs 4, and 11-12
        nan_block2=nan(2*blocks_per_run,size(sub_data,2));
        sub_data=[sub_data(1:(3*blocks_per_run),:);nan_block;sub_data((3*blocks_per_run)+1:end,:);nan_block2];
    end
    sub_prac_data=sub_data;
    clear sub_data;

    %test data
    sub_data=h5read(sub_h5,model1_test_key);
    sub_data(1,:)=[]; %remove constant
    %check for subs with partial data
    if subj_id==2
        nan_block=nan(1*blocks_per_run,size(sub_data,2)); %sub2=missing test run 1
        sub_data=[nan_block;sub_data];
    elseif subj_id==6
        nan_block=nan(2*blocks_per_run,size(sub_data,2)); %sub6=missing test runs 9,10
        sub_data=[sub_data;nan_block];
    elseif subj_id==49
        nan_block=nan(1*blocks_per_run,size(sub_data,2)); %sub49=missing test run 1
        sub_data=[nan_block;sub_data];
    end
    sub_test_data=sub_data;
    clear sub_data;    
    
    %pull out cond/behav data for this sub - Prac session
    sub_prac_ind=find(cell2mat(group_edat_prac(:,sub_col))==subj_id);
    sub_edat_prac=group_edat_prac(sub_prac_ind,:);
    
    %pull out cond/behav data for this sub - Test session
    sub_test_ind=find(cell2mat(group_edat_test(:,sub_col))==subj_id);
    sub_edat_test=group_edat_test(sub_test_ind,:);
    
    %pull out rules for practiced tasks
    %*note order practiced tasks from 1-4 is based on chronological order 
    %that they appeared to subjects in their Practice session, to stay consistent
    %with how they're ordered in the fMRI data
    logics=sub_edat_prac(:,logic_col);
    un_logics=unique(logics,'stable'); %stable=preserve chronological ordering
    sensorys=sub_edat_prac(:,sensory_col);
    un_sensorys=unique(sensorys,'stable');
    motors=sub_edat_prac(:,motor_col);
    un_motors=unique(motors,'stable');
    prac_tasks=[un_logics,un_sensorys,un_motors];
    
    %code tasks per block (only need one rule e.g. logic for this)
    trial_cond=sub_edat_prac(:,logic_col);
    block_cond={};
    trial_count=1;
    for b=1:num_prac_blocks
        block_cond=[block_cond;trial_cond{trial_count}];
        trial_count=trial_count+3;
    end
    
     %pull out behavioral accuracy and RT for Prac session
    sub_acc=cell2mat(sub_edat_prac(:,acc_col));
    sub_rt=sub_edat_prac(:,rt_col);
    %need to convert char nans to number - these appear in RT column as
    %they were manually added in the spreadsheet
    cr_double=zeros(size(sub_rt));
    cr_number=cellfun(@(x) isnumeric(x),sub_rt);
    cr_double(cr_number)=[sub_rt{cr_number}];
    cr_double(~cr_number)=NaN;
    sub_rt=cr_double; 
    
    %average over 3 trials within each block
    sub_acc_block=[];
    sub_rt_block=[];
    trial_count=1;
    for b=1:num_prac_blocks
        sub_acc_block=[sub_acc_block;nanmean(sub_acc(trial_count:trial_count+2))];
        sub_rt_block=[sub_rt_block;nanmean(sub_rt(trial_count:trial_count+2))];
        trial_count=trial_count+3;
    end
    
    %separate blockwise data for each practiced task
    prac_block_data=struct;
    prac_block_acc=struct;
    prac_block_rt=struct;
    for t=1:num_prac_tasks
        ind=strcmp(block_cond,un_logics{t});
        prac_block_data.(['Task',num2str(t)])=sub_prac_data(ind,:);
        %acc
        prac_block_acc.(['Task',num2str(t)])=sub_acc_block(ind,:);
        %rt
        prac_block_rt.(['Task',num2str(t)])=sub_rt_block(ind,:);
        
        %also create separate rt timecourse (for later analyses)
        prac_rep.rt(t,:,i)=sub_rt_block(ind,:)';
    end
    
    %pick out indices for novel tasks in the Test session (used to create
    %rule templates that are fit to the Prac session data)
    test_cond_trial=sub_edat_test(:,cond_col);
    test_cond_block=[];
    trial_count=1;
    for b=1:num_test_blocks
        test_cond_block=[test_cond_block;test_cond_trial(trial_count)];
        trial_count=trial_count+3;
    end
    novel_block_inds=find(strcmp(test_cond_block,'Novel')==1);
    
    %pull out rules for tasks in the Test session
    logics=sub_edat_test(:,logic_col);
    sensorys=sub_edat_test(:,sensory_col);
    motors=sub_edat_test(:,motor_col);
    test_tasks=[logics,sensorys,motors];
    
    %convert task info from trial to block
    test_blocks=[];
    trial_count=1;
    for b=1:num_test_blocks
        test_blocks=[test_blocks;test_tasks(trial_count,:)];
        trial_count=trial_count+3;
    end
    
    %pull out novel block rules
    novel_blocks=test_blocks(novel_block_inds,:);
    
    %pull out novel task fmri betas
    novel_data=sub_test_data(novel_block_inds,:);
    
   %% RULE TYPE QUANTIFICATION 
   
   %loop through regions, to quantify strength of rule type templates 
   %(estimated using novel tasks in the Test session) for each practiced
   %task over time (blocks) in the Prac session
   %also compute associations between resultant rule type timecourses and
   %behavior
   for j=1:num_regions
        %pick out vertex/voxel inds corresponding to this CAB-NP region
        region_num=region_labels(j);
        vert_inds=find(CIFTIlabel==region_num);
                
        %loop through prac tasks
        for tt=1:size(prac_tasks,1)
            %pull out rule
            tt_rules=prac_tasks(tt,:);
            
            %compute rule templates using Test data
            rule_types={'logic','sensory','motor'};
            for rr=1:length(rule_types)

                %find Test session blocks with same rule
                rr_rule=tt_rules{rr};
                inds=strcmp(novel_blocks(:,rr),rr_rule);
                inds=find(double(inds)==1);
                
                %create template by averaging over Test blocks
                temp=nanmean(novel_data(inds,vert_inds),1);
                
                %store in template var 
                eval([rule_types{rr},'_temp=temp;']); %template
                
            end
            
            %conj
            %create conjunction template by multiplying three rule
            %templates together
            conj_temp=logic_temp.*sensory_temp.*motor_temp;
            
            %loop through prac blocks and quantify rule template strength
            %via multiple regression
            for bb=1:(num_prac_blocks/num_prac_tasks)
                y=prac_block_data.(['Task',num2str(tt)])(bb,vert_inds)';
                x=[logic_temp',sensory_temp',motor_temp',conj_temp',ones(length(conj_temp),1)];
                
                %exclude subcortical regions with low number of vertices (these will generate errors)
                if length(conj_temp)<7
                    b=[NaN;NaN;NaN;NaN;NaN]; %7 is the cutoff for min vertices in most regionwise analyses
                else
                    if sum(isnan(y))==length(y)
                        %store NaN betas for subjects with missing Prac blocks (sub 9 and sub 44)
                        b=[NaN;NaN;NaN;NaN;NaN];
                    else
                        %estimate betas via robust regression
                        b=robustfit(x,y,'bisquare',4.685,'off');
                    end
                end
                
                %store betas for all 4 rule terms
                prac_rep.logic_r(tt,bb,j,i)=b(1);
                prac_rep.sensory_r(tt,bb,j,i)=b(2);
                prac_rep.motor_r(tt,bb,j,i)=b(3);
                prac_rep.conj_r(tt,bb,j,i)=b(4);
                
            end

            %BEHAVIORAL ASSOCIATIONS
            %compute Spearman correlation with prac acc/rt
            
            %pull out behavioral timecourses for this prac task
            tt_acc=prac_block_acc.(['Task',num2str(tt)]);
            tt_rt=prac_block_rt.(['Task',num2str(tt)]);
            
            %acc
            r=corr([prac_rep.logic_r(tt,:,j,i)',tt_acc],'Type','Spearman','Rows','Complete');
            prac_rep.logic_acc(tt,j,i)=r(1,2);
            r=corr([prac_rep.sensory_r(tt,:,j,i)',tt_acc],'Type','Spearman','Rows','Complete');
            prac_rep.sensory_acc(tt,j,i)=r(1,2);
            r=corr([prac_rep.motor_r(tt,:,j,i)',tt_acc],'Type','Spearman','Rows','Complete');
            prac_rep.motor_acc(tt,j,i)=r(1,2);
            r=corr([prac_rep.conj_r(tt,:,j,i)',tt_acc],'Type','Spearman','Rows','Complete');
            prac_rep.conj_acc(tt,j,i)=r(1,2);
            
            %rt
            r=corr([prac_rep.logic_r(tt,:,j,i)',tt_rt],'Type','Spearman','Rows','Complete');
            prac_rep.logic_rt(tt,j,i)=r(1,2);
            r=corr([prac_rep.sensory_r(tt,:,j,i)',tt_rt],'Type','Spearman','Rows','Complete');
            prac_rep.sensory_rt(tt,j,i)=r(1,2);
            r=corr([prac_rep.motor_r(tt,:,j,i)',tt_rt],'Type','Spearman','Rows','Complete');
            prac_rep.motor_rt(tt,j,i)=r(1,2);
            r=corr([prac_rep.conj_r(tt,:,j,i)',tt_rt],'Type','Spearman','Rows','Complete');
            prac_rep.conj_rt(tt,j,i)=r(1,2);
        end
        
         disp(['time taken for region',num2str(j),'=',num2str(toc)]);
        
   end
   
   %save
   sub_file=[subdir,'sub',num2str(subj_id),'.mat'];
   save(sub_file,'prac_rep','-v7.3');
   
   disp(['time taken for subject',num2str(i),'=',num2str(toc)]);
    
end


end