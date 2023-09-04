% clear all

dirct_optim = 'UIO Physiology Dropbox Dropbox/Lab Data/Andreas Lande/Analysis/Paper 1/timeVsPos_glm/paramOptim';
dirct_final = 'UIO Physiology Dropbox Dropbox/Lab Data/Andreas Lande/Analysis/Paper 1/timeVsPos_glm/GLM_all_lcss';
% Recordings to include
data.sessionIDs = { 'm1410_20191023_01', ... # 1 OLD (10-14 ish days training with cue setup)
                    'm1411_20191023_01', ... # 2 OLD (10-14 ish days training with cue setup)
                    'm1412_20191023_01', ... # 3 OLD (10-14 ish days training with cue setup)
                    'm1415_20191023_01', ... # 4 OLD (11 ish days training with cue setup)
                    'm1420_20200225_01', ... # 5 NEW (125 ym depth) (15 days training with cue setup)
                    'm1420_20200304_01', ... # 6 NEW (180 ym depth) (20 days training with cue setup)
                    'm1422_20200225_01', ... # 7 NEW (130 ym depth) (9 days training with cue setup)
                    'm1424_20200311_01', ... # 8 NEW (Posterior FOV) (6 days training with cue setup)
                    'm1424_20200315_01', ... # 9 NEW (Anterior FOV) (9 days training with cue setup) %'m1425_20200311_01', ... # 10 NEW (11 days training on wheel)
                    'm1438_20210606_01', ... # 11 NEW (5 days training with landmarks)
                    };
 
                
% run analysis with sb.glm.run_glm_pos if results are not available 
% run analysis with sb.glm.run_glm_time if results are not available 

 parameter_optimisation_pos()
 parameter_optimisation_time()
 compare_time_vs_pos()
                
               
