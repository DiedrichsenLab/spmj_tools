function [et1, et2, tert] = spmj_et1_et2_tert(dataDir, subj_name, sn)
% spmj_et1_et2_tert calculates MRI timing parameters for field map correction.
%
% This function reads the necessary MRI acquisition parameters from the
% JSON metadata associated with the subject's field map and first functional
% MRI run. It then computes and returns the first and second echo times
% (et1 and et2) in milliseconds, as well as the total effective EPI readout
% time (tert) considering GRAPPA acceleration.
%
% Inputs:
%   dataDir: The root directory of the BIDS-structured MRI dataset.
%   subj_name: The subject's identifier (equal to subj_id in participants.tsv).
%   sn: The subject's number in participants.tsv.
%
% Outputs:
%   et1: The first echo time in milliseconds.
%   et2: The second echo time in milliseconds.
%   tert: The total effective EPI readout time in milliseconds.

fmapDir = fullfile(dataDir, "BIDS", subj_name, "fmap");
funcDir = fullfile(dataDir, "BIDS", subj_name, "func");

phasediff = jsondecode(fileread(fullfile(fmapDir, ...
    ['sub-' num2str(sn) '_phasediff.json'])));
func1st = jsondecode(fileread(fullfile(funcDir, ...
    ['sub-' num2str(sn) '_task-task_run-01_bold.json'])));

%% compute tert
% total EPI readout time = = echo spacing (in ms) * base resolution 
% (also knows as number of echos). If you use GRAPPA acceleration, 
% you need to divide the total number of echos by two:
base_resolution = func1st.BaseResolution;
echo_spacing = func1st.EffectiveEchoSpacing * 1000;                        % compute echo spacing in milliseconds
tert = echo_spacing * base_resolution / 2;

%% retrieve et1 et2 (in milliseconds)
et1 = phasediff.EchoTime1 * 1000;
et2 = phasediff.EchoTime2 * 1000;









