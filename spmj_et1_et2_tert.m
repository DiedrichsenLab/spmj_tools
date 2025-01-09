function [et1, et2, tert] = spmj_et1_et2_tert(fmapDir, funcDir, subj_name)
% spmj_et1_et2_tert calculates MRI timing parameters for field map correction.
%
% This function reads the necessary MRI acquisition parameters from the
% JSON metadata associated with the subject's field map and first functional
% MRI run. It then computes and returns the first and second echo times
% (et1 and et2) in milliseconds, as well as the total effective EPI readout
% time (tert) considering GRAPPA acceleration.
%
% Inputs:
%   fmapDir: path to fmap in BIDS direcotry.
%   funcDir: path to func in BIDS direcotry.
%   subj_name: The subject's name as it appears in the DICOM server.
%
% Outputs:
%   et1: The first echo time in milliseconds.
%   et2: The second echo time in milliseconds.
%   tert: The total effective EPI readout time in milliseconds.
%
% by Marco Emanuele, March 2024

% fmapDir = fullfile(baseDir, "BIDS", subj_name, "fmap");
% funcDir = fullfile(baseDir, "BIDS", subj_name, "func");

phasediff = jsondecode(fileread(fullfile(fmapDir, ...
    ['sub-' subj_name '_phasediff.json'])));
func1st = jsondecode(fileread(fullfile(funcDir, ...
    ['sub-' subj_name '_task-task_run-01_bold.json'])));

%% compute tert
% total EPI readout time = = echo spacing (in ms) * base resolution 
% (also knows as number of echos). If you use GRAPPA acceleration, 
% you need to divide the total number of echos by two:
num_of_echoes = func1st.PhaseEncodingSteps;
echo_spacing = func1st.EffectiveEchoSpacing * 1000;                        % compute echo spacing in milliseconds
tert = echo_spacing * num_of_echoes;                                 % for GRAPPA sequence, EffectiveEchoSpacing is the echospacing divided by the acceleration factor

%% retrieve et1 et2 (in milliseconds)
et1 = phasediff.EchoTime1 * 1000;
et2 = phasediff.EchoTime2 * 1000;









