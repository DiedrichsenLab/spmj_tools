function array=dotstr2array(instr)
% function array = dotstr2array(str)
% Makes a string using . to denote a list into an array of integers
% '1.2.4' -> [1,2,4]
elements = strsplit(instr, '.');
array = str2double(elements);
