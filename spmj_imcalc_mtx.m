function spmj_imcalc_mtx(P,Out,formula);
% function spmj_imcalc_mtx(P,Out,formula);
% Alias for imcalc in matric mode
% spm_imcalc_ui([],[],[],{1,[],[],[]});
if nargin<3
    formula=[];
end;
if nargin<2
    Out=[];
end;
if nargin<1
    P=[];
end;
spm_imcalc_ui(P,Out,formula,{1,[],[],[]});