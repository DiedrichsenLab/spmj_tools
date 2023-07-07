function spmj_change_origin(P,origin)
% function spmj_change_origin(P,origin)
% Changes the origin of a list of header files to a new setting 
% Useful for viewing files in MRICro 
if (ischar(P))
    A{1}=P;
    P=A;
end;
for i=1:length(P)
    [hdr,swap]=spm_read_hdr(P{i});
    hdr.hist.origin(1:3)=origin;
    spmj_write_hdr(P{i},hdr,swap);    
end;
