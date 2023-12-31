function SPM=spmj_basic_design(nr,P,SPMBATCH,opt,anadir)
% function SPM=spmj_basic_design(nr,P,SPMBACTH,opt,anadir)
%
% nr: specifier for analysis type
%   1:One sample t-test                  
%   2:Two sample t-test                  
%   3:Paired t-test                      
%   4:One way Anova                      
%   5:One way Anova (with constant)      
%   6:One way Anova (Within-subjects)    
%   7:Simple regression (correlation)    
%   8:Multiple regression                
%   9:Multiple regression (with constant)
%   10:AnCova   
% 
% SPMBATCH: Struct with fields: 
%   n: Number of levels (1*4 vector)
%   C: Matrix of covariates 
%   Cnames: Names of covariates 
%   masks: cell array of names of mask images 
%
% opt: options for analysis
% opt(1):   Grand Mean Scaling: 1=Scaling of overall grand mean 9=no GMS
% opt(2):   Masking: 1 = none, 2= absolute, 3 = relative
% opt(3):   explicitely mask images: 0 =no, 1 = yes
% opt(4):   Global calculation: 1=omit
% opt(5):   Non Sphericity correction: 0=yes, 1=no
% opt(6):   Replications over: 2=subjects, 1 = conditions
% opt(7):   Correlated repeated measures: 1=yes, 0=no


%-----------------------------------------------------------------------

SCCSid  = '2.49';
SPMid = spm('FnBanner',mfilename,SCCSid);
[Finter,Fgraph,CmdLine] = spm('FnUIsetup','Stats: Setup analysis',0);
spm_help('!ContextHelp',mfilename)
Designs=spm_spm_ui('DesDefs_Stats');    %get all group design structs from SPM

D=Designs(nr);                %select within-subject ANOVA
% subj_dirs=SPMBATCH.dirs;
%=======================================================================
% - C O N F I G U R E   D E S I G N
%=======================================================================
spmj_init_groupana;
%-Set factor names for this design
%-----------------------------------------------------------------------
sCC      = sf_estrrep(sCC,[my_sF',D.sF']);
sCFI     = sf_estrrep(sCFI,[my_sF',D.sF']);
sGloNorm = sf_estrrep(sGloNorm,[my_sF',D.sF']);
sGMsca   = sf_estrrep(sGMsca,[my_sF',D.sF']);

%-Get filenames & factor indicies
%-----------------------------------------------------------------------
% D.n=[1 length(img_nrs) length(subj_dirs) 1];
% AUTOFILL.subj_dirs=subj_dirs;
% AUTOFILL.img_nrs=img_nrs;
% AUTOFILL.img_mask=img_mask;
% [P,I]    = get_files_indices(AUTOFILL,D.sF,D.n,D.b.aTime);
if (isempty(P))
    P=spm_select(inf,'image','select_scans',[],anadir);
    P=cellstr(P);
end;
D.n=SPMBATCH.n;
nScan    = length(P);						%-#obs
I=ones(nScan,4);                            % for now 
%-Additional design parameters
%-----------------------------------------------------------------------
bL       = any(diff(I,1),1); 	%-Multiple factor levels?
	    % NB: bL(2) might be thrown by user specified f1 levels
	    %     (D.b.aTime & D.n(2)>1) - assumme user is consistent?
bFI      = [bL(1),bL(2:3)&~bL(4),bL(4),bL([2,3])&bL(4)];
	    %-Allowable interactions for covariates
	    %-Only offer interactions with multi-level factors, and
	    % don't offer by F2|F3 if bL(4)!

%-Build Condition (H) and Block (B) partitions
%=======================================================================
eval(['[H,Hnames] = spm_DesMtx(',D.Hform,');'])
if rank(H)==nScan, error('unestimable condition effects'), end
eval(['[B,Bnames] = spm_DesMtx(',D.Bform,');'])
if rank(B)==nScan, error('unestimable block effects'), end

%-Drop a constant H partition if B partition can model constant
if size(H,2)>0 & all(H(:)==1) & (rank([H B])==rank(B))
	H = []; Hnames = {};
	warning('Dropping redundant constant H partition')
end


%-Covariate partition(s): interest (C) & nuisance (G) excluding global
%=======================================================================
if (isfield(SPMBATCH,'C'))
    nC = [size(SPMBATCH.C,2) 0];			%-Default #covariates
    C  = {SPMBATCH.C,[]}; Cnames = {SPMBATCH.Cnames,{}};	%-Covariate DesMtx partitions & names
    xC = [];			%-Struct array to hold raw covariates
else
    nC = D.nC;			%-Default #covariates
    C  = {[],[]}; Cnames = {{},{}};	%-Covariate DesMtx partitions & names
    xC = [];			%-Struct array to hold raw covariates
end;


dcname = {'CovInt','NusCov'};	%-Default root names for covariates
dstr   = {'covariate','nuisance variable'};

GUIpos = spm_input('!NextPos');
nc     = [0,0];
for i  = 1:2			% 1:covariates of interest, 2:nuisance variables

    if isinf(nC(i)), nC(i)=spm_input(['# ',dstr{i},'s'],GUIpos,'w1'); end

    while nc(i) < nC(i)

	%-Create prompt, get covariate, get covariate name
        %---------------------------------------------------------------
	    if nC(i)==1, str=dstr{i}; else, str=sprintf('%s %d',dstr{i},nc(i)+1); end
        try
            c=C{i}(:,nc(i)+1);
        catch 
            c = spm_input(str,GUIpos,'r',[],[nScan,Inf]);
            if any(isnan(c(:))), break, end		%-NaN is dummy value to exit
        end;
	    nc(i)  = nc(i)+1;			%-#Covariates (so far)
	    if nC(i)>1,	tstr = sprintf('%s^{%d}',dcname{i},nc(i));
	    else,		tstr = dcname{i}; end
        try
            cname=Cnames{i}{nc(i)};
        catch     
            cname  = spm_input([str,' name?'],'+1','s',tstr);
        end;
       	rc     = c;				%-Save covariate value
	    rcname = cname;				%-Save covariate name

        %-Interaction option? (if single covariate vector entered)?
        %---------------------------------------------------------------
        if size(c,2) == 1
       	    if length(D.iCFI{i})>1
       		%-User choice of interaction options, default is negative
       		%-Only offer interactions for appropriate factor combinations
		iCFI = intersect(abs(D.iCFI{i}),find([1,bFI]));
		dCFI = max([1,intersect(iCFI,-D.iCFI{i}(D.iCFI{i}<0))]);
        	iCFI = spm_input([str,': interaction?'],'+1','m',...
			sCFI(iCFI),iCFI,find(iCFI==dCFI));
	    else
		iCFI = abs(D.iCFI{i});		%-AutoSelect default option
	    end
	else
	    iCFI = 1;
	end

        %-Centre covariate(s)? (Default centring to correspond to CFI)
        % Always offer "no centering" as default for design matrix blocks
        %---------------------------------------------------------------
	DiCC = D.iCC{i};
	if size(c,2)>1, DiCC = union(DiCC,-8); end
        if length(DiCC)>1
        	%-User has a choice of centering options
		%-Only offer factor specific for appropriate factor combinations
		iCC = intersect(abs(DiCC),find([1,bFI,1]) );
        	%-Default is max -ve option in D, overridden by iCFI if CFI
		if iCFI == 1, dCC = -DiCC(DiCC<0); else, dCC = iCFI; end
		dCC = max([1,intersect(iCC,dCC)]);
		iCC = spm_input([str,': centre?'],'+1','m',...
			sCC(iCC),iCC,find(iCC==dCC));
        else
        	iCC = abs(DiCC);	%-AutoSelect default option
        end
	%-Centre within factor levels as appropriate
        if any(iCC == [1:7]), c = c - spm_meanby(c,eval(CCforms{iCC})); end

        %-Do any interaction (only for single covariate vectors)
        %---------------------------------------------------------------
        if iCFI > 1				%-(NB:iCFI=1 if size(c,2)>1)
       		tI        = [eval(CFIforms{iCFI,1}),c];
		tConst    = CFIforms{iCFI,2};
		tFnames   = [eval(CFIforms{iCFI,3}),{cname}];
		[c,cname] = spm_DesMtx(tI,tConst,tFnames);
	elseif size(c,2)>1			%-Design matrix block
		[null,cname] = spm_DesMtx(c,'X',cname);
	else
		cname = {cname};
	end

	%-Store raw covariate details in xC struct for reference
	%-Pack c into appropriate DesMtx partition
        %---------------------------------------------------------------
	%-Construct description string for covariate
	str = {sprintf('%s: %s',str,rcname)};
	if size(rc,2)>1, str = {sprintf('%s (block of %d covariates)',...
		str{:},size(rc,2))}; end
	if iCC < 8, str=[str;{['used centered ',sCC{iCC}]}]; end
	if iCFI> 1, str=[str;{['fitted as interaction ',sCFI{iCFI}]}]; end

	tmp       = struct(	'rc',rc,	'rcname',rcname,...
				'c',c,		'cname',{cname},...
				'iCC',iCC,	'iCFI',iCFI,...
				'type',i,...
				'cols',[1:size(c,2)] + ...
						size([H,C{1}],2) +  ...
						size([B,C{2}],2)*(i-1),...
				'descrip',{str}				);
	if isempty(xC), xC = tmp; else, xC = [xC,tmp]; end
	C{i}(:,nc(i))=c; 
	Cnames{i}{nc(i),1} = cname{1};

    end	% (while)

end % (for)
clear c tI tConst tFnames
spm_input('!SetNextPos',GUIpos);

%-Unpack into C & G design matrix sub-partitions
G = C{2}; Gnames = Cnames{2};
C = C{1}; Cnames = Cnames{1};


%-Options...
%=======================================================================
%-Global normalization options                                 (GloNorm)
%-----------------------------------------------------------------------
if length(D.iGloNorm)>1
	%-User choice of global normalisation options, default is negative
	%-Only offer factor specific for appropriate factor combinations
	iGloNorm = intersect(abs(D.iGloNorm),find([1,bFI,1,1]));
	dGloNorm = max([0,intersect(iGloNorm,-D.iGloNorm(D.iGloNorm<0))]);
	iGloNorm = spm_input('GloNorm: Select global normalisation','+1','m',...
	    	sGloNorm(iGloNorm),iGloNorm,find(iGloNorm==dGloNorm));
else
	iGloNorm = abs(D.iGloNorm);
end


%-Grand mean scaling options                                     (GMsca)
%-----------------------------------------------------------------------
if iGloNorm==8
	iGMsca=8;	%-grand mean scaling implicit in PropSca GloNorm
elseif length(D.iGMsca)==1
	iGMsca = abs(D.iGMsca);
else
	%-User choice of grand mean scaling options
	%-Only offer factor specific for appropriate factor combinations
	iGMsca = intersect(abs(D.iGMsca),find([1,bFI,0,1]));
        %-Default is max -ve option in D, overridden by iGloNorm if AnCova
        if iGloNorm==9, dGMsca=-D.iGMsca(D.iGMsca<0); else, dGMsca=iGloNorm; end
	dGMsca = max([0,intersect(iGMsca,dGMsca)]);
	%iGMsca = spm_input('GMsca: grand mean scaling','+1','m',...
	%    	sGMsca(iGMsca),iGMsca,find(iGMsca==dGMsca));
    iGMsca = opt(1);
end


%-Value for PropSca / GMsca                                         (GM)
%-----------------------------------------------------------------------
if iGMsca == 9                          %-Not scaling (GMsca or PropSca)
	GM = 0;                         %-Set GM to zero when not scaling
else                                    %-Ask user value of GM
	if iGloNorm==8
		str = 'PropSca global mean to';
	else
		str = [strrep(sGMsca{iGMsca},'scaling of','scale'),' to'];
	end
	GM = spm_input(str,'+1','r',D.GM,1);
	%-If GM is zero then don't GMsca! or PropSca GloNorm
	if GM==0, iGMsca=9; if iGloNorm==8, iGloNorm=9; end, end
end

%-Sort out description strings for GloNorm and GMsca
%-----------------------------------------------------------------------
sGloNorm = sGloNorm{iGloNorm};
sGMsca   = sGMsca{iGMsca};
if iGloNorm==8
	sGloNorm = sprintf('%s to %-4g',sGloNorm,GM);
elseif iGMsca<8
	sGMsca   = sprintf('%s to %-4g',sGMsca,GM);
end


%-Global centering (for AnCova GloNorm)                             (GC)
%-----------------------------------------------------------------------
%-Specify the centering option for the global covariate for AnCova
%-Basically, if 'GMsca'ling then should centre to GM (iGC=11). Otherwise,
% should centre in similar fashion to AnCova (i.e. by the same factor(s)),
% such that models are seperable (iGC=10). This is particularly important
% for subject specific condition effects if then passed on to a second-level
% model. (See also spm_adjmean_ui.m) SPM96 (& earlier) used to just centre
% GX around its (overall) mean (iGC=1).

%-This code allows more general options to be specified (but is a bit complex)
%-Setting D.iGC=[-10,-11] gives the standard choices above

%-If not doing AnCova then GC is irrelevant
if ~any(iGloNorm == [1:7])
	iGC = 12;
	gc  = [];
else
	%-Annotate options 10 & 11 with specific details
	%---------------------------------------------------------------
	%-Tag '(as implied by AnCova)' with actual AnCova situation
	sCC{10} = [sCC{iGloNorm},' (<= ',sGloNorm,')'];
	%-Tag 'GM' case with actual GM & GMsca case
	sCC{11} = sprintf('around GM=%g (i.e. %s after grand mean scaling)',...
		GM,strrep(sCC{iGMsca},'around ',''));

	%-Constuct vector of allowable iGC
	%---------------------------------------------------------------
	%-Weed out redundent factor combinations from pre-set allowable options
	iGC = intersect(abs(D.iGC),find([1,bFI,1,1,1,1]));
	%-Omit 'GM' option if didn't GMsca (iGMsca~=8 'cos doing AnCova)
	if any(iGMsca==[8,9]), iGC = setdiff(iGC,11); end
	%-Omit 'GM' option if same as '(as implied by AnCova)'
	if iGloNorm==iGMsca, iGC = setdiff(iGC,11); end

	%-If there's a choice, set defaults (if any), & get answer
	%---------------------------------------------------------------
	if length(iGC)>1
		dGC = max([0,intersect(iGC,-D.iGC(D.iGC<0))]);
		str = 'Centre global covariate';
		if iGMsca<8, str = [str,' (after grand mean scaling)']; end
		iGC = spm_input(str,'+1','m',sCC(iGC),iGC,find(iGC==dGC));
	elseif isempty(iGC)
		error('Configuration error: empty iGC')
	end

	%-If 'user specified' then get value
	%---------------------------------------------------------------
	if iGC==9
		gc     = spm_input('Centre globals around','+0','r',D.GM,1);
		sCC{9} = sprintf('%s of %g',sCC{iGC},gc);
	else
		gc  = 0;
	end
end


%-Thresholds & masks defining voxels to analyse                   (MASK)
%=======================================================================
GUIpos = spm_input('!NextPos');

%-Analysis threshold mask
%-----------------------------------------------------------------------
%-Work out available options:
% -Inf=>None, real=>absolute, complex=>proportional, (i.e. times global)
M_T = D.M_.T; if isempty(M_T), M_T = [-Inf, 100, 0.8*sqrt(-1)]; end
M_T = {	'none',		M_T(min(find(isinf(M_T))));...
	'absolute',	M_T(min(find(isfinite(M_T)&(M_T==real(M_T)))));...
	'relative',	M_T(min(find(isfinite(M_T)&(M_T~=real(M_T)))))	};

%-Work out available options
%-If there's a choice between proportional and absolute then ask
%-----------------------------------------------------------------------
q = ~[isempty(M_T{1,2}), isempty(M_T{2,2}), isempty(M_T{3,2})];

if all(q(2:3))
	tmp = opt(2); %spm_input('Threshold masking',GUIpos,'b',M_T(q,1),find(q));
	q(setdiff([1:3],tmp))=0;
end

%-Get mask value - note that at most one of q(2:3) is true
%-----------------------------------------------------------------------
if ~any(q)				%-Oops - nothing specified!
	M_T = -Inf;
elseif all(q==[1,0,0])			%-no threshold masking
	M_T = -Inf;
else					%-get mask value
	if q(1),	args = {'br1','None',-Inf,abs(M_T{1+find(q(2:3)),2})};
	else,		args = {'r',abs(M_T{1+find(q(2:3)),2})}; end
	if q(2)
		M_T = spm_input('threshold',GUIpos,args{:});
	elseif q(3)
		M_T = spm_input('threshold (relative to global)',GUIpos,...
								args{:});
		if isfinite(M_T) & isreal(M_T), M_T=M_T*sqrt(-1); end
	else
		error('Shouldn''t get here!')
	end
end

%-Make a description string
%-----------------------------------------------------------------------
if isinf(M_T)
	xsM.Analysis_threshold = 'None (-Inf)';
elseif isreal(M_T)
	xsM.Analysis_threshold = sprintf('images thresholded at %6g',M_T);
else
	xsM.Analysis_threshold = sprintf(['images thresholded at %6g ',...
		'times global'],imag(M_T));
end


%-Implicit masking: Ignore zero voxels in low data-types?
%-----------------------------------------------------------------------
% (Implicit mask is NaN in higher data-types.)
type = getfield(spm_vol(P{1,1}),'dt');
type = type(1); 
if ~spm_type(type,'nanrep')
	switch D.M_.I
	case Inf,    M_I = spm_input('Implicit mask (ignore zero''s)?',...
			'+1','y/n',[1,0],1);		%-Ask
	case {0,1}, M_I = D.M_.I;			%-Pre-specified
	otherwise,  error('unrecognised D.M_.I type')
	end

	if M_I, xsM.Implicit_masking = 'Yes: zero''s treated as missing';
	else,   xsm.Implicit_masking = 'No'; end
else
	M_I = 1;
	xsM.Implicit_masking = 'Yes: NaN''s treated as missing';
end


%-Explicit mask images (map them later...)
%-----------------------------------------------------------------------
switch(D.M_.X)
case Inf,   M_X = opt(3); %spm_input('explicitly mask images?','+1','y/n',[1,0],2);
case {0,1}, M_X = D.M_.X;
otherwise,  error('unrecognised D.M_.X type')
end
% if M_X, M_P = spm_get(Inf,'*.img',{'select mask images'}); else, M_P = {}; end
if M_X, M_P = SPMBATCH.masks; else, M_P = {}; end


%-Global calculation                                            (GXcalc)
%=======================================================================
iGXcalc = abs(D.iGXcalc);
%-Only offer "omit" option if not doing any GloNorm, GMsca or PropTHRESH
if ~(iGloNorm==9 & iGMsca==9 & (isinf(M_T)|isreal(M_T)))
	iGXcalc = intersect(iGXcalc,[2:size(sGXcalc,1)]);
end
if isempty(iGXcalc)
	error('no GXcalc options')
elseif length(iGXcalc)>1
	%-User choice of global calculation options, default is negative
	dGXcalc = max([1,intersect(iGXcalc,-D.iGXcalc(D.iGXcalc<0))]);
	iGXcalc = opt(4); %spm_input('Global calculation','+1','m',sGXcalc(iGXcalc),iGXcalc,find(iGXcalc==dGXcalc));
else
	iGXcalc = abs(D.iGXcalc);
end

if iGXcalc==2				%-Get user specified globals
	g = spm_input('globals','+0','r',[],[nScan,1]);
end
sGXcalc = sGXcalc{iGXcalc};


% Non-sphericity correction
%=======================================================================

% if there are multilevel factors, ask for correction
%-----------------------------------------------------------------------
if length(find(max(I) > 1)) > 1
	xVi.iid  = opt(5); %spm_input('non-sphericity correction?','+1','y/n',[0,1],1);
else
	xVi.iid  = 1;
end


if xVi.iid

	% i.i.d. assumptions where xVi.V = 1
	%---------------------------------------------------------------
	xVi.V    = speye(nScan);

else
	% otherwise, we have repeated measures design 
	%===============================================================
	nL      = max(I);		% number of levels
	mL      = find(nL > 1);		% multilevel factors
	xVi.I   = I;
	xVi.sF  = D.sF;
	xVi.var = sparse(1,4);
	xVi.dep = sparse(1,4);


	% eliminate replication factor from mL
	%---------------------------------------------------------------
	for i = 1:4
		mstr{i} = sprintf('%s (%i)',D.sF{i},nL(i));
	end
	str   = 'replications are over?';
	rep   = opt(6);     %spm_input(str,'+1','m',mstr(mL),1:length(mL));

	% and ask whether repeated measures are independent
	%---------------------------------------------------------------
	str   = 'correlated repeated measures';
	dep   = opt(7); %spm_input(str,'+1','b',{'yes','no'},[1 0],1);


	%-Place covariance components Q{:} in xVi.Vi
	%---------------------------------------------------------------
	mL(rep)     = [];
	xVi.var(mL) = 1;
	xVi.dep(mL) = dep;
	xVi         = spm_non_sphericity(xVi);

end


%=======================================================================
% - C O N F I G U R E   D E S I G N
%=======================================================================
spm('FigName','Stats: configuring',Finter,CmdLine);
spm('Pointer','Watch');


%-Images & image info: Map Y image files and check consistency of
% dimensions and orientation / voxel size
%=======================================================================
fprintf('%-40s: ','Mapping files')                                   %-#
VY    = spm_vol(char(P));


%-Check compatability of images (Bombs for single image)
%-----------------------------------------------------------------------
if any(any(diff(cat(1,VY(:).dim),1,1),1)) 
	error('images do not all have the same dimensions')
end
if any(any(any(diff(cat(3,VY(:).mat),1,3),3)))
	error('images do not all have same orientation & voxel size')
end

fprintf('%30s\n','...done')                                          %-#


%-Global values, scaling and global normalisation
%=======================================================================
%-Compute global values
%-----------------------------------------------------------------------
switch iGXcalc, case 1
	%-Don't compute => no GMsca (iGMsca==9) or GloNorm (iGloNorm==9)
	g = [];
case 2
	%-User specified globals
case 3
	%-Compute as mean voxel value (within per image fullmean/8 mask)
	g     = zeros(nScan,1 );
	fprintf('%-40s: %30s','Calculating globals',' ')             %-#
	for i = 1:nScan
		str = sprintf('%3d/%-3d',i,nScan);
		fprintf('%s%30s',sprintf('\b')*ones(1,30),str)%-#
		g(i) = spm_global(VY(i));
	end
	fprintf('%s%30s\n',sprintf('\b')*ones(1,30),'...done')       %-#
otherwise
	error('illegal iGXcalc')
end
rg = g;


fprintf('%-40s: ','Design configuration')                            %-#


%-Scaling: compute global scaling factors gSF required to implement proportional
% scaling global normalisation (PropSca) or grand mean scaling (GMsca),
% as specified by iGMsca (& iGloNorm)
%-----------------------------------------------------------------------
switch iGMsca, case 8
	%-Proportional scaling global normalisation
	if iGloNorm~=8, error('iGloNorm-iGMsca(8) mismatch for PropSca'), end
	gSF    = GM./g;
	g      = GM*ones(nScan,1);
case {1,2,3,4,5,6,7}
	%-Grand mean scaling according to iGMsca
	gSF    = GM./spm_meanby(g,eval(CCforms{iGMsca}));
	g      = g.*gSF;
case 9
	%-No grand mean scaling
	gSF    = ones(nScan,1);
otherwise
	error('illegal iGMsca')
end


%-Apply gSF to memory-mapped scalefactors to implement scaling
%-----------------------------------------------------------------------
for i = 1:nScan
	VY(i).pinfo(1:2,:) = VY(i).pinfo(1:2,:)*gSF(i);
end


%-AnCova: Construct global nuisance covariates partition (if AnCova)
%-----------------------------------------------------------------------
if any(iGloNorm == [1:7])

	%-Centre global covariate as requested
	%---------------------------------------------------------------
	switch iGC, case {1,2,3,4,5,6,7}	%-Standard sCC options
		gc = spm_meanby(g,eval(CCforms{iGC}));
	case 8					%-No centering
		gc = 0;
	case 9					%-User specified centre
		%-gc set above
	case 10					%-As implied by AnCova option
		gc = spm_meanby(g,eval(CCforms{iGloNorm}));
	case 11					%-Around GM
		gc = GM;
	otherwise				%-unknown iGC
		error('unexpected iGC value')
	end
	
	
	%-AnCova - add scaled centred global to DesMtx `G' partition
	%---------------------------------------------------------------
	rcname     = 'global'; 
	tI         = [eval(CFIforms{iGloNorm,1}),g - gc];
	tConst     = CFIforms{iGloNorm,2};
	tFnames    = [eval(CFIforms{iGloNorm,3}),{rcname}];
	[f,gnames]  = spm_DesMtx(tI,tConst,tFnames);
	clear tI tConst tFnames
	
	%-Save GX info in xC struct for reference
	%---------------------------------------------------------------
	str     = {sprintf('%s: %s',dstr{2},rcname)};
	if any(iGMsca==[1:7]), str=[str;{['(after ',sGMsca,')']}]; end
	if iGC ~= 8, str=[str;{['used centered ',sCC{iGC}]}]; end
	if iGloNorm > 1
		str=[str;{['fitted as interaction ',sCFI{iGloNorm}]}]; 
	end
	tmp  = struct(	'rc',rg.*gSF,		'rcname',rcname,...
			'c',f,			'cname'	,{gnames},...
			'iCC',iGC,		'iCFI'	,iGloNorm,...
			'type',			3,...
			'cols',[1:size(f,2)] + size([H C B G],2),...
				'descrip',		{str}		);

	G = [G,f]; Gnames = [Gnames; gnames];
	if isempty(xC), xC = tmp; else, xC = [xC,tmp]; end


elseif iGloNorm==8 | iGXcalc>1

	%-Globals calculated, but not AnCova: Make a note of globals
	%---------------------------------------------------------------
	if iGloNorm==8
		str = { 'global values: (used for proportional scaling)';...
			'("raw" unscaled globals shown)'};
	elseif isfinite(M_T) & ~isreal(M_T)
		str = { 'global values: (used to compute analysis threshold)'};
	else
		str = { 'global values: (computed but not used)'};
	end

	rcname ='global';
	tmp     = struct(	'rc',rg,	'rcname',rcname,...
				'c',{[]},	'cname'	,{{}},...
				'iCC',0,	'iCFI'	,0,...
				'type',		3,...
				'cols',		{[]},...
				'descrip',	{str}			);

	if isempty(xC), xC = tmp; else, xC = [xC,tmp]; end
end


%-Save info on global calculation in xGX structure
%-----------------------------------------------------------------------
xGX = struct(...
	'iGXcalc',iGXcalc,	'sGXcalc',sGXcalc,	'rg',rg,...
	'iGMsca',iGMsca,	'sGMsca',sGMsca,	'GM',GM,'gSF',gSF,...
	'iGC',	iGC,		'sGC',	sCC{iGC},	'gc',	gc,...
	'iGloNorm',iGloNorm,	'sGloNorm',sGloNorm);



%-Construct masking information structure and compute actual analysis
% threshold using scaled globals (rg.*gSF)
%-----------------------------------------------------------------------
if isreal(M_T),	M_TH =      M_T  * ones(nScan,1);	%-NB: -Inf is real
else,		M_TH = imag(M_T) * (rg.*gSF); end

if ~isempty(M_P)
	VM  = spm_vol(char(M_P));
	xsM.Explicit_masking = [{'Yes: mask images :'};{VM.fname}'];
else
	VM  = [];
	xsM.Explicit_masking = 'No';
end
xM     = struct('T',M_T, 'TH',M_TH, 'I',M_I, 'VM',{VM}, 'xs',xsM);


%-Construct full design matrix (X), parameter names and structure (xX)
%=======================================================================
X      = [H C B G];
tmp    = cumsum([size(H,2), size(C,2), size(B,2), size(G,2)]);
xX     = struct(	'X',		X,...
			'iH',		[1:size(H,2)],...
			'iC',		[1:size(C,2)] + tmp(1),...
			'iB',		[1:size(B,2)] + tmp(2),...
			'iG',		[1:size(G,2)] + tmp(3),...
			'name',		{[Hnames; Cnames; Bnames; Gnames]},...
			'I',		I,...
			'sF',		{D.sF});


%-Design description (an nx2 cellstr) - for saving and display
%=======================================================================
tmp = {	sprintf('%d condition, +%d covariate, +%d block, +%d nuisance',...
		size(H,2),size(C,2),size(B,2),size(G,2));...
	sprintf('%d total, having %d degrees of freedom',...
		size(X,2),rank(X));...
	sprintf('leaving %d degrees of freedom from %d images',...
		size(X,1)-rank(X),size(X,1))				};
xsDes = struct(	'Design',			{D.DesName},...
		'Global_calculation',		{sGXcalc},...
		'Grand_mean_scaling',		{sGMsca},...
		'Global_normalisation',		{sGloNorm},...
		'Parameters',			{tmp}			);


fprintf('%30s\n','...done')                                          %-#



%-Assemble SPM structure
%=======================================================================
SPM.xY.P	= P;			% filenames
SPM.xY.VY	= VY;			% mapped data
SPM.nscan	= size(xX.X,1);		% scan number
SPM.xX		= xX;			% design structure
SPM.xC		= xC;			% covariate structure
SPM.xGX		= xGX;			% global structure
SPM.xVi		= xVi;			% non-sphericity structure
SPM.xM		= xM;			% mask structure
SPM.xsDes	= xsDes;		% description
SPM.SPMid	= SPMid;		% version

%-Save SPM.mat and set output argument
%-----------------------------------------------------------------------
fprintf('%-40s: ','Saving SPM configuration')   %-#

cd(SPMBATCH.parentdir);
if (exist([anadir])~=7)
    mkdir(anadir);
end

save([anadir,'/SPM'],'SPM');
SPM.swd=[SPMBATCH.parentdir '/' anadir];
delete([SPM.swd,'/mask.img']); %delete 'mask.img' from parent workdir, otherwise SPM will ask if should be overwritten

fprintf('%30s\n','...SPM.mat saved')                                 %-#
varargout = {SPM};

%-Display Design report
%=======================================================================
fprintf('%-40s: ','Design reporting')                                %-#
fname     = cat(1,{SPM.xY.VY.fname}');
spm_DesRep('DesMtx',SPM.xX,fname,SPM.xsDes)
fprintf('%30s\n','...done')     


%-End: Cleanup GUI
%=======================================================================
spm_clf(Finter)
spm('Pointer','Arrow')
fprintf('%-40s: %30s\n','Completed',spm('time'))                     %-#
spm('FigName','Stats: configured',Finter,CmdLine);
spm('Pointer','Arrow')
fprintf('\n\n')



function varargout=get_files_indices(varargin)

%=======================================================================
% - Get files and factor indices
%=======================================================================
% [P,I] = spm_spm_ui('Files&Indices',DsF,Dn,DbaTime,nV)
% DbaTime=D.b.aTime; Dn=D.n; DsF=D.sF;
if nargin<5, nV = 1; else, nV = varargin{5}; end
if nargin<4, DbaTime = 1; else, DbaTime = varargin{4}; end
if nargin<3, Dn  = [Inf,Inf,Inf,Inf]; else, Dn=varargin{3}; end
if nargin<2, DsF = {'Fac1','Fac2','Fac3','Fac4'}; else, DsF=varargin{2}; end
AF=varargin{1};
%-Initialise variables
%-----------------------------------------------------------------------
i4 = [];		% factor 4 index (usually group)
i3 = [];		% factor 3 index (usually subject), per f4
i2 = [];		% factor 2 index (usually condition), per f3/f4
i1 = [];		% factor 1 index (usually replication), per f2/f3/f4
P  = {};		% cell array of string filenames

%-Accrue filenames and factor level indicator vectors
%-----------------------------------------------------------------------
bMV = nV>1;
if isinf(Dn(4)), n4 = spm_input(['#',DsF{4},'''s'],'+1','n1');
	else, n4 = Dn(4); end
bL4 = n4>1;

ti2 = '';
GUIpos = spm_input('!NextPos');
for j4  = 1:n4
    spm_input('!SetNextPos',GUIpos);
    sF4P=''; if bL4, sF4P=[DsF{4},' ',int2str(j4),': ']; end
    if isinf(Dn(3)), n3=spm_input([sF4P,'#',DsF{3},'''s'],'+1','n1');
	    else, n3 = Dn(3); end
    bL3 = n3>1;
    
    if DbaTime & Dn(2)>1
	%disp('NB:selecting in time order - manually specify conditions')
	%-NB: This means f2 levels might not be 1:n2
	GUIpos2 = spm_input('!NextPos');
	for j3 = 1:n3
	    sF3P=''; if bL3, sF3P=[DsF{3},' ',int2str(j3),': ']; end
	    str = [sF4P,sF3P];
	    tP  = {};
	    n21 = Dn(2)*Dn(1);
	    for v=1:nV
	    	vstr=''; if bMV, vstr=sprintf(' (var-%d)',v); end
	    	for tsi=1:length(AF.img_nrs)
                ttP{tsi} =[AF.subj_dirs{j3},filesep,sprintf('%s_%04i.img',AF.img_mask,AF.img_nrs(tsi))]; 
            end
            n21 = length(ttP);
	    	tP  = [tP,ttP];
	    end
	    ti2 = spm_input([str,' ',DsF{2},'?'],GUIpos2,'c',ti2',n21,Dn(2));
	    %-Work out i1 & check
	    [tl2,null,j] = unique(ti2);
	    tn1 = zeros(size(tl2)); ti1 = zeros(size(ti2));
	    for i=1:length(tl2)
		    tn1(i)=sum(j==i); ti1(ti2==tl2(i))=1:tn1(i); end
	    if isfinite(Dn(1)) & any(tn1~=Dn(1))
		%-#i1 levels mismatches specification in Dn(1)
		error(sprintf('#%s not %d as pre-specified',DsF{1},Dn(1)))
	    end
	    P  = [P;tP];
	    i4 = [i4; j4*ones(n21,1)];
	    i3 = [i3; j3*ones(n21,1)];
	    i2 = [i2; ti2];
	    i1 = [i1; ti1];
	end

    else

	if isinf(Dn(2))
	    n2 = spm_input([sF4P,'#',DsF{2},'''s'],'+1','n1');
	else
	    n2 = Dn(2);
	end
	bL2 = n2>1;

	if n2==1 & Dn(1)==1 %-single scan per f3 (subj)
	    %disp('NB:single scan per f3')
	    str = [sF4P,'select images, ',DsF{3},' 1-',int2str(n3)];
	    tP = {};
	    for v=1:nV
	    	vstr=''; if bMV, vstr=sprintf(' (var-%d)',v); end
            for j3=1:n3,
  	    	    for tsi=1:length(AF.img_nrs)
                    ttP{tsi,1} =[AF.subj_dirs{j3},filesep,sprintf('%s_%04i.img',AF.img_mask,AF.img_nrs(tsi))]; 
                end
      	    	%ttP = spm_get(n3,'.img',{[str,vstr]});
	    	    tP = [tP,ttP];

            end

	    end
	    P   = [P;tP];
	    i4  = [i4; j4*ones(n3,1)];
	    i3  = [i3; [1:n3]'];
	    i2  = [i2; ones(n3,1)];
	    i1  = [i1; ones(n3,1)];
	else
	    %-multi scan per f3 (subj) case
	    %disp('NB:multi scan per f3')
	    for j3 = 1:n3
		sF3P=''; if bL3, sF3P=[DsF{3},' ',int2str(j3),': ']; end
		if Dn(1)==1
			%-No f1 (repl) within f2 (cond)
			%disp('NB:no f1 within f2')
			str = [sF4P,sF3P,'select images: ',DsF{2},...
				 ' 1-',int2str(n2)];
			tP = {};
			for v=1:nV
       	    	for tsi=1:length(AF.img_nrs)
                    ttP{tsi,1} =[AF.subj_dirs{j3},filesep,sprintf('%s_%04i.img',AF.img_mask,AF.img_nrs(tsi))]; 
                end
				tP = [tP,ttP];
			end
			P   = [P;tP];
			i4  = [i4; j4*ones(n2,1)];
			i3  = [i3; j3*ones(n2,1)];
			i2  = [i2; [1:n2]'];
			i1  = [i1; ones(n2,1)];
		else
		    %-multi f1 (repl) within f2 (cond)
		    %disp('NB:f1 within f2')
		    for j2  = 1:n2
			sF2P='';
			if bL2, sF2P=[DsF{2},' ',int2str(j2),': ']; end
			str = [sF4P,sF3P,sF2P,' select images...'];
			tP  = {};
			n1  = Dn(1);
			for v=1:nV
       	    	for tsi=1:length(AF.img_nrs)
                    ttP{tsi} =[AF.subj_dirs{j3},filesep,sprintf('%s_%04i.img',AF.img_mask,AF.img_nrs(tsi))]; 
                end
				n1  = length(ttP);
				tP  = [tP,ttP];
			end
			P   = [P;tP];
			i4  = [i4; j4*ones(n1,1)];
			i3  = [i3; j3*ones(n1,1)];
			i2  = [i2; j2*ones(n1,1)];
			i1  = [i1; [1:n1]'];
		    end                         % (for j2)
		end                             % (if Dn(1)==1)
	    end                                 % (for j3)
    end                                     % (if  n2==1 &...)
    end                                         % (if DbaTime & Dn(2)>1)
end                                             % (for j4)
varargout = {P,[i1,i2,i3,i4]};


function str = sf_estrrep(str,srstr)
%=======================================================================
for i = 1:size(srstr,1)
	str = strrep(str,srstr{i,1},srstr{i,2});
end