%-Generic factor names
my_sF = {'sF1','sF2','sF3','sF4'};

%-Covariate by factor interaction options
sCFI = {'<none>';...							%-1
	'with sF1';'with sF2';'with sF3';'with sF4';...			%-2:5
	'with sF2 (within sF4)';'with sF3 (within sF4)'};		%-6,7

%-DesMtx argument components for covariate by factor interaction options
% (Used for CFI's Covariate Centering (CC), GMscale & Global normalisation)
CFIforms = {	'[]',		'C',	'{}';...			%-1
		'I(:,1)',	'FxC',	'{D.sF{1}}';...			%-2
		'I(:,2)',	'FxC',	'{D.sF{2}}';...			%-3
		'I(:,3)',	'FxC',	'{D.sF{3}}';...			%-4
		'I(:,4)',	'FxC',	'{D.sF{4}}';...			%-5
		'I(:,[4,2])',	'FxC',	'{D.sF{4},D.sF{2}}';...		%-6
		'I(:,[4,3])',	'FxC',	'{D.sF{4},D.sF{3}}'	};	%-7

%-Centre (mean correction) options for covariates & globals            (CC)
% (options 9-12 are for centering of global when using AnCova GloNorm) (GC)
sCC = {		'around overall mean';...				%-1
		'around sF1 means';...					%-2
		'around sF2 means';...					%-3
		'around sF3 means';...					%-4
		'around sF4 means';...					%-5
		'around sF2 (within sF4) means';...			%-6
		'around sF3 (within sF4) means';...			%-7
		'<no centering>';...					%-8
		'around user specified value';...			%-9
		'(as implied by AnCova)';...				%-10
		'GM';...						%-11
		'(redundant: not doing AnCova)'}';			%-12
%-DesMtx I forms for covariate centering options
CCforms = {'ones(nScan,1)',CFIforms{2:end,1},''}';


%-Global normalization options (options 1-7 match CFIforms)       (GloNorm)
sGloNorm = {	'AnCova';...						%-1
		'AnCova by sF1';...					%-2
		'AnCova by sF2';...					%-3
		'AnCova by sF3';...					%-4
		'AnCova by sF4';...					%-5
		'AnCova by sF2 (within sF4)';...			%-6
		'AnCova by sF3 (within sF4)';...			%-7
		'proportional scaling';...				%-8
		'<no global normalisation>'};				%-9

%-Grand mean scaling options                                        (GMsca)
sGMsca = {	'scaling of overall grand mean';...			%-1
		'scaling of sF1 grand means';...			%-2
		'scaling of sF2 grand means';...			%-3
		'scaling of sF3 grand means';...			%-4
		'scaling of sF4 grand means';...			%-5
		'scaling of sF2 (within sF4) grand means';...		%-6
		'scaling of sF3 (within sF4) grand means';...		%-7
		'(implicit in PropSca global normalisation)';...	%-8
		'<no grand Mean scaling>'	};			%-9
%-NB: Grand mean scaling by subject is redundent for proportional scaling


%-Global calculation options                                       (GXcalc)
sGXcalc  = {	'omit';...						%-1
		'user specified';...					%-2
		'mean voxel value (within per image fullmean/8 mask)'};	%-3

