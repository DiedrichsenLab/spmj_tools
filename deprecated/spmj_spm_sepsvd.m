function [SPM] = spmj_spm_sepsvd(SPM)
% [Re]ML Estimation of a General Linear Model
% FORMAT [SPM] = spmj_spm_sepsvd(SPM)
% Estimation of finite response design, assuming seperability of time
% course and response strength.
% Allows for a seperate hrf time course for each voxel
%
% Required fields of SPM:
%
% xY.VY - nScan x 1 struct array of image handles (see spm_vol)
%         Images must have the same orientation, voxel size and data type
%       - Any scaling should have already been applied via the image handle
%         scalefactors.
%
% xX    - Structure containing design matrix information
%       - Required fields are:
%         xX.X      - Design matrix (raw, not temporally smoothed)
%         xX.name   - cellstr of parameter names corresponding to columns
%                     of design matrix
%       - Optional fields are:
%         xX.K      - cell of session-specific structures (see spm_filter)
%                   - Design & data are pre-multiplied by K
%                     (K*Y = K*X*beta + K*e)
%                   - Note that K should not smooth across block boundaries
%                   - defaults to speye(size(xX.X,1))
%         xX.W      - Optional whitening/weighting matrix used to give
%                     weighted least squares estimates (WLS). If not specified
%                     spm_spm will set this to whiten the data and render
%                     the OLS estimates maximum likelihood
%                     i.e. W*W' = inv(xVi.V).
%
% xVi   - Structure describing intrinsic temporal non-sphericity
%       - Required fields are:
%         xVi.Vi    - array of non-sphericity components
%                   - defaults to {speye(size(xX.X,1))} - i.i.d.
%                   - specifying a cell array of constraints (Qi)
%                     These constraints invoke spm_reml to estimate
%                     hyperparameters assuming V is constant over voxels.
%                     that provide a high precise estimate of xX.V
%       - Optional fields are:
%         xX.V      - Optional non-sphericity matrix.  Cov(e) = sigma^2*V
%                     If not specified spm_spm will compute this using
%                     a 1st pass to identify significant voxels over which
%                     to estimate V.  A 2nd pass is then used to re-estimate
%                     the parameters with WLS and save the ML estimates
%                     (unless xX.W is already specified).
%
% xM    - Structure containing masking information, or a simple column vector
%         of thresholds corresponding to the images in VY [default: -Inf]
%       - If a structure, the required fields are:
%         xM.TH - nVar x nScan matrix of analysis thresholds, one per image
%         xM.I  - Implicit masking (0=>none, 1 => implicit zero/NaN mask)
%         xM.VM - struct array of explicit mask image handles
%       - (empty if no explicit masks)
%               - Explicit mask images are >0 for valid voxels to assess.
%               - Mask images can have any orientation, voxel size or data
%                 type. They are interpolated using nearest neighbour
%                 interpolation to the voxel locations of the data Y.
%       - Note that voxels with constant data (i.e. the same value across
%         scans) are also automatically masked out.
%
% swd   - Directory where the output files will be saved [default: pwd]
%         If exists, it becomes the current working directory.
%
% In addition, global SPM "defaults" variable is used (see spm_defaults):
%
% stats.<modality>.UFp - critical F-threshold for selecting voxels over
%                        which the non-sphericity is estimated (if
%                        required) [default: 0.001]
%
% stats.maxres         - maximum number of residual images for smoothness
%                        estimation
%
% stats.maxmem         - maximum amount of data processed at a time (in bytes)
%
% modality             - SPM modality {'PET','FMRI','EEG'}
%
%__________________________________________________________________________
%
%
%                           ----------------
%
% The following SPM.fields are set by spm_spm (unless specified)
%
%     xVi.V      - estimated non-sphericity trace(V) = rank(V)
%     xVi.h      - hyperparameters  xVi.V = xVi.h(1)*xVi.Vi{1} + ...
%     xVi.Cy     - spatially whitened <Y*Y'> (used by ReML to estimate h)
%     xVi.CY     - <(Y - <Y>)*(Y - <Y>)'>    (used by spm_spm_Bayes)
%
%                           ----------------
%
%     Vbeta     - struct array of beta image handles (relative)
%     VResMS    - file struct of ResMS image handle  (relative)
%     VM        - file struct of Mask  image handle  (relative)
%
%                           ----------------
%
%     xX.W      - if not specified W*W' = inv(x.Vi.V)
%     xX.V      - V matrix (K*W*Vi*W'*K') = correlations after K*W is applied
%     xX.xKXs   - space structure for K*W*X, the 'filtered and whitened'
%                 design matrix
%               - given as spm_sp('Set',xX.K*xX.W*xX.X) - see spm_sp
%     xX.pKX    - pseudoinverse of K*W*X, computed by spm_sp
%     xX.Bcov   - xX.pKX*xX.V*xX.pKX - variance-covariance matrix of
%                 parameter estimates
%                 (when multiplied by the voxel-specific hyperparameter ResMS
%                 of the parameter estimates (ResSS/xX.trRV = ResMS) )
%     xX.trRV   - trace of R*V
%     xX.trRVRV - trace of RVRV
%     xX.erdf   - effective residual degrees of freedom (trRV^2/trRVRV)
%     xX.nKX    - design matrix (xX.xKXs.X) scaled for display
%                 (see spm_DesMtx('sca',... for details)
%
%                           ----------------
%
%     xVol.M    - 4x4 voxel->mm transformation matrix
%     xVol.iM   - 4x4 mm->voxel transformation matrix
%     xVol.DIM  - image dimensions - column vector (in voxels)
%     xVol.XYZ  - 3 x S vector of in-mask voxel coordinates
%     xVol.S    - Lebesgue measure or volume       (in voxels)
%     xVol.R    - vector of resel counts           (in resels)
%     xVol.FWHM - Smoothness of components - FWHM, (in voxels)
%
%                           ----------------
%
%     xCon      - Contrast structure (created by spm_FcUtil.m)
%     xCon.name - Name of contrast
%     xCon.STAT - 'F', 'T' or 'P' - for F/T-contrast ('P' for PPMs)
%     xCon.c    - (F) Contrast weights
%     xCon.X0   - Reduced design matrix (spans design space under Ho)
%                 It is in the form of a matrix (spm99b) or the
%                 coordinates of this matrix in the orthogonal basis
%                 of xX.X defined in spm_sp.
%     xCon.iX0  - Indicates how contrast was specified:
%                 If by columns for reduced design matrix then iX0 contains
%                 the column indices. Otherwise, it's a string containing
%                 the spm_FcUtil 'Set' action: Usually one of {'c','c+','X0'}
%                 (Usually this is the input argument F_iX0.)
%     xCon.X1o  - Remaining design space (orthogonal to X0).
%                 It is in the form of a matrix (spm99b) or the
%                 coordinates of this matrix in the orthogonal basis
%                 of xX.X defined in spm_sp.
%     xCon.eidf - Effective interest degrees of freedom (numerator df)
%     xCon.Vcon - ...for handle of contrast/ESS image (empty at this stage)
%     xCon.Vspm - ...for handle of SPM image (empty at this stage)
%
%                           ----------------
%
%
% The following images are written to file
%
% mask.{img,hdr}                                   - analysis mask image
% 8-bit (uint8) image of zero-s & one's indicating which voxels were
% included in the analysis. This mask image is the intersection of the
% explicit, implicit and threshold masks specified in the xM argument.
% The XYZ matrix contains the voxel coordinates of all voxels in the
% analysis mask. The mask image is included for reference, but is not
% explicitly used by the results section.
%
%                           ----------------
%
% beta_????.{img,hdr}                                 - parameter images
% These are 32-bit (float32) images of the parameter estimates. The image
% files are numbered according to the corresponding column of the
% design matrix. Voxels outside the analysis mask (mask.img) are given
% value NaN.
%
%                           ----------------
% hrf_????.{img,hdr}
% Shows the time course of each voxel: As many as FIR in the original
% design matrix
%
% ResMS.{img,hdr}                    - estimated residual variance image
% This is a 64-bit (float64) image of the residual variance estimate.
% Voxels outside the analysis mask are given value NaN.
%
%                           ----------------
%
% ResI_????.{img,hdr}        - standardised residual (temporary) images
% These are 64-bit (float64) images of standardised residuals. At most
% maxres images will be saved and used by spm_est_smoothness, after which
% they will be deleted.
%
%--------------------------------------------------------------------------
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% (c) Joern Diedrichsen

SVNid     = '$Rev: 3468 $';

%-Say hello
%--------------------------------------------------------------------------
SPMid     = spm('FnBanner',mfilename,SVNid);
Finter    = spm('FigName','Stats: estimation...'); spm('Pointer','Watch');

%-Get SPM.mat[s] if necessary
%--------------------------------------------------------------------------
if nargin == 0
    P     = cellstr(spm_select(Inf,'^SPM\.mat$','Select SPM.mat[s]'));
    for i = 1:length(P)
        swd     = fileparts(P{i});
        load(fullfile(swd,'SPM.mat'));
        SPM.swd = swd;
        spmj_spm_sepsvd(SPM);
    end
    return
end

%-Change to SPM.swd if specified
%--------------------------------------------------------------------------
try
    cd(SPM.swd);
catch
    SPM.swd = pwd;
end

%-Ensure data are assigned
%--------------------------------------------------------------------------
try
    SPM.xY.VY;
catch
    spm('alert!','Please assign data to this design', mfilename);
    spm('FigName','Stats: done',Finter); spm('Pointer','Arrow')
    return
end

%-Delete files from previous analyses
%--------------------------------------------------------------------------
if exist(fullfile(SPM.swd,'mask.img'),'file') == 2
    
    str = {'Current directory contains SPM estimation files:',...
        'pwd = ',SPM.swd,...
        'Existing results will be overwritten!'};
    if spm_input(str,1,'bd','stop|continue',[1,0],1)
        spm('FigName','Stats: done',Finter); spm('Pointer','Arrow')
        return
    else
        warning('Overwriting old results\n\t (pwd = %s) ',SPM.swd);
        try, SPM = rmfield(SPM,'xVol'); end
    end
end

files = {'^mask\..{3}$','^ResMS\..{3}$','^RPV\..{3}$',...
    '^beta_.{4}\..{3}$','^hrf_.{4}\..{3}$','^con_.{4}\..{3}$','^ResI_.{4}\..{3}$',...
    '^ess_.{4}\..{3}$', '^spm\w{1}_.{4}\..{3}$'};

for i=1:length(files)
    j = spm_select('List',SPM.swd,files{i});
    for k=1:size(j,1)
        spm_unlink(deblank(j(k,:)));
    end
end


%==========================================================================
% - A N A L Y S I S   P R E L I M I N A R I E S
%==========================================================================

%-Initialise
%==========================================================================
fprintf('%-40s: %30s','Initialising parameters','...computing');        %-#
xX            = SPM.xX;
[nScan nBetaInitial] = size(xX.X);
nRun=length(SPM.nscan);


%-If xM is not a structure then assume it's a vector of thresholds
%--------------------------------------------------------------------------
try
    xM = SPM.xM;
catch
    xM = -Inf(nScan,1);
end
if ~isstruct(xM)
    xM = struct('T',    [],...
        'TH',   xM,...
        'I',    0,...
        'VM',   {[]},...
        'xs',   struct('Masking','analysis threshold'));
end

%-Check confounds (xX.K) and non-sphericity (xVi)
%--------------------------------------------------------------------------
if ~isfield(xX,'K')
    xX.K  = 1;
end
try
    %-If covariance components are specified use them
    %----------------------------------------------------------------------
    xVi   = SPM.xVi;
catch
    
    %-otherwise assume i.i.d.
    %----------------------------------------------------------------------
    xVi   = struct( 'form',  'i.i.d.',...
        'V',     speye(nScan,nScan));
end


%-Get non-sphericity V
%==========================================================================
try
    %-If xVi.V is specified proceed directly to parameter estimation
    %----------------------------------------------------------------------
    V     = xVi.V;
catch
    % otherwise invoke ReML selecting voxels under i.i.d assumptions
    %----------------------------------------------------------------------
    V     = speye(nScan,nScan);
    xVi=V;      % Just do OLS
end
str   = 'parameter estimation';

%-Get whitening/Weighting matrix: If xX.W exists we will save WLS estimates
%--------------------------------------------------------------------------
try
    %-If W is specified, use it
    %----------------------------------------------------------------------
    W     = xX.W;
catch
    % otherwise make W a whitening filter W*W' = inv(V)
    %------------------------------------------------------------------
    [u s] = spm_svd(xVi.V);
    s     = spdiags(1./sqrt(diag(s)),0,length(s),length(s));
    W     = u*s*u';
    W     = W.*(abs(W) > 1e-6);
    xX.W  = sparse(W);
end


%-Design space and projector matrix [pseudoinverse] for WLS
%==========================================================================
xX.xKXs   = spm_sp('Set',spm_filter(xX.K,W*xX.X));       % KWX
xX.xKXs.X = full(xX.xKXs.X);
xX.pKX    = spm_sp('x-',xX.xKXs);                        % projector
erdf      = spm_SpUtil('trRV',xX.xKXs);                  % Working error df


%-Check Design and basis function
%==========================================================================
if (~strcmp(SPM.xBF.name,'Finite Impulse Response'))
    error('spmj_spm_sepsvd: original SPM must be a FIR design');
end;
nFIR=SPM.xBF.order;
for s=1:nRun
    xX.iCs=zeros(1,nBetaInitial);   % number of the session
    xX.iCc=zeros(1,nBetaInitial);   % number of the condition
    xX.iCb=zeros(1,nBetaInitial);  % number of the basis function
    xX.iNs=zeros(1,nBetaInitial);   % Flag for the regressors of no interest (+ intercept)
end;
% Figure out the indices
for i=1:length(xX.name)
    a=textscan(xX.name{i},'Sn(%d) %s');
    xX.iCs(i)=a{1};
    name=a{2}{1};
    if (name(1)=='C')
        a=textscan(name,'C%d*bf(%d)');
        xX.iCc(i)=a{1};
        xX.iCb(i)=a{2};
    elseif (name(1)=='N')
        xX.iNs(i)=1;
    elseif (strcmp(name,'constant'))
        xX.iNs(i)=2;
    else
        error('Unknown regressor string');
    end;
end;

%-Prepare convolution matrix + final result SPM
%==========================================================================
save SPMfir.mat SPM         % Saftey copy of the finite response function SPM 
i=find(xX.iCb==1 | xX.iNs>0);
Xconv.X=xX.X(:,i);
nBeta=length(i);
Xconv.iCc=xX.iCc(i);
Xconv.iNs=xX.iNs(i);
Xconv.iCs=xX.iCs(i);
Xconv.Xv=Xconv.X; % This is the copy that will be changed on a voxel-by-voxel basis
Xconv.iC=find(Xconv.iCc>0 | Xconv.iNs==1);
Xconv.iB=find(Xconv.iNs==2);                 % constant terms
Xconv.name={xX.name{i}};

%-Image dimensions and data
%==========================================================================
VY       = SPM.xY.VY;
spm_check_orientations(VY);

for i = 1:numel(VY)
    % check files exists and try pwd
    %----------------------------------------------------------------------
    if ~spm_existfile(VY(i).fname)
        [p,n,e]     = fileparts(VY(i).fname);
        VY(i).fname = [n,e];
    end
end

M        = VY(1).mat;
DIM      = VY(1).dim(1:3)';
xdim     = DIM(1); ydim = DIM(2); zdim = DIM(3);
YNaNrep  = spm_type(VY(1).dt(1),'nanrep');


fprintf('%s%30s\n',repmat(sprintf('\b'),1,30),'...done');               %-#

%-Initialise output images
%==========================================================================
fprintf('%-40s: %30s','Output images','...initialising');           %-#

%-Initialise new mask name: current mask & conditions on voxels
%----------------------------------------------------------------------
VM    = struct('fname',  'mask.img',...
    'dim',    DIM',...
    'dt',     [spm_type('uint8') spm_platform('bigend')],...
    'mat',    M,...
    'pinfo',  [1 0 0]',...
    'descrip','spm_spm:resultant analysis mask');
VM    = spm_create_vol(VM);


%-Initialise beta image files
%----------------------------------------------------------------------
Vbeta(1:nBeta) = deal(struct(...
    'fname',    [],...
    'dim',      DIM',...
    'dt',       [spm_type('float32') spm_platform('bigend')],...
    'mat',      M,...
    'pinfo',    [1 0 0]',...
    'descrip',  ''));
for i = 1:nBeta
    Vbeta(i).fname   = sprintf('beta_%04d.img',i);
    Vbeta(i).descrip = sprintf('spm_spm:beta (%04d) - %s',i,xX.name{i});
end
Vbeta = spm_create_vol(Vbeta);


%-Initialise residual sum of squares image file
%----------------------------------------------------------------------
VResMS = struct('fname',    'ResMS.img',...
    'dim',      DIM',...
    'dt',       [spm_type('float64') spm_platform('bigend')],...
    'mat',      M,...
    'pinfo',    [1 0 0]',...
    'descrip',  'spm_spm:Residual sum-of-squares');
VResMS = spm_create_vol(VResMS);


%-Initialise hrf images
%----------------------------------------------------------------------
Vhrf(1:nFIR) = deal(struct(...
    'fname',    [],...
    'dim',      DIM',...
    'dt',       [spm_type('float32') spm_platform('bigend')],...
    'mat',      M,...
    'pinfo',    [1 0 0]',...
    'descrip',  'spm_spm:StandardisedResiduals'));

for i = 1:nFIR
    Vhrf(i).fname   = sprintf('hrf_%04d.img', i);
    Vhrf(i).descrip = sprintf('spm_spm_sepsvd:hrf(%04d)', i);
end
Vhrf = spm_create_vol(Vhrf);
fprintf('%s%30s\n',repmat(sprintf('\b'),1,30),'...initialised');    %-#


%==========================================================================
% - F I T   M O D E L   &   W R I T E   P A R A M E T E R    I M A G E S
%==========================================================================

%-MAXMEM is the maximum amount of data processed at a time (bytes)
%--------------------------------------------------------------------------
MAXMEM = spm_get_defaults('stats.maxmem');
mmv    = MAXMEM/8/nScan;
blksz  = min(xdim*ydim,ceil(mmv));                             %-block size
nbch   = ceil(xdim*ydim/blksz);                                %-# blocks
nbz    = max(1,min(zdim,floor(mmv/(xdim*ydim))));   nbz = 1;   %-# planes
blksz  = blksz * nbz;

%-Initialise variables used in the loop
%==========================================================================
[xords, yords] = ndgrid(1:xdim, 1:ydim);
xords = xords(:)'; yords = yords(:)';           % plane X,Y coordinates
S     = 0;                                      % Volume (voxels)
% s     = 0;                                      % Volume (voxels > UF)
% Cy    = 0;                                      % <Y*Y'> spatially whitened
% CY    = 0;                                      % <(Y - <Y>) * (Y - <Y>)'>
% EY    = 0;                                      % <Y>    for ReML

%-Initialise XYZ matrix of in-mask voxel co-ordinates (real space)
%--------------------------------------------------------------------------
XYZ   = zeros(3,xdim*ydim*zdim);

%-Cycle over bunches blocks within planes to avoid memory problems
%==========================================================================
spm_progress_bar('Init',100,str,'');

for z = 1:nbz:zdim                       %-loop over planes (2D or 3D data)
    
    % current plane-specific parameters
    %----------------------------------------------------------------------
    CrPl    = z:min(z+nbz-1,zdim);       %-plane list
    zords   = CrPl(:)*ones(1,xdim*ydim); %-plane Z coordinates
    CrBl    = [];                        %-parameter estimates
    CrResI  = [];                        %-normalized residuals
    CrResSS = [];                        %-residual sum of squares
    CrHrf   = [];                        %-voxel-based hrf estimates
    Q       = [];                        %-in mask indices for this plane
    
    for bch = 1:nbch                     %-loop over blocks
        
        %-Print progress information in command window
        %------------------------------------------------------------------
        if numel(CrPl) == 1
            str = sprintf('Plane %3d/%-3d, block %3d/%-3d',...
                z,zdim,bch,nbch);
        else
            str = sprintf('Planes %3d-%-3d/%-3d',z,CrPl(end),zdim);
        end
        if z==1&&bch==1, str2=''; else str2=repmat(sprintf('\b'),1,72); end
        fprintf('%s%-40s: %30s',str2,str,' ');                          %-#
        
        %-construct list of voxels in this block
        %------------------------------------------------------------------
        I     = (1:blksz) + (bch - 1)*blksz;       %-voxel indices
        I     = I(I <= numel(CrPl)*xdim*ydim);     %-truncate
        xyz   = [repmat(xords,1,numel(CrPl)); ...
            repmat(yords,1,numel(CrPl)); ...
            reshape(zords',1,[])];
        xyz   = xyz(:,I);                          %-voxel coordinates
        nVox  = size(xyz,2);                       %-number of voxels
        
        %-Get data & construct analysis mask
        %=================================================================
        fprintf('%s%30s',repmat(sprintf('\b'),1,30),'...read & mask data')
        Cm    = true(1,nVox);                      %-current mask
        
        
        %-Compute explicit mask
        % (note that these may not have same orientations)
        %------------------------------------------------------------------
        for i = 1:length(xM.VM)
            
            %-Coordinates in mask image
            %--------------------------------------------------------------
            j = xM.VM(i).mat\M*[xyz;ones(1,nVox)];
            
            %-Load mask image within current mask & update mask
            %--------------------------------------------------------------
            Cm(Cm) = spm_get_data(xM.VM(i),j(:,Cm),false) > 0;
        end
        
        %-Get the data in mask, compute threshold & implicit masks
        %------------------------------------------------------------------
        Y     = zeros(nScan,nVox);
        for i = 1:nScan
            
            %-Load data in mask
            %--------------------------------------------------------------
            if ~any(Cm), break, end                %-Break if empty mask
            Y(i,Cm)  = spm_get_data(VY(i),xyz(:,Cm),false);
            
            Cm(Cm)   = Y(i,Cm) > xM.TH(i);         %-Threshold (& NaN) mask
            if xM.I && ~YNaNrep && xM.TH(i) < 0    %-Use implicit mask
                Cm(Cm) = abs(Y(i,Cm)) > eps;
            end
        end
        
        %-Mask out voxels where data is constant
        %------------------------------------------------------------------
        Cm(Cm) = any(diff(Y(:,Cm),1));
        Y      = Y(:,Cm);                          %-Data within mask
        CrS    = sum(Cm);                          %-# current voxels
        
        
        %==================================================================
        %-Proceed with General Linear Model (if there are voxels)
        %==================================================================
        if CrS
            
            %-Whiten/Weight data and remove filter confounds
            %--------------------------------------------------------------
            fprintf('%s%30s',repmat(sprintf('\b'),1,30),'...filtering');%-#
            
            KWY   = spm_filter(xX.K,W*Y);
            
            %-General linear model: Weighted least squares estimation
            %--------------------------------------------------------------
            fprintf('%s%30s',repmat(sprintf('\b'),1,30),'...estimation');%-#
            
            beta1  = xX.pKX*KWY;                    %-Parameter estimates
            
            % determine the hrf for each voxel
            N=size(KWY,2);
            beta=zeros(size(Xconv.X,2),N);
            HRF=zeros(nFIR,N);
            ResSS=zeros(1,N); 
            CrossSS=zeros(1,N); 
            for n=1:N
                % Caluclate the crossvalidated prediction accuracy 
%                 Xcr.Xv=Xconv.X; 
%                 for s=1:nRun
%                     testB=(xX.iCs==s); 
%                     trainB=(xX.iCs~=s); 
%                     testY=SPM.Sess(s).row; 
%                     trainY=ones(1,nScan);
%                     trainY(testY)=0; 
%                     testC=(Xconv.iCs==s); 
%                     trainC=(Xconv.iCs~=s); 
%                     % 
%                     B=reshape(beta1(xX.iCb>0 & trainB,n),nFIR,sum(xX.iCb(trainB)==1));
%                     [v,a]=svds(B,1);
%                     hrf=sign(mean(v(2:5)))*v;
%                     i1=Xconv.iCc>0;
%                     a=conv2(Xconv.X(:,i1),hrf);
%                     Xcr.Xv(:,i1)=a(1:size(Xconv.X,1),:);
%                     Xcr.xKXs = spm_sp('Set',Xcr.Xv);
%                     Xcr.pKX = spm_sp('x-',Xcr.xKXs);        % Pseudo inverse
%                     bet=Xcr.pKX*KWY(:,n);
%                     yhatN=Xcr.Xv(testY,testC & Xconv.iNs>0)*bet(testC & Xconv.iNs>0,1); 
%                     yhatC=Xcr.Xv(testY,testC & Xconv.iCc>0)*bet(testC & Xconv.iCc>0,1);
%                     
%                     for i=1:max(Xconv.iCc)
%                         mbeta(i,1)=mean(bet(trainC & Xconv.iCc==i)); 
%                     end; 
%                     ycrossC=Xcr.Xv(testY,testC & Xconv.iCc>0)*mbeta; 
%             
%                     t=[1:123];
%                     subplot(2,1,1); 
%                     plot(t,KWY(testY,n),'k',t,yhatN,'r');
%                     subplot(2,1,2); 
%                     plot(t,KWY(testY,n)-yhatN,'k',t,yhatC,'r',t,ycrossC,'g');
%                     keyboard; 
%                     
%                 end; 
                
                % Now get the full set 
                % 
                B=reshape(beta1(xX.iCb>0,n),nFIR,sum(xX.iCb==1));
                [v,a]=svds(B,1);
                HRF(:,n)=sign(mean(v(2:5)))*v;
                i1=Xconv.iCc>0;
                a=conv2(Xconv.X(:,i1),HRF(:,n));
                Xconv.Xv(:,i1)=a(1:size(Xconv.X,1),:);
                Xconv.xKXs = spm_sp('Set',Xconv.Xv);
                Xconv.pKX = spm_sp('x-',Xconv.xKXs);        % Pseudo inverse
                beta(:,n)=Xconv.pKX*KWY(:,n);
                res   = spm_sp('r',Xconv.xKXs,KWY(:,n));       %-Residuals
                ResSS(1,n) = sum(res.^2);                   %-Residual SSQ
            end;
            clear KWY                              %-Clear to save memory
            
            %-sample covariance and mean of Y (all voxels)
            %----------------------------------------------------------
            % CY         = CY + Y*Y';
            % EY         = EY + sum(Y,2);
            
            %-Save betas etc. for current plane as we go along
            %----------------------------------------------------------
            CrBl       = [CrBl,    beta];
            CrResSS    = [CrResSS, ResSS];
            CrHrf      = [CrHrf, HRF];
            clear Y                         %-Clear to save memory
            
        end % (CrS)
        
        %-Append new inmask voxel locations and volumes
        %------------------------------------------------------------------
        XYZ(:,S + (1:CrS)) = xyz(:,Cm);     %-InMask XYZ voxel coords
        Q                  = [Q I(Cm)];     %-InMask XYZ voxel indices
        S                  = S + CrS;       %-Volume analysed (voxels)
        
    end % (bch)
    
    
    %-Plane complete, write plane to image files (unless 1st pass)
    %======================================================================
    
    fprintf('%s%30s',repmat(sprintf('\b'),1,30),'...saving plane'); %-#
    
    jj = NaN(xdim,ydim,numel(CrPl));
    
    %-Write Mask image
    %------------------------------------------------------------------
    if ~isempty(Q), jj(Q) = 1; end
    VM    = spm_write_plane(VM, ~isnan(jj), CrPl);
    
    %-Write beta images
    %------------------------------------------------------------------
    for i = 1:nBeta
        if ~isempty(Q), jj(Q) = CrBl(i,:);
        end
        Vbeta(i) = spm_write_plane(Vbeta(i), jj, CrPl);
    end
    
    %-Write Hrf images
    %------------------------------------------------------------------
    for i = 1:nFIR
        if ~isempty(Q), jj(Q) = CrHrf(i,:);
        end
        Vhrf(i) = spm_write_plane(Vhrf(i), jj, CrPl);
    end
    
    
    %-Write ResSS into ResMS (variance) image scaled by tr(RV) above
    %------------------------------------------------------------------
    if ~isempty(Q), jj(Q) = CrResSS; end
    VResMS  = spm_write_plane(VResMS, jj, CrPl);
    
    
    %-Report progress
    %----------------------------------------------------------------------
    fprintf('%s%30s',repmat(sprintf('\b'),1,30),'...done');             %-#
    spm_progress_bar('Set',100*(bch + nbch*(z - 1))/(nbch*zdim));
    
    
end % (for z = 1:zdim)
fprintf('\n');                                                          %-#
spm_progress_bar('Clear')

%==========================================================================
% - P O S T   E S T I M A T I O N   C L E A N U P
%==========================================================================
if S == 0, spm('alert!','No inmask voxels - empty analysis!'); return; end

%-average sample covariance and mean of Y (over voxels)
%--------------------------------------------------------------------------
% CY   = CY/S;
% EY   = EY/S;
% CY   = CY - EY*EY';


%-Use non-sphericity xVi.V to compute [effective] degrees of freedom
%==========================================================================
Xconv.V         = spm_filter(xX.K,spm_filter(xX.K,W*V*W')');% KWVW'K'
[trRV trRVRV]   = spm_SpUtil('trRV',Xconv.xKXs,Xconv.V);          % trRV (for X)
Xconv.trRV         = trRV;                                     % <R'*y'*y*R>
Xconv.trRVRV       = trRVRV;                                   %-Satterthwaite
Xconv.erdf         = trRV^2/trRVRV;                            % approximation
Xconv.Bcov         = Xconv.pKX*Xconv.V*Xconv.pKX';                      % Cov(beta)

%-Set VResMS scalefactor as 1/trRV (raw voxel data is ResSS)
%--------------------------------------------------------------------------
VResMS.pinfo(1) = 1/Xconv.trRV;
VResMS          = spm_create_vol(VResMS);


%-Compute scaled design matrix for display purposes
%--------------------------------------------------------------------------
Xconv.nKX        = spm_DesMtx('sca',Xconv.X,Xconv.name);


%-Save remaining results files and analysis parameters
%==========================================================================
fprintf('%-40s: %30s','Saving results','...writing');                   %-#

%-place fields in SPM
%--------------------------------------------------------------------------
SPM.xVol.XYZ   = XYZ(:,1:S);        %-InMask XYZ coords (voxels)
SPM.xVol.M     = M;                 %-voxels -> mm
SPM.xVol.iM    = inv(M);            %-mm -> voxels
SPM.xVol.DIM   = DIM;               %-image dimensions
SPM.xVol.FWHM  = NaN;              %-Smoothness data
SPM.xVol.R     = NaN;                 %-Resel counts
SPM.xVol.S     = S;                 %-Volume (voxels)
SPM.xVol.VRpv  = NaN;              %-Filehandle - Resels per voxel

SPM.Vbeta      = Vbeta;             %-Filehandle - Beta
SPM.VResMS     = VResMS;            %-Filehandle - Hyperparameter
SPM.Vhrf       = Vhrf;              % File handel HRF 
SPM.VM         = VM;                %-Filehandle - Mask

SPM.xVi        = xVi;               % non-sphericity structure
% SPM.xVi.CY     = CY;                %-<(Y - <Y>)*(Y - <Y>)'>

SPM.xX         = Xconv;                %-design structure

SPM.xM         = xM;                %-mask structure

SPM.xCon       = struct([]);        %-contrast structure

SPM.SPMid      = SPMid;
SPM.swd        = pwd;


%-Save analysis parameters in SPM.mat file
%--------------------------------------------------------------------------
if spm_matlab_version_chk('7') >=0
    save('SPM','SPM','-V6');
else
    save('SPM','SPM');
end

%==========================================================================
%- E N D: Cleanup GUI
%==========================================================================
fprintf('%s%30s\n',repmat(sprintf('\b'),1,30),'...done')                %-#
spm('FigName','Stats: done',Finter); spm('Pointer','Arrow')
fprintf('%-40s: %30s\n','Completed',spm('time'))                        %-#
fprintf('...use the results section for assessment\n\n')                %-#
