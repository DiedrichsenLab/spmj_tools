function [U] = spm_get_ons(SPM,s)
% returns input [designed effects] structures
% FORMAT [U] = spm_get_ons(SPM,s)
%
% s  - session number (used by batch system)
%
% U     - (1 x n)   struct array of (n) trial-specific structures
%
% 	U(i).name   - cell of names for each input or cause
% 	U(i).u      - inputs or stimulus function matrix
% 	U(i).dt     - time bin (seconds)
% 	U(i).ons    - onsets    (in SPM.xBF.UNITS)
% 	U(i).dur    - durations (in SPM.xBF.UNITS)
%	U(i).P      - parameter struct.
%
% 	    U(i).P(p).name - parameter name
% 	    U(i).P(p).P    - parameter vector
% 	    U(i).P(p).h    - order of polynomial expansion
% 	    U(i).P(p).i    - sub-indices of u pertaining to P
%_______________________________________________________________________
%
%
% SLICE TIMIING
%
% With longs TRs you may want to shift the regressors so that they are
% aligned to a particular slice.  This is effected by resetting the
% values of defaults.stats.fmri.t and defaults.stats.fmri.t0 in
% spm_defaults. defaults.stats.fmri.t is the number of time-bins per
% scan used when building regressors.  Onsets are defined
% in temporal units of scans starting at 0.  defaults.stats.fmri.t0 is
% the first time-bin at which the regressors are resampled to coincide
% with data acquisition.  If defaults.stats.fmri.t0 = 1 then the
% regressors will be appropriate for the first slice.  If you want to
% temporally realign the regressors so that they match responses in the
% middle slice then make defaults.stats.fmri.t0 =
% defaults.stats.fmri.t/2 (assuming there is a negligible gap between
% volume acquisitions. Default values are defaults.stats.fmri.t = 16
% and defaults.stats.fmri.t0 = 1.
%
%
%_______________________________________________________________________
% @(#)spm_get_ons.m	2.39 Karl Friston 03/02/21

%-GUI setup
%-----------------------------------------------------------------------
spm_help('!ContextHelp',mfilename)

% time units
%-----------------------------------------------------------------------
k     = SPM.nscan(s);
T     = SPM.xBF(1).T;
dt    = SPM.xBF(1).dt;
try
	UNITS = SPM.xBF(1).UNITS;
catch
	UNITS = 'scans';
end
switch UNITS

	case 'scans'
	%----------------------------------------------------------------
	TR = T*dt;

	case 'secs'
	%----------------------------------------------------------------
	TR = 1;
end

% get inputs and names (try SPM.Sess(s).U first)
%=======================================================================
try
	U   = SPM.Sess(s).U;
	v   = length(U);
catch

	%-prompt string
	%---------------------------------------------------------------
	str = sprintf('Session %d: trial specification in %s',s,UNITS);
	spm_input(str,1,'d')

	U   = {};
	v   = spm_input('number of conditions/trials',2,'w1');
end

% get trials
%-----------------------------------------------------------------------
for i = 1:v

	% get names
	%---------------------------------------------------------------
	try
		Uname     = U(i).name(1);
	catch
		str       = sprintf('name for condition/trial %d ?',i);
		Uname     = {spm_input(str,3,'s',sprintf('trial %d',i))};
		U(i).name = Uname;
	end

	% get main [trial] effects
	%================================================================

	% onsets
	%---------------------------------------------------------------
	try
		ons = U(i).ons(:);
	catch
		ons = [];
	end
	if ~length(ons)
		str      = ['vector of onsets - ' Uname{1}];
		ons      = spm_input(str,4,'r',' ',[Inf 1]);
		U(i).ons = ons(:);

	end

	% durations
	%---------------------------------------------------------------
	try
		dur = U(i).dur(:);
	catch
		dur = [];
	end
	if ~length(dur)
		str = 'duration[s] (events = 0)';
		while 1
			dur = spm_input(str,5,'r',' ',[Inf 1]);
			if length(dur) == 1
				dur    = dur*ones(size(ons));
			end
			if length(dur) == length(ons), break, end
			str = sprintf('enter a scalar or [%d] vector',...
					length(ons));
		end
		U(i).dur = dur;
	end

	% peri-stimulus times {seconds}
	%---------------------------------------------------------------
	pst   = [1:k]*T*dt - ons(1)*TR;			
	for j = 1:length(ons)
		w      = [1:k]*T*dt - ons(j)*TR;
		v      = find(w >= -1);
		pst(v) = w(v);
	end


	% add parameters x trial interactions
	%================================================================

	% get parameter stucture xP
	%----------------------------------------------------------------
	try 
		xP          = U(i).P;
		Pname       = xP(1).name;
	catch

		Pname       = {'none','time','other'};
		Pname       = spm_input('parametric modulation',6,'b',Pname);

		switch Pname

		case 'none'
		%--------------------------------------------------------
		xP(1).name  = 'none';
		xP(1).P     = ons*TR;
		xP(1).h     = 0;

		case 'time'
		%--------------------------------------------------------
		xP(1).name  = 'time';
		xP(1).P     = ons*TR;

		case 'other'
		%--------------------------------------------------------
		str   = ['# parameters (' Uname{1} ')'];
		for q = 1:spm_input(str,7,'n1',1);

			% get names and parametric variates
			%------------------------------------------------
			str   = sprintf('parameter %d name',q);
			Pname = spm_input(str,7,'s');
			P     = spm_input(Pname,7,'r',[],[length(ons),1]);


			% sub-indices and inputs
			%------------------------------------------------
			xP(q).name  = Pname;
			xP(q).P     = P(:);

		end
		end % switch

	end % try

	% interaction with causes (u) - 1st = main effects
	%----------------------------------------------------------------
	u     = ons.^0;
    if (~strcmp('none',Pname))
	    for q = 1:length(xP)
		    xP(q).i = [1, ([1:xP(q).h] + size(u,2))];
			P   = spm_en(xP(q).P);
 		    u   = [u xP(q).P];
		    str = sprintf('%s x %s',Uname{1},xP(q).name);
		    Uname{end + 1} = str;
        end;
	end

	% and scale so sum(u*dt) = number of events, if event-related 
	%---------------------------------------------------------------
	if ~any(dur)
		u  = u/dt;
	end

	% create stimulus functions (32 bin offset)
	%===============================================================
	ton       = round(ons*TR/dt) + 32;			% onsets
	tof       = round(dur*TR/dt) + ton + 1;			% offset
	sf        = sparse((k*T + 128),size(u,2));
	for j = 1:length(ton)
		sf(ton(j),:) = sf(ton(j),:) + u(j,:);
		sf(tof(j),:) = sf(tof(j),:) - u(j,:);
	end
	sf        = cumsum(sf);					% integrate
	sf        = sf(1:(k*T + 32),:);				% stimulus

	% place in ouputs structure
	%---------------------------------------------------------------
	U(i).name = Uname;		% - input names
	U(i).dt   = dt;			% - time bin {seconds}
	U(i).u    = sf;			% - stimulus function matrix
	U(i).pst  = pst;		% - pst (seconds)
	U(i).P    = xP;			% - parameter struct

end % (v)
