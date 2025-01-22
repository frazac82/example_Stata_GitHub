


/* CS = cause-specific mortality

syntax,	time(varname) failure(varname) aged(varname)											///
		disease(varname) [covs(string)]															///
		dfaged(integer) dftve(integer) dfrp(integer) 											///
		agemin(real) agemax(real) steps(integer) tau(real) [intp(integer 200)] frmft(string)	///
		[fcif(string)] [cifti(real 0.1)]	

yll_cs,	time(varname) 					/// Time into study in years
		failure(varname)				/// Failure variable (i.e., cause-specific deaths as nonzeros) 
		aged(varname)					/// Age for conditional estimations in years (i.e., age at diagnosis of disease and matched non-exposed people)
		disease(varname) 				/// Disease (1 = yes; 0 = no)
		covs(string) 					/// Covariables for the stpm3 Royston-Parmar model (for standardised estimate)
		dfsaged(integer)				/// Degree of freedom for spline of age at diagnosis
		dftve(integer)					/// Degree of freedom for time-varying effect of age at diagnosis, exposure, and age*exposure interaction (i.e., time-varying, non-prop hazard)
		dfrp(integer) 					/// Degree of freedom for Royston-Parmar model
		agemin(real) 					/// Age for conditional estimate in years (i.e., predictions for value of age confounder), minimum
		agemax(real) 					/// Age for conditional estimate in years (i.e., predictions for value of age confounder), maximum
		steps(integer)					/// Steps/deltas between min and max age in years for age-at-diagnosis-specific estimates
		tau(real) 						/// Max age in years for AUC of survival and for CIF (i.e., RMST/RMFT)
		[intp(integer 200)]				/// [Optional] Number of nodes for numerical integration for AUC
		frmft(string)					/// Name of frame with results of YLL (frame is replaced)
		[fcif(string)]					/// [Optional] Name of frame with results of CIF (frame is replaced)
		[cifti(real 0.1)]				/// [Optional] Time intervals in years for CIF curves (i.e., x-axis)
*/



cls
clear all

capture program drop yll_cs
program yll_cs

	version 18.0
	
	syntax,	time(varname) failure(varname) aged(varname)												///
			disease(varname) [covs(string)]																///
			dfsaged(integer) dftve(integer) dfrp(integer)												///
			agemin(real) agemax(real) steps(integer) tau(real) [intp(integer 200)] frmft(string)		///
			[fcif(string)] [cifti(real 0.1)]										
			
			qui tempfile currentdatabase
			qui save `currentdatabase', replace
			
			di ""
			di "-------------------------------------------------------------------------------------------------------------------------------------------"
			di " Regression details"		
			di "-------------------------------------------------------------------------------------------------------------------------------------------"
			di " stpm3 @ns(`aged', df(`dfsaged'))"  		_column(30) "--->  Main effect age at diagnosis" 		
			di " i.`disease'"								_column(30) "--->  Main effect exposure"				
			di " @ns(`aged', df(`dfsaged'))#i.`disease'"	_column(30) "--->  Interaction age*exposure"			
			di " `covs'"									_column(30) "--->  Confounders (if any)"							
			di " ,"											_column(30) "--->  OPTIONS"									
			di " tvc(" 										_column(30) "--->  Time-varying effects..."			
			di " @ns(`aged', df(`dfsaged'))" 				_column(30) "--->  ...for age"							
			di " i.`disease'" 								_column(30) "--->  ...for exposure"					
			di " @ns(`aged', df(`dfsaged'))#i.`disease'"	_column(30) "--->  ...for interaction age*exposure"	
			di " )"											_column(30)	"--->  Close Time-varying effects"						
			di " dftvc(`dftve')" 							_column(30) "--->  DFs for time-varying effects of age, exposure, and age*exposure interaction"
			di " scale(lncumhazard) df(`dfrp')" 			_column(30) "--->  Details of RP model"
			di "-------------------------------------------------------------------------------------------------------------------------------------------"

			qui sum `failure'
			local nkd = r(max)
			
			preserve
			qui clear
			qui set obs 1
			forval ncs = 1/`nkd' {
				qui gen var`ncs'_1  = 1 
				qui gen var`ncs'_0  = 0 
				qui gen var`ncs'_m1 = -1 
			}

			qui egen CONCexp0 = concat(var*_1 var*_0), punct(" ")
			local CONCexp0 = CONCexp0[1]

			qui egen CONCexp1 = concat(var*_0 var*_1), punct(" ")
			local CONCexp1 = CONCexp1[1]

			qui egen CONCexpD = concat(var*_m1 var*_1), punct(" ")
			local CONCexpD = CONCexpD[1]
			restore

			qui levelsof `failure' if `failure'>0, local(cause)
			
			foreach c in `cause' {
				
				stset `time', f(`failure' == `c')
				
				stpm3 @ns(`aged', df(`dfsaged'))								/// Main effect age
						  i.`disease'											/// Main effect exposure
						  @ns(`aged', df(`dfsaged'))#i.`disease'				/// Interaction age*exposure
						  `covs'												/// Confounders
						  ,														/// -- OPTIONS --
						  tvc( 													/// Time-varying effects...
						  @ns(`aged', df(`dfsaged')) 							/// ...for age
						  i.`disease' 											/// ...for exposure
						  @ns(`aged', df(`dfsaged'))#i.`disease'				/// ...for interaction age*exposure
						  )														/// 
						  dftvc(`dftve')										/// DFs for time-varying effects of age at diagnosis, exposure, and age*exposure interaction	
						  scale(lncumhazard) df(`dfrp') 						/*  Details of RP model */
				
				qui estimate store cause`c'
				
			}
			
			di as text ""
			di "-------------------------------------------------------------------------------------------------------------------------------------------"
			di " Prediction details YLL || From diagnosis at age" 		as result " `agemin'" 										///
									   as text " to diagnosis at age" 	as result " `agemax'" 										///
									   as text " every" 	 	  		as result " `steps'" 										///
									   as text " year(s). Max age for survival prediction: " as result "`tau'" as text " years."		
			di "-------------------------------------------------------------------------------------------------------------------------------------------"

			if "`fcif'" != "" {
			di "-------------------------------------------------------------------------------------------------------------------------------------------"
			di " Prediction details CIF || From diagnosis at age" 		as result " `agemin'" 										///
									   as text " to diagnosis at age" 	as result " `agemax'" 										///
									   as text " every" 	 	  		as result " `steps'" 										///
									   as text " year(s). Time intervals for CIF prediction: " as result "`cifti'" as text " years."		
			di "-------------------------------------------------------------------------------------------------------------------------------------------"			
			}
			
			qui estimates dir
			local models `r(names)'
			
			if "`covs'" == "" {
				local tcovs = "if _n == 1"
			}
			else {
				local tcovs = ""
			}
			
			forval jmm = `agemin'(`steps')`agemax' {
				
				preserve
				local max_timest = `tau' - `jmm'
				local tpn = (`max_timest'/`cifti') + 1
				qui range timep 0 `max_timest' `tpn'
				qui gen tau_rmft = `max_timest' in 1
				
				
				*******************************************************************************************************************
				***RMFT************************************************************************************************************ 				   
				standsurv `tcovs', cif rmft crmodels(`models') at1(`aged' `jmm' exp 0) at2(`aged' `jmm' exp 1)					///
								   timevar(tau_rmft) atvars(YLL_*) contrast(difference) contrastvar(DYLL_*)       nodes(`intp') ci
				
				standsurv `tcovs', cif rmft crmodels(`models') at1(`aged' `jmm' exp 0) at2(`aged' `jmm' exp 1) 					///
								   timevar(tau_rmft) atvars(Ynotexp0_*) lincom("`CONCexp0'") lincomvar(TYLL_exp0) nodes(`intp') ci

				standsurv `tcovs', cif rmft crmodels(`models') at1(`aged' `jmm' exp 0) at2(`aged' `jmm' exp 1) 					///
								   timevar(tau_rmft) atvars(Ynotexp1_*) lincom("`CONCexp1'") lincomvar(TYLL_exp1) nodes(`intp') ci
				
				standsurv `tcovs', cif rmft crmodels(`models') at1(`aged' `jmm' exp 0) at2(`aged' `jmm' exp 1) 					///
								   timevar(tau_rmft) atvars(Ynotdif_*)  lincom("`CONCexpD'") lincomvar(DTYLL) 	  nodes(`intp') ci
				*******************************************************************************************************************
				
				
				if "`fcif'" != "" { 
				
				*******************************************************************************************************************
				***CUMULATIVE INCIDENCES******************************************************************************************* 
				standsurv `tcovs', cif crmodels(`models') at1(`aged' `jmm' exp 0) at2(`aged' `jmm' exp 1) 						///
								   timevar(timep) atvars(CIF_*) contrast(difference) contrastvar(DCIF_*)        nodes(`intp') ci
								   
				standsurv `tcovs', cif crmodels(`models') at1(`aged' `jmm' exp 0) at2(`aged' `jmm' exp 1) 						///
								   timevar(timep) atvars(Cnotexp0_*) lincom("`CONCexp0'") lincomvar(TFAIL_exp0) nodes(`intp') ci		
								   
				standsurv `tcovs', cif crmodels(`models') at1(`aged' `jmm' exp 0) at2(`aged' `jmm' exp 1) 						///
								   timevar(timep) atvars(Cnotexp1_*) lincom("`CONCexp1'") lincomvar(TFAIL_exp1) nodes(`intp') ci
								   
				standsurv `tcovs', cif crmodels(`models') at1(`aged' `jmm' exp 0) at2(`aged' `jmm' exp 1) 						///
								   timevar(timep) atvars(Cnotdif_*)  lincom("`CONCexpD'") lincomvar(DTFAIL)     nodes(`intp') ci
				*******************************************************************************************************************
				
				}
				
				if "`fcif'" != "" { 
					qui keep timep* CIF* DCIF* TFAIL* DTFAIL* YLL* DYLL* TYLL* DTYLL*
				}
				else {
					qui keep timep* YLL* DYLL* TYLL* DTYLL*
				}
				
				qui missings dropobs *, force
				gen age_pred = `jmm'
				gen tmax_pred = `tau'
				qui gen age = timep + `jmm'
				order age_pred tmax_pred age, first
				drop timep
				tempfile res_age`jmm'
				qui save `res_age`jmm'', replace
				di as text " Prediction at age at diagnosis = `jmm' years (Max age at diagnosis for prediction = `agemax' years) -- $S_TIME ($S_DATE)"
				restore
				}

				clear
				forval jmm = `agemin'(`steps')`agemax' {
				append using `res_age`jmm''
				}
		 
			 
			 **********************************************************************************************************************
			 ****Results for Excess Years of Life Lost*****************************************************************************
			 **********************************************************************************************************************
			 preserve
			 keep age_pred tmax_pred YLL* DYLL* TYLL* DTYLL*
			 qui missings dropobs YLL* DYLL* TYLL* DTYLL*, force
			 
			 foreach var of varlist YLL* DYLL* TYLL* DTYLL* {
				if strpos("`var'", "_lci") == 0 & strpos("`var'", "_uci") == 0 {
					rename `var' `var'_est
				}
			 }
			 
			 rename YLL_1_* 	YLL_exp0_*
			 rename YLL_2_* 	YLL_exp1_*
			 rename DYLL_1_* 	DYLL_*
			 
			 forval jn = 0/1 {
			 	gen ERL_exp`jn'_est = (tmax_pred - age_pred) - TYLL_exp`jn'_est
				gen ERL_exp`jn'_lci = (tmax_pred - age_pred) - TYLL_exp`jn'_uci
				gen ERL_exp`jn'_uci = (tmax_pred - age_pred) - TYLL_exp`jn'_lci
				}

			foreach var of varlist YLL* {
			 	label variable `var' "Years of Life Lost"
				}
				
			foreach var of varlist DYLL* {
			 	label variable `var' "Excess (Difference) in Years of Life Lost (Exp1 - Exp0)"
				}
	
			foreach var of varlist TYLL* {
			 	label variable `var' "Total Years of Life Lost"
				}				 
			 
			 foreach var of varlist DTYLL* {
			 	label variable `var' "Excess (Difference) in Total Years of Life Lost (Exp1 - Exp0)"
				}			 
			 
			foreach var of varlist ERL* {
			 	label variable `var' "Estimated Residual Years of Life (Life Expectancy)"
				}	
			 
			 label variable age_pred 	"Conditional age"
			 label variable tmax_pred 	"Max age for survival/failure predictions"
			 qui frame copy default "`frmft'", replace
			 restore
			 
			 
			 **********************************************************************************************************************
			 ****Results for Cumulative Incidence Function*************************************************************************
			 **********************************************************************************************************************
			 if "`fcif'" != "" { 
			 
			 preserve
			 drop YLL* DYLL* TYLL* DTYLL*
			 foreach var of varlist CIF* DCIF* TFAIL* DTFAIL* {
				if strpos("`var'", "_lci") == 0 & strpos("`var'", "_uci") == 0 {
					rename `var' `var'_est
				}
			 }
			 
			 rename CIF_1_* 	CIF_exp0_*
			 rename CIF_2_* 	CIF_exp1_*
			 rename DCIF_1_* 	DCIF_*
			 
			foreach var of varlist CIF* {
			 	label variable `var' "Cumulative Incidence"
				}
				
			foreach var of varlist DCIF* {
			 	label variable `var' "Difference Cumulative Incidence (Exp1 - Exp0)"
				}
	
			foreach var of varlist TFAIL* {
			 	label variable `var' "Total Failure"
				}				 
			 
			 foreach var of varlist DTFAIL* {
			 	label variable `var' "Difference Total Failure (Exp1 - Exp0)"
				}			 
			 
			 
			 label variable age_pred 	"Conditional age"
			 label variable tmax_pred 	"Max age for survival/failure predictions"
			 label variable age 		"Age from conditional age (i.e., age at diagnosis) until tmax"
			 qui frame copy default "`fcif'", replace
			 restore
			 
			 }

			 use `currentdatabase', clear
end
