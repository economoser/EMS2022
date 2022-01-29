********************************************************************************
* DESCRIPTION: Produces output tables to be included in the paper.
*
* CREDIT:      Please cite: Engbom, Niklas & Christian Moser & Jan Sauermann,
*              2021. "Firm Pay Dynamics," Working Paper.
*
* TIME STAMP:  January 7, 2022.
********************************************************************************

* sections to run
global table1 = 0 // TABLE 1: Summary stats
*				  // TABLE 2: AKM results -- produced in matlab
*				  // TABLE 3: AKM results by length of panel -- produced in matlab
global table4 = 0 // TABLE 4: Autocorrelation by year
global table5 = 0 // TABLE 5: Share of firms by year
global table6 = 0 // TABLE 6: Second stage
global table7 = 0 // TABLE 7: Autocorrelation by firm age
global writetables = 1 // print tables using MATLAB? (0 = no, 1 = yes)
global graphs = 0 // plot some output


* Some definitions for figures (taken from 9_lifecycle.do)
global SizeTitle = 4
global SizeLabel = 4
global LineWidth = .7
global color1 = "0 50 255"
global color2 = "230 100 200"
global color3 = "0 200 100"
global color4 = "0 0 100" //"250 150 80"
global color5 = "150 50 150" // "150 0 200"
global color6 = "0 100 50"
global LineStyle = "con"
global LinePattern = "solid longdash dash_dot dash shortdash dot"
global MarkerSymbol = "square circle diamond triangle + X"
global MarkerSize = 1.5
global LegendSize = 3


* make sure all export directories exist
foreach folder in csv eps out pdf tex dta {
	local dir_file = "${DIR_EXPORT}/`folder'"
	${cmd_confirm_dir} "`dir_file'"
	if _rc {
		!mkdir "`dir_file'"
	}
}


*** loop through different individual income concepts and minimum firm size thresholds
foreach inc_concept of global inc_concept_list {
	foreach thresh of global AKM_threshold_list {
		disp _newline(5)
		disp "--> inc_concept = `inc_concept', thresh = `thresh'"
		
		
		
		
		*** TABLE 1: SUMMARY STATS
		if $table1 {
			
			* outcomes to store
			local outcomes = "age college female"
			
			* load worker data
			use "${DIR_TEMP}/akm_estimates_fy_`inc_concept'_`thresh'_${year_min}_${year_max}.dta", clear

			gen byte college = educ >= 6 
			gen byte female = (gender==1)
			drop gender			
			foreach var in `outcomes' {
				format `var' %9.6f
			}
			
			local collapselist = ""
			foreach var in `outcomes' {
				local collapselist = "`collapselist' (mean) `var' = `var' (sd) `var'_sd = `var'"
			}
			${gtools}collapse `collapselist', fast
			
			* outsheet to matlab
			local list = ""
			foreach var in `outcomes' {
				local list = "`list' `var' `var'_sd"
			}
			compress
			outsheet `list' using "${DIR_EXPORT}/out/Table1_workers", comma replace nolabel
				 

				 
			* firm characteristics
			local outcomes = "${f_size} ${f_va}_pw ${f_assets} ${f_assets}_pw"
			use "${DIR_TEMP}/akm_estimates_combined_firm_financials_`inc_concept'_`thresh'_${year_min}_${year_max}.dta", clear
			
			foreach var in `outcomes' {
				format `var' %9.6f
			}			
			local collapselist = ""
			foreach var in `outcomes' {
				local collapselist = "`collapselist' (mean) `var' = `var' (sd) `var'_sd = `var'"
			}
			${gtools}collapse `collapselist', fast
			
			* outsheet to matlab
			local list = ""
			foreach var in `outcomes' {
				local list = "`list' `var' `var'_sd"
			}
			compress
			outsheet `list' using "${DIR_EXPORT}/out/Table1_firms", comma replace nolabel
			
			
			
			* load worker data
			local outcomes = "workers workeryears firms firmyears"
			use "${DIR_TEMP}/akm_estimates_fy_`inc_concept'_`thresh'_${year_min}_${year_max}.dta", clear

			* store output
			bys ${id_pers} ${year}: gen double workeryears = (_n == 1)
			bys ${id_pers}: gen double workers = (_n == 1)
			bys ${id_firm} ${year}: gen double firmyears = (_n == 1)
			bys ${id_firm}: gen double firms = (_n == 1)
			
			${gtools}collapse (sum) `outcomes', fast
			
			* outsheet to matlab
			local list = ""
			foreach var in `outcomes' {
				local list = "`list' `var'"
			}
			compress
			outsheet `list' using "${DIR_EXPORT}/out/Table1_observations", comma replace nolabel
			 
		}

		
		
		
		*** TABLE 4: AUTOCORRELATION BY YEAR
		if $table4 {
			
			local firmlist = "${fy_fe} `inc_concept'_demeaned"
			use `firmlist' $id_firm $year N using "${DIR_TEMP}/akm_estimates_combined_firm_`inc_concept'_`thresh'_${year_min}_${year_max}.dta", clear
						
			local year_diff = ${year_max} - ${year_min} + 1
			bys ${id_firm}: gen byte balanced = (_N == `year_diff')
			
			* average employment as weight
			if "${gtools}" == "" bys ${id_firm}: egen float weight = mean(N)
			else gegen float weight = mean(N), by(${id_firm})
			drop N
			
			xtset ${id_firm} ${year}
			local list = ""
			local list2 = ""
			local list3 = ""
			foreach var in `firmlist' {
				foreach version in "" "_w" "_b" "_b_w" {
					forvalues yy = $year_min / $year_max {
						local y = $year_max -`yy'
						forvalues f = 0/`y' {
							if "`version'" == "" local cond = ""
							else if "`version'" == "_b" local cond = " & balanced == 1"
							else if "`version'" == "_w" local cond = " [aw=weight]"
							else if "`version'" == "_b_w" local cond = " & balanced == 1 [aw=weight]"
							qui corr `var' f`f'.`var' if $year == `yy' `cond'
							* JS qui gen float `var'_corr`version'`yy'_`f' = r(rho)
							qui gen float `var'_cor`version'`yy'_`f' = r(rho)
							qui corr `var' f`f'.`var' if $year == `yy' `cond', c
							qui gen float `var'_cov`version'`yy'_`f' = r(cov_12)
							* JS local list3 = "`list3' `var'_corr`version'`yy'_`f' `var'_cov`version'`yy'_`f'"
							local list3 = "`list3' `var'_cor`version'`yy'_`f' `var'_cov`version'`yy'_`f'"
						}
						* JS local list = "`list' `var'_corr`version'`yy'_ `var'_cov`version'`yy'_"
						* JS local list2 = "`list2' `var'_corr`version'`yy' `var'_cov`version'`yy'"
						local list = "`list' `var'_cor`version'`yy'_ `var'_cov`version'`yy'_"
						local list2 = "`list2' `var'_cor`version'`yy' `var'_cov`version'`yy'"
					}
				}
			}
			
			${gtools}collapse `list3', fast
			
			gen byte id = 1		
			${gtools}reshape long `list', i(id) j(lag)
			drop id
			foreach var in `list2' {
				ren `var'_ `var'
			}
			compress

			outsheet lag `list2' ///
					 using "${DIR_EXPORT}/out/Table4", comma replace nolabel
		
		}

		
		
		
		*** TABLE 5: SHARE OF FIRMS BY YEAR
		if $table5 {
		
			* share of survivors
			use $id_firm $year N using "${DIR_TEMP}/akm_estimates_combined_firm_`inc_concept'_`thresh'_${year_min}_${year_max}.dta", clear
			
			* average employment as weight
			if "${gtools}" == "" bys ${id_firm}: egen float weight = mean(N)
			else gegen float weight = mean(N), by(${id_firm})
			drop N
			
			gen active = 1
			reshape wide active, i( ${id_firm} ) j( ${year} )
			forvalues yy = $year_min / $year_max { 
				forvalues y = $year_min / $year_max { 
					if `y' >= `yy' gen firms_`yy'_`y' = (active`y'==1) if active`yy' == 1
				}
			}
			drop active*
		
			local list = ""
			foreach version in "" "_w" {
				forvalues yy = $year_min / $year_max {
					forvalues y = $year_min / $year_max {
						if `y' >= `yy' {
							if "`version'" == "" local cond = ""
							else if "`version'" == "_w" local cond = " [aw=weight]"
							qui sum firms_`yy'_`y' `cond'
							qui gen active`version'_`yy'_`y' = r(mean)
							local list = "`list' active`version'_`yy'_`y'"
						}
					}
				}
			}
			
			${gtools}collapse (firstnm) `list', fast
			compress
			outsheet `list' ///
					 using "${DIR_EXPORT}/out/Table5", comma replace nolabel
		
					 
		}
	
	
	
	
		*** TABLE 6: SECOND STAGE
		if $table6 {
				
			local firmoutcomes = "${f_va}_pw"
						
			use "${DIR_TEMP}/akm_estimates_combined_firm_financials_`inc_concept'_`thresh'_${year_min}_${year_max}.dta" if year >= 1997, clear
	
			*** run regressions in levels
			foreach var in `firmoutcomes' {
				keep if `var' < .
			}
			
			* compute employment-weight as average firm employment
			if "${gtools}" == "" bys ${id_firm}: egen float weight = mean(N)
			else gegen float weight = mean(N), by(${id_firm})
			drop N
			
			* univariate regressions
			local firmcollapse = ""
			local reshapelist = ""
			foreach var in `firmoutcomes' {
				xtset ${id_firm} ${year}
				forvalues i=1/4 {
					gen f`var'`i' = f`i'.`var'
				}
				forvalues i=0/10 {
					gen l`var'`i' = l`i'.`var'
				}				
				foreach version in "uw" "w" {
					if "`version'" == "uw" qui reghdfe fye f`var'* l`var'*, a(${year} ${id_firm}) cluster(${id_firm} ${year})
					else if "`version'" == "w" qui reghdfe fye f`var'* l`var'* [aw=weight], a(${year} ${id_firm}) cluster(${id_firm} ${year})
					forvalues i=1/4 {
						gen float `var'_`version'_betaf`i' = _b[f`var'`i']
						gen float `var'_`version'_sef`i' = _se[f`var'`i']
						local firmcollapse = "`firmcollapse' `var'_`version'_betaf`i' `var'_`version'_sef`i'"
					}					
					forvalues i=0/10 {
						gen float `var'_`version'_betal`i' = _b[l`var'`i']
						gen float `var'_`version'_sel`i' = _se[l`var'`i']
						local firmcollapse = "`firmcollapse' `var'_`version'_betal`i' `var'_`version'_sel`i'"
					}
					gen float `var'_`version'_r2 = e(r2)
					gen float `var'_`version'_n = e(N)
					local firmcollapse = "`firmcollapse' `var'_`version'_r2 `var'_`version'_n"
					local reshapelist = "`reshapelist' `var'_`version'_beta `var'_`version'_se"
					qui gen dtemp = `var'-l.`var'
					if "`version'" == "uw" qui sum dtemp
					else if "`version'" == "w" qui sum dtemp [aw=weight]
					local std = r(sd)
					drop dtemp
					* need to also get autocorrelation of driving variable
					forvalues i=1/4 {
						if "`version'" == "uw" qui reghdfe `var' f`var'`i', a(${year}) cluster(${id_firm} ${year})
						else if "`version'" == "w" qui reghdfe `var' f`var'`i' [aw=weight], a(${year}) cluster(${id_firm} ${year})
						gen float `var'_`version'_indepf`i' = _b[f`var'`i'] * `std'
						local firmcollapse = "`firmcollapse' `var'_`version'_indepf`i'"
					}
					forvalues i=0/10 {
						if "`version'" == "uw" qui reghdfe `var' l`var'`i', a(${year}) cluster(${id_firm} ${year})
						else if "`version'" == "w" qui reghdfe `var' l`var'`i' [aw=weight], a(${year}) cluster(${id_firm} ${year})
						gen float `var'_`version'_indepl`i' = _b[l`var'`i'] * `std'
						local firmcollapse = "`firmcollapse' `var'_`version'_indepl`i'"
					}
					local reshapelist = "`reshapelist' `var'_`version'_indep"
				}
				drop f`var'* l`var'*
			}
			
			${gtools}collapse (firstnm) `firmcollapse', fast
			
			gen id = 1
			reshape long `reshapelist', i(id) j(lag) string
			drop id
			gen time = .
			forvalues lag = 1/4 {
				replace time = -`lag' if lag == "f`lag'"
			}
			forvalues lag = 0/10 {
				replace time = `lag' if lag == "l`lag'"
			}
			sort time
			drop lag
			
			compress
			save "${DIR_EXPORT}/dta/Table6dynamic", replace

			
			
			local firmoutcomes = "${f_va}_pw"
			
			* simulate the cumulative impact of a one st.d. increase in the independent 
			* variables, with an autocorrelation as in the data
			foreach var in `firmoutcomes' {
				foreach version in "uw" "w" {
					gen `var'_`version' = `var'_`version'_indep * (time >= 0) 
					* the cumulative response is the elasticity of pay in year t 
					* to productivity in year 0, but also the elasticity of pay
					* in year t to productivity in year 1, etc
					* contemporaneous impact if captured by the point estimate
					cap drop fye_`var'_`version'
					gen fye_`var'_`version' = 0
					tsset time
					forvalues i=0/10 {
						forvalues j=0/`i' {
							local k = `i'-`j'
							qui replace fye_`var'_`version' = fye_`var'_`version' + l`j'.`var'_`version'_beta * l`k'.`var'_`version' if time == `i'
						}
					}
				}
			}
			sort time
			* plot graph
			foreach var in `firmoutcomes' {
				foreach version in "uw" "w" {
					tw ///
						(connected `var'_`version' time if time >= 0, sort lcolor(black) lpattern(solid) mcolor(black) msymbol(O)) ///
						, graphregion(color(white)) plotregion(lcolor(black) margin(l=0 r=0 b=0 t=0)) ///
						xlabel(, gmin gmax grid gstyle(dot)) ylabel(, gmin gmax grid gstyle(dot)) ///
						xtitle("Years since impulse") ytitle("Productivity impulse") ///
						legend(order(1 "Productivity") cols(1) on region(lcolor(white))) ///
						name(`var'_`version', replace)
					graph_export_eps_pdf "6_impulse_`var'_`version'"
				}
			}
			* plot graph
			foreach var in `firmoutcomes' {
				foreach version in "uw" "w" {
					tw ///
						(connected `var'_`version'_beta fye_`var'_`version' time if time >= 0, sort lcolor(blue red green orange) lpattern(longdash dash) mcolor(blue red) msymbol(D T)) ///
						, graphregion(color(white)) plotregion(lcolor(black) margin(l=0 r=0 b=0 t=0)) ///
						xlabel(, gmin gmax grid gstyle(dot)) ylabel(, gmin gmax grid gstyle(dot)) ///
						xtitle("Years since impulse") ytitle("Effect on firm pay") ///
						legend(order(1 "Marginal effect" 2 "Cumulative effect") cols(2) region(lcolor(white))) ///
						name(pay_`var'_`version', replace)
					graph_export_eps_pdf "6_impulse_pay_`var'_`version'"
				}
			}
			
			
			local firmoutcomes = "${f_size} ${f_assets} ${f_va}_pw"
						
			use "${DIR_TEMP}/akm_estimates_combined_firm_financials_`inc_concept'_`thresh'_${year_min}_${year_max}.dta" if year >= 1997, clear
	
			* deal with left-censoring
// 			replace ${f_age} = min( ${f_age} , 12 )
		
			*** run regressions in levels
			foreach var in `firmoutcomes' {
				keep if `var' < .
			}
			
			* compute employment-weight as average firm employment
			if "${gtools}" == "" bys ${id_firm}: egen float weight = mean(N)
			else gegen float weight = mean(N), by(${id_firm})
			drop N
			
			* univariate regressions
			local firmcollapse = ""
			foreach var in `firmoutcomes' {
				foreach version in "uw" "w" {
					if "`version'" == "uw" qui reghdfe fye `var', a(${year}) cluster(${id_firm} ${year})
					else if "`version'" == "w" qui reghdfe fye `var' [aw=weight], a(${year}) cluster(${id_firm} ${year})
					gen float `var'_`version'_beta = _b[`var']
					gen float `var'_`version'_se = _se[`var']
					gen float `var'_`version'_r2 = e(r2)
					gen float `var'_`version'_n = e(N)
					local firmcollapse = "`firmcollapse' `var'_`version'_beta `var'_`version'_se `var'_`version'_r2 `var'_`version'_n"
				}
			}
			* joint regressions
			foreach version in "uw" "w" {
				if "`version'" == "uw" qui reghdfe fye `firmoutcomes', a(${year}) cluster(${id_firm} ${year})
				else if "`version'" == "w" qui reghdfe fye `firmoutcomes'  [aw=weight], a(${year}) cluster(${id_firm} ${year})
				foreach var in `firmoutcomes' {
					gen float joint_`var'_`version'_beta = _b[`var']
					gen float joint_`var'_`version'_se = _se[`var']
					local firmcollapse = "`firmcollapse' joint_`var'_`version'_beta joint_`var'_`version'_se"
				}
				gen float joint_`version'_r2 = e(r2)
				gen float joint_`version'_n = e(N)				
				local firmcollapse = "`firmcollapse' joint_`version'_r2 joint_`version'_n"
			}
			
			*** run regressions in changes at various horizons
			foreach lag in 1 3 5 10 {
				
				* winsorize or not
				foreach winsor in 0 1 5 10 {
				
					* create long diferences
					xtset ${id_firm} ${year}
					local laglist = ""
					foreach var in `firmoutcomes' {
						gen d`var' = `var'-l`lag'.`var'
						local laglist = "`laglist' d`var'"
					}
					gen dfye = fye-l`lag'.fye
					
					* keep only when all long differences exist
					foreach var in `firmoutcomes' {
						foreach var2 in `firmoutcomes' {
							replace d`var' = . if d`var2' == .
						}
					}
					
					* potentially winsorize independent variables
					foreach var in `firmoutcomes' {
						qui sum d`var', d
						local r1 = r(p`winsor')
						local r = 100-`winsor'
						local r2 = r(p`r')
						replace d`var' = max( min( d`var' , `r2' ) , `r1' )
					}
					
					* univariate regressions
					foreach var in `firmoutcomes' {
						foreach version in "uw" "w" {
							if "`version'" == "uw" qui reghdfe dfye d`var', a(${year}) cluster(${id_firm} ${year})
							else if "`version'" == "w" qui reghdfe dfye d`var' [aw=weight], a(${year}) cluster(${id_firm} ${year})
							gen float `var'_`version'_d`lag'_w`winsor'_beta = _b[d`var']
							gen float `var'_`version'_d`lag'_w`winsor'_se = _se[d`var']
							gen float `var'_`version'_d`lag'_w`winsor'_r2 = e(r2)
							gen float `var'_`version'_d`lag'_w`winsor'_n = e(N)
							local firmcollapse = "`firmcollapse' `var'_`version'_d`lag'_w`winsor'_beta `var'_`version'_d`lag'_w`winsor'_se `var'_`version'_d`lag'_w`winsor'_r2 `var'_`version'_d`lag'_w`winsor'_n"
						}
					}

					* joint regressions
					foreach version in "uw" "w" {
						if "`version'" == "uw" qui reghdfe dfye `laglist', a(${year}) cluster(${id_firm} ${year})
						else if "`version'" == "w" qui reghdfe dfye `laglist'  [aw=weight], a(${year}) cluster(${id_firm} ${year})
						foreach var in `firmoutcomes' {
							gen float joint_`var'_`version'_d`lag'_w`winsor'_beta = _b[d`var']
							gen float joint_`var'_`version'_d`lag'_w`winsor'_se = _se[d`var']
							local firmcollapse = "`firmcollapse' joint_`var'_`version'_d`lag'_w`winsor'_beta joint_`var'_`version'_d`lag'_w`winsor'_se"
						}
						gen float joint_`version'_d`lag'_w`winsor'_r2 = e(r2)
						gen float joint_`version'_d`lag'_w`winsor'_n = e(N)				
						local firmcollapse = "`firmcollapse' joint_`version'_d`lag'_w`winsor'_r2 joint_`version'_d`lag'_w`winsor'_n"
					}
					drop d*
				}
				
			}

			${gtools}collapse (firstnm) `firmcollapse', fast

			compress
			outsheet `firmcollapse' ///
					 using "${DIR_EXPORT}/out/Table6", comma replace nolabel
			
		}

		
		
		
		*** TABLE 7: AUTOCORRELATION BY FIRM AGE
		if $table7 {
			
			* load data
			local f_age_min = 0
			local f_age_max = 19
			use ///
				fye $id_firm $f_age N ///
				if ${f_age} >= `f_age_min' & ${f_age} <= `f_age_max' ///
				using "${DIR_TEMP}/akm_estimates_combined_firm_`inc_concept'_`thresh'_${year_min}_${year_max}.dta", clear
			
			local year_diff = `f_age_max' - `f_age_min' + 1
			bys ${id_firm}: gen byte balanced = (_N == `year_diff')
			
			* average employment as weight
			if "${gtools}" == "" bys ${id_firm}: egen float weight = mean(N)
			else gegen float weight = mean(N), by(${id_firm})
			drop N
			xtset ${id_firm} ${f_age}
			
			local list = ""
			local list2 = ""
			foreach version in "" "_w" "_balanced" "_balanced_w" {
				forvalues a = 0/19 {
					local aa = 19-`a'
					forvalues f = 0/`aa' {
						if "`version'" == "" local cond = ""
						else if "`version'" == "_balanced" local cond = " & balanced == 1"
						else if "`version'" == "_w" local cond = " [aw=weight]"
						else if "`version'" == "_balanced_w" local cond = " & balanced == 1 [aw=weight]"
						qui corr fye f`f'.fye if ${f_age} == `a' `cond'
						qui gen float autocorr`version'`a'_`f' = r(rho)
						qui corr fye f`f'.fye if ${f_age} == `a' `cond', c
						qui gen float autocov`version'`a'_`f' = r(cov_12)
					}
					local list = "`list' autocorr`version'`a'_ autocov`version'`a'_"
					local list2 = "`list2' autocorr`version'`a' autocov`version'`a'"
				}
			}
			
			${gtools}collapse autocorr* autocov*, fast
			
			gen byte id = 1		
			${gtools}reshape long `list', i(id) j(lag)
			drop id
			foreach var in `list2' {
				ren `var'_ `var'
			}
			outsheet lag `list2' ///
					 using "${DIR_EXPORT}/out/Table7", comma replace nolabel
		
		}

		
		
		
		*** WRITE TABLES IN MATLAB
		if $writetables {
		
			*** PRINT TABLES
			if inlist("`c(os)'", "MacOSX", "Unix") {
				!${FILE_MATLAB} -nosplash -noFigureWindows -nodesktop -nodisplay <"${DIR_CODE}/TablesForPaper.m" // XXX OLD SETTINGS: -nodesktop -nodisplay
			}
			else if "`c(os)'" == "Windows" {
				!matlab -nosplash -noFigureWindows -batch -wait "run '${DIR_CODE}/TablesForPaper.m'" // OLD: !matlab -nosplash -noFigureWindows -r "run '${DIR_CODE}/FUN_CONNECTED.m'; exit"
			}
			
		}
		

		
		
		if $graphs {
			
			
			insheet using "${DIR_EXPORT}/out/Table4.out", clear
				* JS foreach version in _balanced_w _balanced _w "" {
				foreach version in _b_w _b _w "" {
				foreach var in fye_cor fye_cov {
					egen t`var'`version' = rowmean(`var'`version'*)
					drop `var'`version'*
				}
			}
			* JS foreach version in _balanced_w _balanced _w "" {
			foreach version in _b_w _b _w "" {
				foreach var in fye_cor fye_cov {
					ren t`var'`version' `var'`version'
				}
			}
			
			foreach version in "" "_w" {
				foreach var in fye_cor fye_cov {
					* JS if "`var'" == "autocorr" local ytitle = "Autocorrelation"
					* JS if "`var'" == "autocov" local ytitle = "Autocovariance"
					if "`var'" == "fye_cor" local ytitle = "Autocorrelation"
					if "`var'" == "fye_cov" local ytitle = "Autocovariance"
					* JS 					tw (${LineStyle} `var'`version' `var'_balanced`version' lag, sort lcolor(blue red green orange) lpattern(solid longdash dash shortdash) mcolor(blue red green orange) msymbol(O D T S)) ///
*					tw (li `var'`version' `var'_b`version' lag, sort lcolor(blue red green orange) lpattern(solid longdash dash shortdash) mcolor(blue red green orange) msymbol(O D T S)) ///
					tw (li fye_cor fye_cor_b lag, sort lcolor(blue red green orange) lpattern(solid longdash dash shortdash) mcolor(blue red green orange) msymbol(O D T S)) ///
						, graphregion(color(white)) plotregion(lcolor(black) margin(l=0 r=0 b=0 t=0)) ///
						xlabel(, gmin gmax grid gstyle(dot)) ylabel(0(.25)1, gmin gmax grid gstyle(dot)) ///
						xtitle("Lag length (years)") ytitle("`ytitle'") ///
						legend(order(1 "Unbalanced panel" 2 "Balanced panel") cols(2) region(lcolor(none) fcolor(none)) size(${LegendSize}) margin(0 0 0 0)) ///
						name("`var'`version'", replace)
						graph export "${DIR_EXPORT}/eps/6_tables_`var'`version'.eps", replace
				}
			}
		
		}
		
	}
}
