********************************************************************************
* DESCRIPTION: Analyzes time series of the firm-year pay distribution.
*
* CREDIT:      Please cite: Engbom, Niklas & Christian Moser & Jan Sauermann,
*              2021. "Firm Pay Dynamics," Working Paper.
*
* OUTPUT:      - Dataset with mean and percentiles of three firm pay measures (6 versions: unbalanced vs. semi-balanced vs. balanced x weighted vs. unweighted).
*              - Graph showing evolution of means of three firm pay measures (6 versions: unbalanced vs. semi-balanced vs. balanced x weighted vs. unweighted).
*              - Graph showing evolution of normalized means of three firm pay measures (6 versions: unbalanced vs. semi-balanced vs. balanced x weighted vs. unweighted).
*              - Graph showing evolution of percentiles of mean firm pay (6 versions: unbalanced vs. semi-balanced vs. balanced x weighted vs. unweighted).
*              - Graph showing evolution of percentiles of firm fixed effects (6 versions: unbalanced vs. semi-balanced vs. balanced x weighted vs. unweighted).
*              - Graph showing evolution of percentiles of firm-year fixed effects (6 versions: unbalanced vs. semi-balanced vs. balanced x weighted vs. unweighted).
*              - Graph showing evolution of normalized percentiles of mean firm pay (6 versions: unbalanced vs. semi-balanced vs. balanced x weighted vs. unweighted).
*              - Graph showing evolution of normalized percentiles of firm fixed effects (6 versions: unbalanced vs. semi-balanced vs. balanced x weighted vs. unweighted).
*              - Graph showing evolution of normalized percentiles of firm-year fixed effects (6 versions: unbalanced vs. semi-balanced vs. balanced x weighted vs. unweighted).
*
* TIME STAMP:  October 23, 2021.
********************************************************************************


*** loop through different individual income concepts and minimum firm size thresholds
foreach inc_concept of global inc_concept_list {
	foreach thresh of global AKM_threshold_list {
	
		disp _newline(5)
		disp "--> inc_concept = `inc_concept', thresh = `thresh'"

		
		*** sanity checks
		* load combined data
		use ${year} ${fy_fe} N using "${DIR_TEMP}/akm_estimates_combined_firm_`inc_concept'_`thresh'_${year_min}_${year_max}.dta", clear
		
		* set first year for normalization
		sum ${year}, meanonly
		local year_norm = r(min)
		
		* check that the variance of the "year" effects is small
		disp _newline(5)
		disp `"--> Check that weighted variance of "year" effects is small:"'
		if "${gtools}" == "" {
			gen float ${fy_fe}_weighted = ${fy_fe}*N
			bys ${year}: egen long N_sum = total(N)
			bys ${year}: egen float ${fy_fe}_weighted_sum = total(${fy_fe}_weighted)
			gen float ${fy_fe}_mean = ${fy_fe}_weighted_sum/N_sum
			drop ${fy_fe}_weighted ${fy_fe}_weighted_sum N_sum
		}
		else gegen float ${fy_fe}_mean = mean(${fy_fe}) [aw=N], by(${year})
		sum ${fy_fe}_mean, d
		drop ${fy_fe}_mean
		
		
		*** plot evolution of variance of wage components in firm-year FEs model
		* load data
		use `inc_concept' `inc_concept'${demeaned} ${pe} ${fy_fe} ${year} using "${DIR_TEMP}/akm_estimates_fy_`inc_concept'_`thresh'_${year_min}_${year_max}.dta", clear
		
		* compute covariance between person FE and firm-year FE each year
		gen float cov_pe_fe = .
		forvalues yy = $year_min/$year_max {
			corr ${pe} ${fy_fe} if ${year} == `yy', cov
			replace cov_pe_fe = r(cov_12) if ${year} == `yy'
		}

		* collapse to yearly level
		local collapse_str_sd = ""
		foreach var of varlist `inc_concept'* $pe $fy_fe {
			local collapse_str_sd = "`collapse_str_sd' (sd) `var'_var=`var'"
		}
		local collapse_str_p = ""
		foreach var of varlist `inc_concept'* {
			foreach p in 5 10 25 50 75 90 95 {
				local collapse_str_p = "`collapse_str_p' (p`p') `var'_p`p'=`var'"
			}
		}
		
		${gtools}collapse ///
			(firstnm) cov_pe_fe ///
			`collapse_str_sd' ///
			`collapse_str_p' ///
			, by(${year}) fast
		foreach var of varlist *_var {
			replace `var' = `var'^2
		}
		replace cov_pe_fe = 2*cov_pe_fe
		
		* normalize percentiles to 0 in first year
		foreach var of varlist `inc_concept'* {
			sum `var' if ${year} == `year_norm', meanonly
			gen float `var'_norm = `var' - r(mean)
		}
		
		* export data
		export delim using "${DIR_EXPORT}/csv/8_time_series_var_perc.csv", replace
		
		* plot evolution of variance of log earnings
		tw ///
			(connected `inc_concept'_var ${year}, sort lcolor(blue) lpattern(l _ -) mcolor(blue) msymbol(Oh)) ///
			, graphregion(color(white)) plotregion(lcolor(black) margin(l=0 r=0 b=0 t=0)) ///
			xlabel(, gmin gmax grid gstyle(dot)) ylabel(0(.1).5, format(%9.1f) gmin gmax grid gstyle(dot)) ///
			yscale(range(0 .5)) ///		
			xtitle("") ytitle("Variance of log earnings") ///
			legend(order(1 "Earnings") cols(1) region(lcolor(white))) ///
			name(ts_`inc_concept'_var, replace)
		graph_export_eps_pdf "8_time_series_`inc_concept'_var"
		
		* plot evolution of variance of (demeaned) log earnings
		tw ///
			(connected `inc_concept'${demeaned}_var ${year}, sort lcolor(blue) lpattern(l _ -) mcolor(blue) msymbol(Oh)) ///
			, graphregion(color(white)) plotregion(lcolor(black) margin(l=0 r=0 b=0 t=0)) ///
			xlabel(, gmin gmax grid gstyle(dot)) ylabel(, gmin gmax grid gstyle(dot)) ///
			xtitle("") ytitle("Variance of demeaned log earnings") ///
			legend(order(1 "Earnings") cols(1) region(lcolor(white))) ///
			name(ts_`inc_concept'${demeaned}_var, replace)
		graph_export_eps_pdf "8_time_series_`inc_concept'${demeaned}_var"
	
		
		* plot evolution of real earnings percentiles
		tw ///
			(connected `inc_concept'_p5 `inc_concept'_p10 `inc_concept'_p25 `inc_concept'_p50 `inc_concept'_p75 `inc_concept'_p90 `inc_concept'_p95 ${year}, sort lcolor(blue red green orange cyan magenta lime gold purple) lpattern(solid longdash dash shortdash dot longdash_dot dash_dot shortdash_dot -.) mcolor(blue red green orange cyan magenta lime gold purple) msymbol(O D T S Oh Dh Th Sh X)) ///
			, graphregion(color(white)) plotregion(lcolor(black) margin(l=0 r=0 b=0 t=0)) ///
			xlabel(, gmin gmax grid gstyle(dot)) ylabel(, gmin gmax grid gstyle(dot)) ///
			xtitle("") ytitle("Normal. percentiles of mean earnings (`year_norm' = 0.0)") ///
			legend(order(1 "P5" 2 "P10" 3 "P25" 4 "P50" 5 "P75" 6 "P90" 7 "P95") cols(4) region(lcolor(white)))  ///
			name(ts_perc_`inc_concept', replace)
		graph_export_eps_pdf "8_time_series_perc_`inc_concept'"
		
		* plot evolution of normalized real earnings percentiles
		tw ///
			(connected `inc_concept'_p5_norm `inc_concept'_p10_norm `inc_concept'_p25_norm `inc_concept'_p50_norm `inc_concept'_p75_norm `inc_concept'_p90_norm `inc_concept'_p95_norm ${year}, sort lcolor(blue red green orange cyan magenta lime gold purple) lpattern(solid longdash dash shortdash dot longdash_dot dash_dot shortdash_dot -.) mcolor(blue red green orange cyan magenta lime gold purple) msymbol(O D T S Oh Dh Th Sh X)) ///
			, graphregion(color(white)) plotregion(lcolor(black) margin(l=0 r=0 b=0 t=0)) ///
			xlabel(, gmin gmax grid gstyle(dot)) ylabel(, gmin gmax grid gstyle(dot)) ///
			xtitle("") ytitle("Normal. percentiles of mean earnings (`year_norm' = 0.0)") ///
			legend(order(1 "P5" 2 "P10" 3 "P25" 4 "P50" 5 "P75" 6 "P90" 7 "P95") cols(4) region(lcolor(white)))  ///
			name(ts_perc_norm_`inc_concept', replace)
		graph_export_eps_pdf "8_time_series_perc_norm_`inc_concept'"
		
		* plot graph of AKM decomposition
		tw ///
			(connected `inc_concept'${demeaned}_var ${pe}_var ${fy_fe}_var cov_pe_fe ${year}, sort lcolor(blue red green orange) lpattern(solid longdash dash shortdash) mcolor(blue red green orange) msymbol(O D T S)) ///
			, graphregion(color(white)) plotregion(lcolor(black) margin(l=0 r=0 b=0 t=0)) ///
			xlabel(, gmin gmax grid gstyle(dot)) ylabel(, gmin gmax grid gstyle(dot)) ///
			xtitle("") ytitle("Variance of pay") ///
			legend(order(1 "Total" 2 "Person" 3 "Firm-year" 4 "2*Covariance") cols(4) region(lcolor(white))) ///
			name(akm, replace)
		graph_export_eps_pdf "8_time_series_AKM"
		

		*** plot evolution of variance of AKM components in the sub-period specification
		* JS foreach y1 in 1985 1993 2000 2008 {
		foreach y1 in 1985 1993 2001 2009 {
			local y2 = `y1'+7 
			
			* load data
			use "${DIR_TEMP}/akm_estimates_fy_`inc_concept'_`thresh'_`y1'_`y2'.dta", clear
			
			* compute covariance between worker and firm FEs each year
			gen float cov_pe_fe = .
			forvalues yy = `y1'/`y2' {
				corr ${pe} ${fy_fe} if ${year} == `yy', cov
				replace cov_pe_fe = r(cov_12) if ${year} == `yy'
			}

			* collapse to yearly time series
			${gtools}collapse ///
				(firstnm) cov_pe_fe ///
				(sd) `inc_concept'${demeaned}_var=`inc_concept'${demeaned} ${fy_fe}_var=${fy_fe} ${pe}_var=${pe} ///
				, by(${year}) fast
			foreach var of varlist *_var {
				replace `var' = `var'^2
			}
			replace cov_pe_fe = 2*cov_pe_fe
			compress
			save "${DIR_TEMP}/akm_temp`y1'", replace
		}
		
		* append data
		* JS foreach y1 in 1985 1993 2000 2008 {
		foreach y1 in 1985 1993 2001 2009 { 
			if `y1' == 1985 use "${DIR_TEMP}/akm_temp`y1'", clear
			else append using "${DIR_TEMP}/akm_temp`y1'"
			rm "${DIR_TEMP}/akm_temp`y1'.dta"
		}
		
		* plot graph
		tw ///
			(connected `inc_concept'${demeaned}_var ${pe}_var ${fy_fe}_var cov_pe_fe ${year}, sort lcolor(blue red green orange) lpattern(solid longdash dash shortdash) mcolor(blue red green orange) msymbol(O D T S)) ///
			, graphregion(color(white)) plotregion(lcolor(black) margin(l=0 r=0 b=0 t=0)) ///
			xlabel(1980(10)2020, gmin gmax grid gstyle(dot)) ylabel(, gmin gmax grid gstyle(dot)) ///
			xtitle("") ytitle("Variance of pay") ///
			legend(order(1 "Total" 2 "Person" 3 "Firm-year" 4 "2*Covariance") cols(4) region(lcolor(white))) ///
			name(akm_sub, replace)
		graph_export_eps_pdf "8_time_series_AKM_subperiod"
		
		
	
		*** plot evolution of permanent vs. transitory firm pay inequality at the worker level
		* load data
		use "${DIR_TEMP}/akm_estimates_fy_`inc_concept'_`thresh'_${year_min}_${year_max}.dta", clear
		
		* compute first autocovariance as measure of permanent inequality in firm-year pay
		xtset ${id_pers} ${year}
		foreach var in ${fy_fe} {
			gen float `var'_d1 = `var'-L.`var'
			gen float `var'_d5 = `var'-L5.`var'
			gen float L_`var' = L.`var'			
			gen float `var'_autocov_perm = .
			local year_min_plus_1 = ${year_min} + 1
			forvalues yy = `year_min_plus_1'/${year_max} {
				corr `var' L_`var' if ${year} == `yy', cov
				replace `var'_autocov_perm = r(cov_12) if ${year} == `yy'
			}
			drop L_`var'			
		}
		
		* collapse to yearly time series
		${gtools}collapse ///
			(firstnm) ${fy_fe}_autocov_perm ///
			(sd) ${fy_fe}_var=${fy_fe} ${fy_fe}_d1_var=${fy_fe}_d1 ${fy_fe}_d5_var=${fy_fe}_d5 ///
			, by(${year}) fast
		foreach var of varlist *_var {
			replace `var' = `var'^2
		}
		
		* compute temporary component
		gen float ${fy_fe}_autocov_temp = ${fy_fe}_var - ${fy_fe}_autocov_perm
				
		* plot graph
		tw ///
			(connected ${fy_fe}_var ${fy_fe}_autocov_perm ${fy_fe}_autocov_temp ${year}, sort lcolor(blue red green orange) lpattern(solid longdash dash shortdash) mcolor(blue red green orange) msymbol(O D T S)) ///
			, graphregion(color(white)) plotregion(lcolor(black) margin(l=0 r=0 b=0 t=0)) ///
			xlabel(, gmin gmax grid gstyle(dot)) ylabel(, gmin gmax grid gstyle(dot)) ///
			xtitle("") ytitle("Variance of firm pay") ///
			legend(order(1 "Total" 2 "Permanent" 3 "Transitory") cols(4) region(lcolor(white))) ///
			name(permanent_transitory, replace)
		graph_export_eps_pdf "8_time_series_workers_permanent_transitory"
		
		tw ///
			(connected ${fy_fe}_d1_var ${fy_fe}_d5_var ${year}, sort lcolor(blue red green orange) lpattern(solid longdash dash shortdash) mcolor(blue red green orange) msymbol(O D T S)) ///
			, graphregion(color(white)) plotregion(lcolor(black) margin(l=0 r=0 b=0 t=0)) ///
			xlabel(, gmin gmax grid gstyle(dot)) ylabel(, gmin gmax grid gstyle(dot)) ///
			xtitle("") ytitle("Variance of changes in firm pay") ///
			legend(order(1 "1-year" 2 "5-year") cols(2) region(lcolor(white))) ///
			name(permanent_transitory, replace)
		graph_export_eps_pdf "8_time_series_workers_5_vs_1_year"		
		
		
		*** plot evolution of permanent vs. transitory firm pay inequality at the firm level
		* load data
		use ${id_firm} ${year} ${fy_fe} N using "${DIR_TEMP}/akm_estimates_combined_firm_`inc_concept'_`thresh'_${year_min}_${year_max}.dta", clear
		
		* compute first autocovariance as measure of permanent inequality in firm-year pay: Let fe = fep+fet, then  Covar(fe1, fe2) = Covar(fep1+fet1, fep2+fet2) = Covar(fep1 ,fep2) + covar(fet1, fet2) = Covar(fep1 ,fep2).
		xtset ${id_firm} ${year}
		foreach var in ${fy_fe} {
			gen float `var'_d1 = `var'-L.`var'
			gen float `var'_d5 = `var'-L5.`var'		
			gen float L_`var' = L.`var'
			gen float `var'_autocov_perm = .
			gen float `var'_autocorr = .
			local year_min_plus_1 = ${year_min} + 1
			forvalues yy = `year_min_plus_1'/${year_max} {
				corr `var' L_`var' if ${year} == `yy' [aw=N], cov
				replace `var'_autocov_perm = r(cov_12) if ${year} == `yy'
				corr `var' L_`var' if ${year} == `yy' [aw=N]
				replace `var'_autocorr = r(rho) if ${year} == `yy'
			}
			drop L_`var'			
		}
		
		* collapse to yearly time series
		${gtools}collapse ///
			(firstnm) ${fy_fe}_autocov_perm ${fy_fe}_autocorr ///
			(sd) ${fy_fe}_var=${fy_fe} ${fy_fe}_d1_var=${fy_fe}_d1 ${fy_fe}_d5_var=${fy_fe}_d5 ///
			[aw=N] ///
			, by(${year}) fast
		foreach var of varlist *_var {
			replace `var' = `var'^2
		}
		
		* compute temporary component
		gen float ${fy_fe}_autocov_temp = ${fy_fe}_var - ${fy_fe}_autocov_perm
		
		* plot graph
		tw ///
			(connected ${fy_fe}_var ${fy_fe}_autocov_perm ${fy_fe}_autocov_temp ${year}, sort lcolor(blue red green orange) lpattern(solid longdash dash shortdash) mcolor(blue red green orange) msymbol(O D T S)) ///
			, graphregion(color(white)) plotregion(lcolor(black) margin(l=0 r=0 b=0 t=0)) ///
			xlabel(, gmin gmax grid gstyle(dot)) ylabel(, gmin gmax grid gstyle(dot)) ///
			xtitle("") ytitle("Variance of firm pay") ///
			legend(order(1 "Total" 2 "Permanent" 3 "Transitory") cols(4) region(lcolor(white))) ///
			name(permanent_transitory, replace)
		graph_export_eps_pdf "8_time_series_firms_permanent_transitory"
		
		tw ///
			(connected ${fy_fe}_d1_var ${fy_fe}_d5_var ${year}, sort lcolor(blue red green orange) lpattern(solid longdash dash shortdash) mcolor(blue red green orange) msymbol(O D T S)) ///
			, graphregion(color(white)) plotregion(lcolor(black) margin(l=0 r=0 b=0 t=0)) ///
			xlabel(, gmin gmax grid gstyle(dot)) ylabel(, gmin gmax grid gstyle(dot)) ///
			xtitle("") ytitle("Variance of changes in firm pay") ///
			legend(order(1 "1-year" 2 "5-year") cols(2) region(lcolor(white))) ///
			name(permanent_transitory, replace)
		graph_export_eps_pdf "8_time_series_firms_5_vs_1_year"			
		
		tw ///
			(connected ${fy_fe}_autocorr ${year}, sort lcolor(blue red green orange) lpattern(solid longdash dash shortdash) mcolor(blue red green orange) msymbol(O D T S)) ///
			, graphregion(color(white)) plotregion(lcolor(black) margin(l=0 r=0 b=0 t=0)) ///
			xlabel(, gmin gmax grid gstyle(dot)) ylabel(, gmin gmax grid gstyle(dot)) ///
			xtitle("") ytitle("Autocorrelation of firm pay") ///
			legend(order(1 "Autocorrelation") cols(1) region(lcolor(white)) on) ///
			name(permanent_transitory, replace)
		graph_export_eps_pdf "8_time_series_firms_autocorrelation"
		
		
		*** plot time series of mean, variance, and percentiles
		* load data
		use ${id_firm} ${year} `inc_concept' ${f_fe} ${fy_fe} N using "${DIR_TEMP}/akm_estimates_combined_firm_`inc_concept'_`thresh'_${year_min}_${year_max}.dta", clear
		
		* demean firm pay variables each year
		foreach var in `inc_concept' ${f_fe} ${fy_fe} {
			if "${gtools}" == "" {
				gen float `var'_weighted = `var'*N
				bys ${year}: egen long N_sum = total(N)
				bys ${year}: egen long `var'_weighted_sum = total(`var'_weighted)
				gen float `var'_mean = `var'_weighted_sum/N_sum
				drop `var'_weighted `var'_weighted_sum N_sum
			}
			else gegen float `var'_mean = mean(`var') [aw=N], by(${year})
			replace `var' = `var' - `var'_mean
			drop `var'_mean
		}
		
		* save temporary data
		prog_comp_desc_sum_save "${DIR_TEMP}/temp_time_series.dta"

		* loop through different selection criteria and weights
		foreach balance in "unbalanced" "semi_balanced" "balanced" {
			foreach weight in "unweighted" "weighted" {
			
				* load temporary data
				local vars_list = "${id_firm} ${year} ${f_fe} ${fy_fe} `inc_concept'"
				if "`weight'" == "unweighted" use `vars_list' using "${DIR_TEMP}/temp_time_series.dta", clear
				else if "`weight'" == "weighted" use `vars_list' N using "${DIR_TEMP}/temp_time_series.dta", clear
				
				* use weights
				if "`weight'" == "unweighted" local aw = ""
				else if "`weight'" == "weighted" local aw = "[aw=N]"
				
				* impose balancedness of panel
				if "`balance'" == "semi_balanced" {
					${gtools}levelsof ${year}, local(year_list)
					local N_years: word count `year_list'
					bys ${id_firm}: keep if _N >= `N_years'/2
				}
				else if "`balance'" == "balanced" {
					${gtools}levelsof ${year}, local(year_list)
					local N_years: word count `year_list'
					bys ${id_firm}: keep if _N == `N_years'
				}
				
				* proceed only if there are observations left in the data (i.e., if some firms exist for all years in the data)
				if _N > 0 {
					
					* collapse to quantiles of variable on x-axis
					${gtools}collapse ///
						/* JS (mean) `inc_concept'_mean ${f_fe}_mean ${fy_fe}_mean */ ///
						(mean) `inc_concept'_mean=`inc_concept' ${f_fe}_mean=${f_fe} ${fy_fe}_mean=${fy_fe} ///
						(sd) `inc_concept'_var=`inc_concept' ${f_fe}_var=${f_fe} ${fy_fe}_var=${fy_fe} ///
						(p5) `inc_concept'_p5=`inc_concept' ${f_fe}_p5=${f_fe} ${fy_fe}_p5=${fy_fe} ///
						(p10) `inc_concept'_p10=`inc_concept' ${f_fe}_p10=${f_fe} ${fy_fe}_p10=${fy_fe} ///
						(p25) `inc_concept'_p25=`inc_concept' ${f_fe}_p25=${f_fe} ${fy_fe}_p25=${fy_fe} ///
						(p50) `inc_concept'_p50=`inc_concept' ${f_fe}_p50=${f_fe} ${fy_fe}_p50=${fy_fe} ///
						(p75) `inc_concept'_p75=`inc_concept' ${f_fe}_p75=${f_fe} ${fy_fe}_p75=${fy_fe} ///
						(p90) `inc_concept'_p90=`inc_concept' ${f_fe}_p90=${f_fe} ${fy_fe}_p90=${fy_fe} ///
						(p95) `inc_concept'_p95=`inc_concept' ${f_fe}_p95=${f_fe} ${fy_fe}_p95=${fy_fe} ///
						(count) N=`inc_concept' ///
						`aw', by(${year}) fast
					
					* construct variances out of standard deviations
					foreach var of varlist *_var {
						replace `var' = `var'^2
					}
					
					* normalize income measures to = 0 in first year
					foreach var of varlist `inc_concept'* ${f_fe}* ${fy_fe}* {
						sum `var' if ${year} == `year_norm', meanonly
						gen float `var'_norm = `var' - r(mean)
					}
					
					* export collapsed data
					export delim using "${DIR_EXPORT}/csv/8_time_series_`balance'_`weight'.csv", replace
					
					* plot evolution of means
					tw ///
						(connected `inc_concept'_mean ${f_fe}_mean ${fy_fe}_mean ${year}, sort lcolor(blue red green) lpattern(l _ -) mcolor(blue red green) msymbol(O D T)) ///
						, graphregion(color(white)) plotregion(lcolor(black) margin(l=0 r=0 b=0 t=0)) ///
						xlabel(, gmin gmax grid gstyle(dot)) ylabel(, gmin gmax grid gstyle(dot)) ///
						xtitle("") ytitle("Mean income measure") ///
						legend(order(1 "Earnings" 2 "Firm FEs" 3 "Firm-year FEs") cols(3) region(lcolor(white))) ///
						name(ts_mean, replace)
					graph_export_eps_pdf "8_time_series_mean_`balance'_`weight'"
					
					* plot normalized evolution of means
					tw ///
						(connected `inc_concept'_mean_norm ${f_fe}_mean_norm ${fy_fe}_mean_norm ${year}, sort lcolor(blue red green) lpattern(l _ -) mcolor(blue red green) msymbol(Oh Dh Th)) ///
						, graphregion(color(white)) plotregion(lcolor(black) margin(l=0 r=0 b=0 t=0)) ///
						xlabel(, gmin gmax grid gstyle(dot)) ylabel(, gmin gmax grid gstyle(dot)) ///
						xtitle("") ytitle("Normalized mean income measure (`year_norm' = 0.0)") ///
						legend(order(1 "Earnings" 2 "Firm FEs" 3 "Firm-year FEs") cols(3) region(lcolor(white))) ///
						name(ts_mean_norm, replace)
					graph_export_eps_pdf "8_time_series_mean_norm_`balance'_`weight'"
					
					* plot evolution of firm pay variances
					tw ///
						(connected `inc_concept'_var ${f_fe}_var ${fy_fe}_var ${year}, sort lcolor(blue red green) lpattern(l _ -) mcolor(blue red green) msymbol(O D T)) ///
						, graphregion(color(white)) plotregion(lcolor(black) margin(l=0 r=0 b=0 t=0)) ///
						xlabel(, gmin gmax grid gstyle(dot)) ylabel(, gmin gmax grid gstyle(dot)) ///
						xtitle("") ytitle("Variance of income measure") ///
						legend(order(1 "Earnings" 2 "Firm FEs" 3 "Firm-year FEs") cols(3) region(lcolor(white))) ///
						name(ts_var_firm, replace)
					graph_export_eps_pdf "8_time_series_var_firm_`balance'_`weight'"
					
					* plot normalized evolution of normalized firm pay variances
					tw ///
						(connected `inc_concept'_var_norm ${f_fe}_var_norm ${fy_fe}_var_norm ${year}, sort lcolor(blue red green) lpattern(l _ -) mcolor(blue red green) msymbol(Oh Dh Th)) ///
						, graphregion(color(white)) plotregion(lcolor(black) margin(l=0 r=0 b=0 t=0)) ///
						xlabel(, gmin gmax grid gstyle(dot)) ylabel(, gmin gmax grid gstyle(dot)) ///
						xtitle("") ytitle("Normalized variance of income measure (`year_norm' = 0.0)") ///
						legend(order(1 "Earnings" 2 "Firm FEs" 3 "Firm-year FEs") cols(3) region(lcolor(white))) ///
						name(ts_var_firm_norm, replace)
					graph_export_eps_pdf "8_time_series_var_firm_norm_`balance'_`weight'"
					
					* plot evolution of percentiles of firm-level mean earnings
					tw ///
						(connected `inc_concept'_p5 `inc_concept'_p10 `inc_concept'_p25 `inc_concept'_p50 `inc_concept'_p75 `inc_concept'_p90 `inc_concept'_p95 ${year}, sort lcolor(blue red green orange cyan magenta lime gold purple) lpattern(solid longdash dash shortdash dot longdash_dot dash_dot shortdash_dot -.) mcolor(blue red green orange cyan magenta lime gold purple) msymbol(O D T S Oh Dh Th Sh X)) ///
						, graphregion(color(white)) plotregion(lcolor(black) margin(l=0 r=0 b=0 t=0)) ///
						xlabel(, gmin gmax grid gstyle(dot)) ylabel(, gmin gmax grid gstyle(dot)) ///
						xtitle("") ytitle("Percentiles of mean earnings") ///
						legend(order(1 "P5" 2 "P10" 3 "P25" 4 "P50" 5 "P75" 6 "P90" 7 "P95") cols(4) region(lcolor(white))) ///
						name(ts_perc_firm_`inc_concept', replace)
					graph_export_eps_pdf "8_time_series_perc_firm_`inc_concept'_`balance'_`weight'"
					
					* plot evolution of percentiles of firm FEs
					tw ///
						(connected ${f_fe}_p5 ${f_fe}_p10 ${f_fe}_p25 ${f_fe}_p50 ${f_fe}_p75 ${f_fe}_p90 ${f_fe}_p95 ${year}, sort lcolor(blue red green orange cyan magenta lime gold purple) lpattern(solid longdash dash shortdash dot longdash_dot dash_dot shortdash_dot -.) mcolor(blue red green orange cyan magenta lime gold purple) msymbol(O D T S Oh Dh Th Sh X)) ///
						, graphregion(color(white)) plotregion(lcolor(black) margin(l=0 r=0 b=0 t=0)) ///
						xlabel(, gmin gmax grid gstyle(dot)) ylabel(, gmin gmax grid gstyle(dot)) ///
						xtitle("") ytitle("Percentiles of firm FEs") ///
						legend(order(1 "P5" 2 "P10" 3 "P25" 4 "P50" 5 "P75" 6 "P90" 7 "P95") cols(4) region(lcolor(white))) ///
						name(ts_perc_firm_f_fe, replace)
					graph_export_eps_pdf "8_time_series_perc_firm_${f_fe}_`balance'_`weight'"
					
					* plot evolution of percentiles of firm-year FEs
					tw ///
						(connected ${fy_fe}_p5 ${fy_fe}_p10 ${fy_fe}_p25 ${fy_fe}_p50 ${fy_fe}_p75 ${fy_fe}_p90 ${fy_fe}_p95 ${year}, sort lcolor(blue red green orange cyan magenta lime gold purple) lpattern(solid longdash dash shortdash dot longdash_dot dash_dot shortdash_dot -.) mcolor(blue red green orange cyan magenta lime gold purple) msymbol(O D T S Oh Dh Th Sh X)) ///
						, graphregion(color(white)) plotregion(lcolor(black) margin(l=0 r=0 b=0 t=0)) ///
						xlabel(, gmin gmax grid gstyle(dot)) ylabel(, gmin gmax grid gstyle(dot)) ///
						xtitle("") ytitle("Percentiles of firm-year FEs") ///
						legend(order(1 "P5" 2 "P10" 3 "P25" 4 "P50" 5 "P75" 6 "P90" 7 "P95") cols(4) region(lcolor(white))) ///
						name(ts_perc_firm_fy_fe, replace)
					graph_export_eps_pdf "8_time_series_perc_firm_${fy_fe}_`balance'_`weight'"
					
					* plot evolution of percentiles of normalized mean earnings
					tw ///
						(connected `inc_concept'_p5_norm `inc_concept'_p10_norm `inc_concept'_p25_norm `inc_concept'_p50_norm `inc_concept'_p75_norm `inc_concept'_p90_norm `inc_concept'_p95_norm ${year}, sort lcolor(blue red green orange cyan magenta lime gold purple) lpattern(solid longdash dash shortdash dot longdash_dot dash_dot shortdash_dot -.) mcolor(blue red green orange cyan magenta lime gold purple) msymbol(O D T S Oh Dh Th Sh X)) ///
						, graphregion(color(white)) plotregion(lcolor(black) margin(l=0 r=0 b=0 t=0)) ///
						xlabel(, gmin gmax grid gstyle(dot)) ylabel(, gmin gmax grid gstyle(dot)) ///
						xtitle("") ytitle("Normalized percentiles of mean earnings (`year_norm' = 0.0)") ///
						legend(order(1 "P5" 2 "P10" 3 "P25" 4 "P50" 5 "P75" 6 "P90" 7 "P95") cols(4) region(lcolor(white))) ///
						name(ts_perc_firm_`inc_concept'_norm, replace)
					graph_export_eps_pdf "8_time_series_perc_firm_`inc_concept'_norm_`balance'_`weight'"
					
					* plot evolution of percentiles of normalized firm FEs
					tw ///
						(connected ${f_fe}_p5_norm ${f_fe}_p10_norm ${f_fe}_p25_norm ${f_fe}_p50_norm ${f_fe}_p75_norm ${f_fe}_p90_norm ${f_fe}_p95_norm ${year}, sort lcolor(blue red green orange cyan magenta lime gold purple) lpattern(solid longdash dash shortdash dot longdash_dot dash_dot shortdash_dot -.) mcolor(blue red green orange cyan magenta lime gold purple) msymbol(O D T S Oh Dh Th Sh X)) ///
						, graphregion(color(white)) plotregion(lcolor(black) margin(l=0 r=0 b=0 t=0)) ///
						xlabel(, gmin gmax grid gstyle(dot)) ylabel(, gmin gmax grid gstyle(dot)) ///
						xtitle("") ytitle("Normalized percentiles of firm FEs (`year_norm' = 0.0)") ///
						legend(order(1 "P5" 2 "P10" 3 "P25" 4 "P50" 5 "P75" 6 "P90" 7 "P95") cols(4) region(lcolor(white))) ///
						name(ts_perc_firm_${f_fe}_norm, replace)
					graph_export_eps_pdf "8_time_series_perc_firm_${f_fe}_norm_`balance'_`weight'"
					
					* plot evolution of percentiles of normalized firm-year FEs
					tw ///
						(connected ${fy_fe}_p5_norm ${fy_fe}_p10_norm ${fy_fe}_p25_norm ${fy_fe}_p50_norm ${fy_fe}_p75_norm ${fy_fe}_p90_norm ${fy_fe}_p95_norm ${year}, sort lcolor(blue red green orange cyan magenta lime gold purple) lpattern(solid longdash dash shortdash dot longdash_dot dash_dot shortdash_dot -.) mcolor(blue red green orange cyan magenta lime gold purple) msymbol(O D T S Oh Dh Th Sh X)) ///
						, graphregion(color(white)) plotregion(lcolor(black) margin(l=0 r=0 b=0 t=0)) ///
						xlabel(, gmin gmax grid gstyle(dot)) ylabel(, gmin gmax grid gstyle(dot)) ///
						xtitle("") ytitle("Normalized percentiles of firm-year FEs (`year_norm' = 0.0)") ///
						legend(order(1 "P5" 2 "P10" 3 "P25" 4 "P50" 5 "P75" 6 "P90" 7 "P95") cols(4) region(lcolor(white))) ///
						name(ts_perc_firm_${fy_fe}_norm, replace)
					graph_export_eps_pdf "8_time_series_perc_firm_${fy_fe}_norm_`balance'_`weight'"
				}
			}
		}

		
		*** additional, formatted plots
		* loop through different selection criteria and weights
		foreach balance in "unbalanced" "semi_balanced" "balanced" {
			foreach weight in "unweighted" "weighted" {
				
				* read data
				import delim using "${DIR_EXPORT}/csv/8_time_series_`balance'_`weight'.csv", clear
				
				* plot evolution of percentiles of firm-year FEs
				tw ///
					(connected ${fy_fe}_p5 ${fy_fe}_p10 ${fy_fe}_p25 ${fy_fe}_p50 ${fy_fe}_p75 ${fy_fe}_p90 ${fy_fe}_p95 ${year}, sort lcolor(blue red green orange cyan magenta lime gold purple) lpattern(solid longdash dash shortdash dot longdash_dot dash_dot shortdash_dot -.) mcolor(blue red green orange cyan magenta lime gold purple) msymbol(O D T S Oh Dh Th Sh X)) ///
					, graphregion(color(white)) plotregion(lcolor(black) margin(l=0 r=0 b=0 t=0)) ///
					xlabel(1985(5)2015, format(%4.0f) gmin gmax grid gstyle(dot)) ylabel(-.5(.1).5, format(%2.1f) gmin gmax grid gstyle(dot)) ///
					xtitle("") ytitle("Percentiles of firm-year FEs") ///
					legend(order(1 "P5" 2 "P10" 3 "P25" 4 "P50" 5 "P75" 6 "P90" 7 "P95") cols(4) region(lcolor(white))) ///
					name(c_${fy_fe}_`balance'_`weight', replace)
				graph_export_eps_pdf "8_time_series_perc_${fy_fe}_`balance'_`weight'"
				
				* plot evolution of percentiles of normalized firm-year FEs
				tw ///
					(connected ${fy_fe}_p5_norm ${fy_fe}_p10_norm ${fy_fe}_p25_norm ${fy_fe}_p50_norm ${fy_fe}_p75_norm ${fy_fe}_p90_norm ${fy_fe}_p95_norm ${year}, sort lcolor(blue red green orange cyan magenta lime gold purple) lpattern(solid longdash dash shortdash dot longdash_dot dash_dot shortdash_dot -.) mcolor(blue red green orange cyan magenta lime gold purple) msymbol(O D T S Oh Dh Th Sh X)) ///
					, graphregion(color(white)) plotregion(lcolor(black) margin(l=0 r=0 b=0 t=0)) ///
					xlabel(1985(5)2015, format(%4.0f) gmin gmax grid gstyle(dot)) ylabel(-.2(.05).2, format(%3.2f) gmin gmax grid gstyle(dot)) ///
					xtitle("") ytitle("Normalized percentiles of firm-year FEs (`year_norm' = 0.0)") ///
					legend(order(1 "P5" 2 "P10" 3 "P25" 4 "P50" 5 "P75" 6 "P90" 7 "P95") cols(4) region(lcolor(white))) ///
					name(c_${fy_fe}_n_`balance'_`weight', replace)
				graph_export_eps_pdf "8_time_series_perc_${fy_fe}_norm_`balance'_`weight'"
				
				if "`balance'" == "unbalanced" & "`weight'" == "weighted" {
				
					* plot evolution of percentiles of mean earnings
					tw ///
						(connected `inc_concept'_p5 `inc_concept'_p10 `inc_concept'_p25 `inc_concept'_p50 `inc_concept'_p75 `inc_concept'_p90 `inc_concept'_p95 ${year}, sort lcolor(blue red green orange cyan magenta lime gold purple) lpattern(solid longdash dash shortdash dot longdash_dot dash_dot shortdash_dot -.) mcolor(blue red green orange cyan magenta lime gold purple) msymbol(O D T S Oh Dh Th Sh X)) ///
						, graphregion(color(white)) plotregion(lcolor(black) margin(l=0 r=0 b=0 t=0)) ///
						xlabel(1985(5)2015, format(%4.0f) gmin gmax grid gstyle(dot)) ylabel(-.5(.1).5, format(%2.1f) gmin gmax grid gstyle(dot)) ///
						xtitle("") ytitle("Percentiles of mean earnings") ///
						legend(order(1 "P5" 2 "P10" 3 "P25" 4 "P50" 5 "P75" 6 "P90" 7 "P95") cols(4) region(lcolor(white))) ///
						name(c_`inc_concept'_`balance'_`weight', replace)
					graph_export_eps_pdf "8_time_series_perc_`inc_concept'_`balance'_`weight'"
					
					* plot evolution of percentiles of normalized earnings
					tw ///
						(connected `inc_concept'_p5_norm `inc_concept'_p10_norm `inc_concept'_p25_norm `inc_concept'_p50_norm `inc_concept'_p75_norm `inc_concept'_p90_norm `inc_concept'_p95_norm ${year}, sort lcolor(blue red green orange cyan magenta lime gold purple) lpattern(solid longdash dash shortdash dot longdash_dot dash_dot shortdash_dot -.) mcolor(blue red green orange cyan magenta lime gold purple) msymbol(O D T S Oh Dh Th Sh X)) ///
						, graphregion(color(white)) plotregion(lcolor(black) margin(l=0 r=0 b=0 t=0)) ///
						xlabel(1985(5)2015, format(%4.0f) gmin gmax grid gstyle(dot)) ylabel(-.2(.05).2, format(%3.2f) gmin gmax grid gstyle(dot)) ///
						xtitle("") ytitle("Normalized percentiles of mean earnings (`year_norm' = 0.0)") ///
						legend(order(1 "P5" 2 "P10" 3 "P25" 4 "P50" 5 "P75" 6 "P90" 7 "P95") cols(4) region(lcolor(white))) ///
						name(c_`inc_concept'_n_`balance'_`weight', replace)
					graph_export_eps_pdf "8_time_series_perc_`inc_concept'_norm_`balance'_`weight'"
				}
			}
		}

		* remove unused files
		rm "${DIR_TEMP}/temp_time_series.dta"
	}
}
