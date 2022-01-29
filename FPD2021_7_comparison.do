********************************************************************************
* DESCRIPTION: Compares firm-level mean earnings, firm FEs, and firm-year FEs.
*
* CREDIT:      Please cite: Engbom, Niklas & Christian Moser & Jan Sauermann,
*              2021. "Firm Pay Dynamics," Working Paper.
*
* OUTPUT:      - 6 regression results in log file
*              - 3 graphs showing means
*              - 3 graphs showing quantiles
*              - 3 graphs showing mean quantiles
*              - 3 graphs showing quantiles of quantiles
*
* TIME STAMP:  October 23, 2021.
********************************************************************************


*** loop through different individual income concepts and minimum firm size thresholds
foreach inc_concept of global inc_concept_list {
	foreach thresh of global AKM_threshold_list {
		disp _newline(5)
		disp "--> inc_concept = `inc_concept', thresh = `thresh'"
		
		
		*** graphic comparison between firm-year FEs and firm FEs
		* load combined data
		use ${id_firm}_original ${year} `inc_concept'_demeaned ${f_fe} ${fy_fe} N using "${DIR_TEMP}/akm_estimates_combined_firm_`inc_concept'_`thresh'_${year_min}_${year_max}.dta", clear
		
		* de-mean firm pay variables each year
		foreach var in `inc_concept'_demeaned ${f_fe} ${fy_fe} {
			if "${gtools}" == "" {
				gen float `var'_weighted = `var'*N
				bys ${year}: egen long N_sum = total(N)
				bys ${year}: egen float `var'_weighted_sum = total(`var'_weighted)
				gen float `var'_mean = `var'_weighted_sum/N_sum
				drop `var'_weighted `var'_weighted_sum N_sum
			}
			else gegen float `var'_mean = mean(`var') [aw=N], by(${year})
			replace `var' = `var' - `var'_mean
			drop `var'_mean
		}

		* regression analysis
		foreach var_dep in `inc_concept'_demeaned ${f_fe} ${fy_fe} {
			foreach var_indep of varlist `inc_concept'_demeaned ${f_fe} ${fy_fe} {
				if "`var_dep'" != "`var_indep'" reghdfe `var_dep' `var_indep' [aw=N], noabsorb vce(cluster ${id_firm}_original)
			}
		}
		
		* save temporary data
		prog_comp_desc_sum_save "${DIR_TEMP}/temp_comp.dta"

		* loop through x-axis variables to create dataset and graph comparing income measures
		foreach var_x in `inc_concept'_demeaned ${f_fe} ${fy_fe} {
			
			* load temporary data
			use `inc_concept'_demeaned ${f_fe} ${fy_fe} N using "${DIR_TEMP}/temp_comp.dta", clear
			
			* generate quantiles of variable on x-axis
			foreach var of varlist `inc_concept'_demeaned ${f_fe} ${fy_fe} {
				if "${gtools}" == "" xtile `var'_q = `var' [aw=N], n(${n_bins_comparison})
				else gquantiles `var'_q = `var' [aw=N], xtile n(${n_bins_comparison})
				sum `var'_q, meanonly
				replace `var'_q = 100*(`var'_q - r(min))/(r(max) - r(min))
				local var_l: var label `var'
				label var `var'_q "Quantiles of `var_l'"
			}
			gen float q = `var_x'_q
			label var q "Quantiles of variable on x-axis"
			
			* collapse to quantiles of variable on x-axis
			${gtools}collapse ///
				(mean) `inc_concept'_demeaned_mean=`inc_concept'_demeaned `inc_concept'_demeaned_q_mean=`inc_concept'_demeaned_q ${f_fe}_mean=${f_fe} ${f_fe}_q_mean=${f_fe}_q ${fy_fe}_mean=${fy_fe} ${fy_fe}_q_mean=${fy_fe}_q ///
				(p10) `inc_concept'_demeaned_p10=`inc_concept'_demeaned `inc_concept'_demeaned_q_p10=`inc_concept'_demeaned_q ${f_fe}_p10=${f_fe} ${f_fe}_q_p10=${f_fe}_q ${fy_fe}_p10=${fy_fe} ${fy_fe}_q_p10=${fy_fe}_q ///
				(p50) `inc_concept'_demeaned_p50=`inc_concept'_demeaned `inc_concept'_demeaned_q_p50=`inc_concept'_demeaned_q ${f_fe}_p50=${f_fe} ${f_fe}_q_p50=${f_fe}_q ${fy_fe}_p50=${fy_fe} ${fy_fe}_q_p50=${fy_fe}_q ///
				(p90) `inc_concept'_demeaned_p90=`inc_concept'_demeaned `inc_concept'_demeaned_q_p90=`inc_concept'_demeaned_q ${f_fe}_p90=${f_fe} ${f_fe}_q_p90=${f_fe}_q ${fy_fe}_p90=${fy_fe} ${fy_fe}_q_p90=${fy_fe}_q ///
				(rawsum) N ///
				[aw=N] ///
				, by(q) fast
			
			* export collapsed data
			export delim using "${DIR_EXPORT}/csv/7_comp_var_x_`var_x'.csv", replace
			
			* prepare plots
			local dep_var = subinstr("`inc_concept'_demeaned ${f_fe} ${fy_fe}", "`var_x'", "", .)
			local counter = 1
			foreach var of local dep_var {
				local var_`counter' = "`var'"
				if "`var'" == "`inc_concept'_demeaned" local legend_`counter' = "Mean earnings"
				else if "`var'" == "${f_fe}" local legend_`counter' = "Firm FEs"
				else if "`var'" == "${fy_fe}" local legend_`counter' = "Firm-year FEs"
				local ++counter
			}
			if "`var_x'" == "`inc_concept'_demeaned" local indep_var_title = "Mean earnings"
			else if "`var_x'" == "${f_fe}" local indep_var_title = "Firm fixed effects"
			else if "`var_x'" == "${fy_fe}" local indep_var_title = "Firm-year fixed effects"
			
			* plot mean of variables on y-axis vs. mean of variable on x-axis
			tw ///
				(scatter `var_1'_mean `var_x'_mean, mcolor(blue) msymbol(O)) ///
				(scatter `var_2'_mean `var_x'_mean, mcolor(red) msymbol(Oh)) ///
				(lfit `var_1'_mean `var_x'_mean [aw=N], lcolor(blue)) ///
				(lfit `var_2'_mean `var_x'_mean [aw=N], lcolor(red)) ///
				(line `var_x'_mean `var_x'_mean, sort lcolor(black)) ///
				, graphregion(color(white)) plotregion(lcolor(black) margin(l=0 r=0 b=0 t=0)) ///
				xlabel(, gmin gmax grid gstyle(dot)) ylabel(, gmin gmax grid gstyle(dot)) ///
				xtitle("`indep_var_title'") ytitle("Mean firm pay measure") ///
				legend(order(1 "`legend_1'" 2 "`legend_2'" 5 "45-degree line") cols(3) region(lcolor(white))) ///
				name(comp_`var_x'_mean, replace)
			graph_export_eps_pdf "7_comp_`var_x'_mean"
			
			* plot quantiles of variables on y-axis vs. quantiles of variable on x-axis
			tw ///
				(scatter `var_1'_p10 `var_1'_p50 `var_1'_p90 `var_x'_mean, mcolor(blue midblue eltblue) msymbol(O D T)) ///
				(scatter `var_2'_p10 `var_2'_p50 `var_2'_p90 `var_x'_mean, mcolor(red orange_red orange) msymbol(Oh Dh Th)) ///
				(line `var_x'_mean `var_x'_mean, sort lcolor(black)) ///
				, graphregion(color(white)) plotregion(lcolor(black) margin(l=0 r=0 b=0 t=0)) ///
				xlabel(, gmin gmax grid gstyle(dot)) ylabel(, gmin gmax grid gstyle(dot)) ///
				xtitle("`indep_var_title'") ytitle("Quantiles of firm pay measure") ///
				legend(order(1 "P10 `legend_1'" 2 "P50 `legend_1'" 3 "P90 `legend_1'" 4 "P10 `legend_2'" 5 "P50 `legend_2'" 6 "P90 `legend_2'" 7 "45-degree line") cols(3) region(lcolor(white))) ///
				name(comp_`var_x'_perc, replace)
			graph_export_eps_pdf "7_comp_`var_x'_perc"
			
			* plot mean of quantiles of variables on y-axis vs. mean of quantiles of variable on x-axis
			tw ///
				(scatter `var_1'_q_mean `var_x'_q_mean, mcolor(blue) msymbol(O)) ///
				(scatter `var_2'_q_mean `var_x'_q_mean, mcolor(red) msymbol(Oh)) ///
				(lfit `var_1'_q_mean `var_x'_q_mean [aw=N], lcolor(blue)) ///
				(lfit `var_2'_q_mean `var_x'_q_mean [aw=N], lcolor(red)) ///
				(line `var_x'_q_mean `var_x'_q_mean, sort lcolor(black)) ///
				, graphregion(color(white)) plotregion(lcolor(black) margin(l=0 r=0 b=0 t=0)) ///
				xlabel(, gmin gmax grid gstyle(dot)) ylabel(, gmin gmax grid gstyle(dot)) ///
				xtitle("`indep_var_title' quantiles") ytitle("Mean quantiles of firm pay measure") ///
				legend(order(1 "`legend_1'" 2 "`legend_2'" 5 "45-degree line") cols(3) region(lcolor(white))) ///
				name(comp_`var_x'_mean_q, replace)
			graph_export_eps_pdf "7_comp_`var_x'_mean_q"
			
			* plot quantiles of quantiles of variables on y-axis vs. quantiles of quantiles of variable on x-axis
			tw ///
				(scatter `var_1'_q_p10 `var_1'_q_p50 `var_1'_q_p90 `var_x'_q_mean, mcolor(blue midblue eltblue) msymbol(O D T)) ///
				(scatter `var_2'_q_p10 `var_2'_q_p50 `var_2'_q_p90 `var_x'_q_mean, mcolor(red orange_red orange) msymbol(Oh Dh Th)) ///
				(line `var_x'_q_mean `var_x'_q_mean, sort lcolor(black)) ///
				, graphregion(color(white)) plotregion(lcolor(black) margin(l=0 r=0 b=0 t=0)) ///
				xlabel(, gmin gmax grid gstyle(dot)) ylabel(, gmin gmax grid gstyle(dot)) ///
				xtitle("`indep_var_title' quantiles") ytitle("Quantiles of quantiles of firm pay measure") ///
				legend(order(1 "P10 `legend_1'" 2 "P50 `legend_1'" 3 "P90 `legend_1'" 4 "P10 `legend_2'" 5 "P50 `legend_2'" 6 "P90 `legend_2'" 7 "45-degree line") cols(3) region(lcolor(white))) ///
				name(comp_`var_x'_perc_q, replace)
			graph_export_eps_pdf "7_comp_`var_x'_perc_q"
		}
		
		* remove redundant files
		rm "${DIR_TEMP}/temp_comp.dta"
		
		
// 		*** XXX OLD (February 2020):
// 		* loop through firm pay measures
// 		foreach var in earn_mean f_fe fy_fe {
			
// 			* read data
// 			use "${DIR_EXPORT}/csv/7_comp_var_x_`var_x'.csv", clear

// 			* prepare plots
// 			local dep_var = subinstr("earn_mean f_fe fy_fe", "`var'", "", .)
// 			local counter = 1
// 			foreach var2 of local dep_var {
// 				local var_`counter' = "`var2'"
// 				if "`var2'" == "earn_mean" {
// 					local legend_`counter' = "Mean earnings"
// 					local legend_`counter'_small = "mean earnings"
// 				}
// 				else if "`var2'" == "f_fe" {
// 					local legend_`counter' = "Firm FEs"
// 					local legend_`counter'_small = "firm FEs"
// 				}
// 				else if "`var2'" == "fy_fe" {
// 					local legend_`counter' = "Firm-year FEs"
// 					local legend_`counter'_small = "firm-year FEs"
// 				}
// 				local ++counter
// 			}
// 			if "`var'" == "earn_mean" local indep_var_title = "Firm-level mean earnings"
// 			else if "`var'" == "f_fe" local indep_var_title = "Firm FEs"
// 			else if "`var'" == "fy_fe" local indep_var_title = "Firm-year FEs"

// 			* plot mean of variables on y-axis vs. mean of variable on x-axis
// 			gen empty = .
// 			reg `var_1'_mean `var'_mean
// 			local b_1: di %4.3f `=_b[`var'_mean]'
// 			local se_1: di %4.3f `=_se[`var'_mean]'
// 			reg `var_2'_mean `var'_mean
// 			local b_2: di %4.3f `=_b[`var'_mean]'
// 			local se_2: di %4.3f `=_se[`var'_mean]'
// 			tw ///
// 				(scatter `var_1'_mean `var'_mean, mcolor(blue) msymbol(O)) ///
// 				(lfit `var_1'_mean `var'_mean [aw=empw_sum], lcolor(blue)) ///
// 				(scatter `var_2'_mean `var'_mean, mcolor(red) msymbol(Oh)) ///
// 				(lfit `var_2'_mean `var'_mean [aw=empw_sum], lcolor(red)) ///
// 				(function y = x, range(-.75 1) lcolor(black)) ///
// 				, text(-.3 -.13 "`legend_1' slope (s.e.) = `b_1' (`se_1')", color(blue) place(se)) ///
// 				text(-.4 -.13 "`legend_2' slope (s.e.) = `b_2' (`se_2')", color(red) place(se)) ///
// 				graphregion(color(white)) plotregion(lcolor(black) margin(l=0 r=0 b=0 t=0)) ///
// 				xlabel(-.75(.25)1, format(%3.2f) gmin gmax grid gstyle(dot)) ylabel(-.75(.25)1, format(%3.2f) gmin gmax grid gstyle(dot)) ///
// 				xtitle("`indep_var_title'") ytitle("Mean firm pay measure") ///
// 				legend(order(1 "`legend_1'" 2 "Best fit for `legend_1_small'" 3 "`legend_2'" 4 "Best fit for `legend_2_small'") cols(2) region(lcolor(white))) ///
// 				name(comp_`var'_mean, replace)
// 			graph_export_eps_pdf "7_comp_`var'_mean"

// 			* plot percentiles of variables on y-axis vs. quantiles of variable on x-axis
// 			tw ///
// 				(scatter `var_1'_p10 `var_1'_p50 `var_1'_p90 `var'_mean, mcolor(blue midblue eltblue) msymbol(O D T)) ///
// 				(scatter `var_2'_p10 `var_2'_p50 `var_2'_p90 `var'_mean, mcolor(red orange_red orange) msymbol(Oh Dh Th)) ///
// 				(function y = x, range(-.75 1) lcolor(black)) ///
// 				, graphregion(color(white)) plotregion(lcolor(black) margin(l=0 r=0 b=0 t=0)) ///
// 				xlabel(-.75(.25)1, format(%3.2f) gmin gmax grid gstyle(dot)) ylabel(-.75(.25)1, format(%3.2f) gmin gmax grid gstyle(dot)) ///
// 				xtitle("`indep_var_title'") ytitle("Quantiles of firm pay measure") ///
// 				legend(order(1 "P10 `legend_1_small'" 2 "P50 `legend_1_small'" 3 "P90 `legend_1_small'" 4 "P10 `legend_2_small'" 5 "P50 `legend_2_small'" 6 "P90 `legend_2_small'") cols(3) region(lcolor(white))) ///
// 				name(comp_`var'_perc, replace)
// 			graph_export_eps_pdf "7_comp_`var'_perc"
// 		}
		
		
		*** compare rolling time windows with firm FEs vs. firm-year FE model
		cap confirm file "${DIR_TEMP}/akm_estimates_combined_firm_`inc_concept'_`thresh'_${year_min_data}_${year_max_data}.dta"
		if _rc {
		 disp as error "USER ERROR: File ${DIR_TEMP}/akm_estimates_combined_firm_`inc_concept'_`thresh'_${year_min_data}_${year_max_data}.dta does not exist."
		 error 1
		}
		else {
			
			* start postfile
			postutil clear
			postfile postfile_comp_periods ///
				period_1_y1 period_1_y2 period_2_y1 period_2_y2 period_length ///
				b_level se_level r2_level b_level_period se_level_period r2_level_period r2_within_level_period rho_level rho_level_period N_level ///
				b_diff se_diff r2_diff b_diff_period se_diff_period r2_diff_period r2_within_diff_period rho_diff rho_diff_period N_diff ///
				using "${DIR_TEMP}/comp_periods.dta", replace
			
			* loop through period lengths
			foreach period_length in 2 4 8 16 {
				
				* prepare data with firm FEs from both periods
				local year_max_data_adj = ${year_max_data} - `period_length' + 1
				local period = 0
				foreach y1 in $year_min_data `year_max_data_adj' {
					local ++period
					local y2 = `y1' + `period_length' - 1
					cap confirm file "${DIR_TEMP}/akm_estimates_f_`inc_concept'_`thresh'_`y1'_`y2'.dta"
					if _rc {
						disp as error "USER ERROR: File ${DIR_TEMP}/akm_estimates_f_`inc_concept'_`thresh'_`y1'_`y2'.dta does not exist."
						error 1
					}
					else {
						use ${id_firm}_original ${year} ${f_fe} using "${DIR_TEMP}/akm_estimates_f_`inc_concept'_`thresh'_`y1'_`y2'.dta", clear
						${gtools}collapse (firstnm) ${f_fe}, by(${id_firm}_original ${year}) fast
						prog_comp_desc_sum_save "${DIR_TEMP}/temp_comp_`period'.dta"
					}
				}
				
				* load data with firm-year FEs from full period
				local period_1_y1 = $year_min_data
				local period_1_y2 = $year_min_data + `period_length' - 1
				local period_2_y1 = ${year_max_data} - `period_length' + 1
				local period_2_y2 = $year_max_data
				use ${id_firm}_original ${year} ${fy_fe} N if inrange(${year}, `period_1_y1', `period_1_y2') | inrange(${year}, `period_2_y1', `period_2_y2') using "${DIR_TEMP}/akm_estimates_combined_firm_`inc_concept'_`thresh'_${year_min_data}_${year_max_data}.dta", clear
				
				* merge in firm FEs from above
				gen float ${f_fe} = .
				label var ${f_fe} "Predicted AKM firm-year FE"
				gen byte period = .
				label var period "Period"
				forval period = 1/2 {
					merge 1:1 ${id_firm}_original ${year} using "${DIR_TEMP}/temp_comp_`period'.dta", update keepusing(${f_fe}) keep(master match_update) nogen
					rm "${DIR_TEMP}/temp_comp_`period'.dta"
					replace period = `period' if ${f_fe} < . & period == .
				}
				
				* keep only firms with both pay measures
				keep if ${f_fe} < . & ${fy_fe} < .
				
				* collapse to firm-period level
				foreach var of varlist * {
					local `var'_l: var label `var'
				}
				${gtools}collapse (firstnm) ${f_fe} (mean) ${fy_fe} (rawsum) N [aw=N], by(${id_firm}_original period) fast
				foreach var of varlist * {
					label var `var' "``var'_l'"
				}
				
				* take correlation of levels
				corr ${fy_fe} ${f_fe} [aw=N]
				local rho_level = r(rho)
				sum N if ${fy_fe} < . & ${f_fe} < .
				local N_level = r(sum)
				local vars_list = "${fy_fe} ${f_fe}"
				foreach var of varlist `vars_list' {
					if "${gtools}" == "" {
						gen float `var'_weighted = `var'*N
						bys period: egen long N_sum = total(N)
						bys period: egen float `var'_weighted_sum = total(`var'_weighted)
						gen float `var'_mean = `var'_weighted_sum/N_sum
						drop `var'_weighted `var'_weighted_sum N_sum
					}
					else gegen float `var'_mean = mean(`var') [aw=N], by(period)
					gen float `var'_demeaned = `var' - `var'_mean
					drop `var'_mean
				}
				corr ${fy_fe}_demeaned ${f_fe}_demeaned [aw=N]
				local rho_level_period = r(rho)
				drop ${fy_fe}_demeaned ${f_fe}_demeaned
				
				* regress firm-year FEs in levels on firm FEs in levels
				reghdfe ${fy_fe} ${f_fe} [aw=N], noabsorb vce(cluster ${id_firm}_original)
				local b_level = _b[${f_fe}]
				local se_level = _se[${f_fe}]
				local r2_level = e(r2)
				reghdfe ${fy_fe} ${f_fe} [aw=N], a(period) vce(cluster ${id_firm}_original)
				local b_level_period = _b[${f_fe}]
				local se_level_period = _se[${f_fe}]
				local r2_level_period = e(r2)
				local r2_within_level_period = e(r2_within)
				
				* create panel weights
				if "${gtools}" == "" bys ${id_firm}_original: egen long N_total = total(N)
				else gegen long N_total = total(N), by(${id_firm}_original)
				label var N_total "Number of observations across periods"
				drop N
				
				* set panel
				egen long ${id_firm} = group(${id_firm}_original)
				label var ${id_firm} "Firm ID (numeric)"
				drop ${id_firm}_original
				xtset ${id_firm} period
				
				* create differenced firm pay measures
				foreach var of varlist $f_fe $fy_fe {
					gen float D_`var' = D.`var'
					label var D_`var' "Differenced ``var'_l'"
				}
				keep if D_$f_fe < . & $fy_fe < .
				
				* take correlation of differences
				corr ${fy_fe} ${f_fe} [aw=N_total]
				local rho_diff = r(rho)
				sum N_total if ${fy_fe} < . & ${f_fe} < .
				local N_diff = r(sum)
				local vars_list = "${fy_fe} ${f_fe}"
				foreach var of varlist `vars_list' {
					if "${gtools}" == "" {
						gen float `var'_weighted = `var'*N_total
						bys period: egen long N_sum = total(N_total)
						bys period: egen float `var'_weighted_sum = total(`var'_weighted)
						gen float `var'_mean = `var'_weighted_sum/N_sum
						drop `var'_weighted `var'_weighted_sum N_sum
					}
					else gegen float `var'_mean = mean(`var') [aw=N_total], by(period)
					gen float `var'_demeaned = `var' - `var'_mean
					drop `var'_mean
				}
				corr ${fy_fe}_demeaned ${f_fe}_demeaned [aw=N_total]
				local rho_diff_period = r(rho)
				drop ${fy_fe}_demeaned ${f_fe}_demeaned
				
				* regress differenced firm-year FEs on differenced firm FEs
				reghdfe D_${fy_fe} D_${f_fe} [aw=N_total], noabsorb vce(cluster ${id_firm})
				local b_diff = _b[D_${f_fe}]
				local se_diff = _se[D_${f_fe}]
				local r2_diff = e(r2)
				reghdfe D_${fy_fe} D_${f_fe} [aw=N_total], a(period) vce(cluster ${id_firm})
				local b_diff_period = _b[D_${f_fe}]
				local se_diff_period = _se[D_${f_fe}]
				local r2_diff_period = e(r2)
				local r2_within_diff_period = e(r2_within)
				
				* write results to postfile
				post postfile_comp_periods ///
					(`period_1_y1') (`period_1_y2') (`period_2_y1') (`period_2_y2') (`period_length') ///
					(`b_level') (`se_level') (`r2_level') (`b_level_period') (`se_level_period') (`r2_level_period') (`r2_within_level_period') (`rho_level') (`rho_level_period') (`N_level') ///
					(`b_diff') (`se_diff') (`r2_diff') (`b_diff_period') (`se_diff_period') (`r2_diff_period') (`r2_within_diff_period') (`rho_diff') (`rho_diff_period') (`N_diff')
			}
			postclose postfile_comp_periods
			
			* format, save, and export postfile
			use "${DIR_TEMP}/comp_periods.dta", clear
			label var period_1_y1 "Period 1 start year"
			label var period_1_y2 "Period 1 end year"
			label var period_2_y1 "Period 2 start year"
			label var period_2_y2 "Period 2 end year"
			label var period_length "Period length"
			label var b_level "Point estimate in levels (reg fy_fe f_fe)"
			label var se_level "Standard error in levels (reg fy_fe f_fe)"
			label var r2_level "R^2 in levels (reg fy_fe f_fe)"
			label var b_level_period "Point estimate in levels w/ period FE (reg fy_fe f_fe)"
			label var se_level_period "Standard error in levels w/ period FE (reg fy_fe f_fe)"
			label var r2_level_period "R^2 in levels w/ period FE (reg fy_fe f_fe)"
			label var r2_within_level_period "Within R^2 in levels w/ period FE (reg fy_fe f_fe)"
			label var rho_level "Standard error in levels (reg fy_fe f_fe)"
			label var rho_level_period "Standard error in levels w/ period FEs (reg fy_fe f_fe)"
			label var N_level "Number of worker-years in levels"
			label var b_diff "Point estimate in differences (reg D_fy_fe D_f_fe)"
			label var se_diff "Standard error in differences (reg D_fy_fe D_f_fe)"
			label var r2_diff "R^2 in differences (reg D_fy_fe D_f_fe)"
			label var b_diff_period "Point estimate in differences w/ period FE (reg D_fy_fe D_f_fe)"
			label var se_diff_period "Standard error in differences w/ period FE (reg D_fy_fe D_f_fe)"
			label var r2_diff_period "R^2 in differences w/ period FE (reg D_fy_fe D_f_fe)"
			label var r2_within_diff_period "Within R^2 in differences w/ period FE (reg D_fy_fe D_f_fe)"
			label var rho_diff "Standard error in differences (reg D_fy_fe D_f_fe)"
			label var rho_diff_period "Standard error in differences w/ period FEs (reg D_fy_fe D_f_fe)"
			label var N_diff "Number of worker-years in differences"
			save "${DIR_TEMP}/comp_periods.dta", replace
			export delim using "${DIR_EXPORT}/csv/7_comp_periods.csv", replace
		}
		
		
		*** compare firm-year FEs over shorter time windows vs. firm-year FE model over full time window
		cap confirm file "${DIR_TEMP}/akm_estimates_combined_firm_`inc_concept'_`thresh'_${year_min_data}_${year_max_data}.dta"
		if _rc {
			disp as error "USER ERROR: File ${DIR_TEMP}/akm_estimates_combined_firm_`inc_concept'_`thresh'_${year_min_data}_${year_max_data}.dta does not exist."
			error 1
		}
		else {
			
			* start postfile
			postutil clear
			postfile postfile_comp_short_full ///
				b_level se_level r2_level b_level_year se_level_year r2_level_year r2_within_level_year rho_level rho_level_year N_level ///
				b_diff se_diff r2_diff b_diff_year se_diff_year r2_diff_year r2_within_diff_year rho_diff rho_diff_year N_diff ///
				b_fe se_fe r2_fe r2_within_fe b_fe_year se_fe_year r2_fe_year r2_within_fe_year N_fe ///
				using "${DIR_TEMP}/comp_short_full.dta", replace
			
			* prepare data with firm FEs from both periods
			local y1_1 = $year_min_data
			local y1_2 = $year_min_data + 8
			local y1_3 = $year_max_data - 15
			local y1_4 = $year_max_data - 7
			local period = 0
			foreach y1 in `y1_1' `y1_2' `y1_3' `y1_4' {
				local ++period
				disp "PERIOD = `period'"
				local y2 = `y1' + 7
				cap confirm file "${DIR_TEMP}/akm_estimates_f_`inc_concept'_`thresh'_`y1'_`y2'.dta"
				if _rc {
					disp as error "USER ERROR: File ${DIR_TEMP}/akm_estimates_f_`inc_concept'_`thresh'_`y1'_`y2'.dta does not exist."
					error 1
				}
				else {
					use ${id_firm}_original ${year} ${fy_fe} using "${DIR_TEMP}/akm_estimates_fy_`inc_concept'_`thresh'_`y1'_`y2'.dta", clear
					${gtools}collapse (firstnm) ${fy_fe}, by(${id_firm}_original ${year}) fast
					prog_comp_desc_sum_save "${DIR_TEMP}/temp_comp_`period'.dta"
				}
			}
			
			* load data with firm-year FEs from full period
			use ${id_firm}_original ${year} ${fy_fe} N using "${DIR_TEMP}/akm_estimates_combined_firm_`inc_concept'_`thresh'_${year_min_data}_${year_max_data}.dta", clear
			rename ${fy_fe} ${fy_fe}_full
			label var ${fy_fe}_full "Pred. AKM firm-year FE (full window)"
			
			* merge in firm-year FEs from above
			gen float ${fy_fe} = .
			gen byte period = .
			label var period "Period"
			forval period = 1/4 {
				merge 1:1 ${id_firm}_original ${year} using "${DIR_TEMP}/temp_comp_`period'.dta", update keepusing(${fy_fe}) keep(master match_update) nogen
				rm "${DIR_TEMP}/temp_comp_`period'.dta"
				replace period = `period' if ${fy_fe} < . & period == .
			}
			rename ${fy_fe} ${fy_fe}_short
			label var ${fy_fe}_short "Pred. AKM firm-year FE (short window)"
			
			* keep only firms with both pay measures
			keep if ${fy_fe}_short < . & ${fy_fe}_full < .
			
			* set panel
			egen long ${id_firm} = group(${id_firm}_original)
			label var ${id_firm} "Firm ID (numeric)"
			drop ${id_firm}_original
			xtset ${id_firm} ${year}
			
			* create panel weights
			gen long N_current_plus_lag = N + L.N
			label var N_current_plus_lag "Number of current + lagged observations"
			
			* create differenced firm pay measures
			local vars_list = "${fy_fe}_short ${fy_fe}_full"
			foreach var of varlist `vars_list' {
				gen float D_`var' = D.`var'
				local `var'_l: var label `var'
				label var D_`var' "Differenced ``var'_l'"
			}
			
			* regress firm-year FEs from short time window on firm-year FEs from full time window
			reghdfe ${fy_fe}_short ${fy_fe}_full [aw=N], noabsorb vce(cluster ${id_firm})
			local b_level = _b[${fy_fe}_full]
			local se_level = _se[${fy_fe}_full]
			local r2_level = e(r2)
			local N_level = e(N_full)
			reghdfe ${fy_fe}_short ${fy_fe}_full [aw=N], a(${year}) vce(cluster ${id_firm})
			local b_level_year = _b[${fy_fe}_full]
			local se_level_year = _se[${fy_fe}_full]
			local r2_level_year = e(r2)
			local r2_within_level_year = e(r2_within)

			* take correlation of levels
			corr ${fy_fe}_short ${fy_fe}_full [aw=N]
			local rho_level = r(rho)
			local vars_list = "${fy_fe}_short ${fy_fe}_full"
			foreach var of varlist `vars_list' {
				if "${gtools}" == "" {
					gen float `var'_weighted = `var'*N
					bys ${year}: egen long N_sum = total(N)
					bys ${year}: egen float `var'_weighted_sum = total(`var'_weighted)
					gen float `var'_mean = `var'_weighted_sum/N_sum
					drop `var'_weighted `var'_weighted_sum N_sum
				}
				else gegen float `var'_mean = mean(`var') [aw=N], by(${year})
				gen float `var'_demeaned = `var' - `var'_mean
				drop `var'_mean
			}
			corr ${fy_fe}_short_demeaned ${fy_fe}_full_demeaned [aw=N]
			local rho_level_year = r(rho)
			drop ${fy_fe}_short_demeaned ${fy_fe}_full_demeaned
			
			* regress differenced firm-year FEs from short time window on differenced firm-year FEs from full time window
			reghdfe D_${fy_fe}_short D_${fy_fe}_full [aw=N_current_plus_lag], noabsorb vce(cluster ${id_firm})
			local b_diff = _b[D_${fy_fe}_full]
			local se_diff = _se[D_${fy_fe}_full]
			local r2_diff = e(r2)
			local N_diff = e(N_full)
			reghdfe D_${fy_fe}_short D_${fy_fe}_full [aw=N_current_plus_lag], a(${year}) vce(cluster ${id_firm})
			local b_diff_year = _b[D_${fy_fe}_full]
			local se_diff_year = _se[D_${fy_fe}_full]
			local r2_diff_year = e(r2)
			local r2_within_diff_year = e(r2_within)
			
			* take correlation of levels
			corr D_${fy_fe}_short D_${fy_fe}_full [aw=N_current_plus_lag]
			local rho_diff = r(rho)
			local vars_list = "D_${fy_fe}_short D_${fy_fe}_full"
			foreach var of varlist `vars_list' {
				if "${gtools}" == "" {
					gen float `var'_weighted = `var'*N_current_plus_lag
					bys ${year}: egen long N_sum = total(N_current_plus_lag)
					bys ${year}: egen float `var'_weighted_sum = total(`var'_weighted)
					gen float `var'_mean = `var'_weighted_sum/N_sum
					drop `var'_weighted `var'_weighted_sum N_sum
				}
				else gegen float `var'_mean = mean(`var') [aw=N_current_plus_lag], by(${year})
				gen float `var'_demeaned = `var' - `var'_mean
				drop `var'_mean
			}
			corr D_${fy_fe}_short_demeaned D_${fy_fe}_full_demeaned [aw=N_current_plus_lag]
			local rho_diff_year = r(rho)
			drop D_${fy_fe}_short_demeaned D_${fy_fe}_full_demeaned
			
			* regress firm-year FEs from short time window on firm-year FEs from full time window, with firm FE
			reghdfe ${fy_fe}_short ${fy_fe}_full [aw=N], a(${id_firm}) vce(cluster ${id_firm})
			local b_fe = _b[${fy_fe}_full]
			local se_fe = _se[${fy_fe}_full]
			local r2_fe = e(r2)
			local r2_within_fe = e(r2_within)
			local N_fe = e(N_full)
			reghdfe ${fy_fe}_short ${fy_fe}_full [aw=N], a(${id_firm} ${year}) vce(cluster ${id_firm})
			local b_fe_year = _b[${fy_fe}_full]
			local se_fe_year = _se[${fy_fe}_full]
			local r2_fe_year = e(r2)
			local r2_within_fe_year = e(r2_within)
			
			* write results to postfile
			post postfile_comp_short_full ///
				(`b_level') (`se_level') (`r2_level') (`b_level_year') (`se_level_year') (`r2_level_year') (`r2_within_level_year') (`rho_level') (`rho_level_year') (`N_level') ///
				(`b_diff') (`se_diff') (`r2_diff') (`b_diff_year') (`se_diff_year') (`r2_diff_year') (`r2_within_diff_year') (`rho_diff') (`rho_diff_year') (`N_diff') ///
				(`b_fe') (`se_fe') (`r2_fe') (`r2_within_fe') (`b_fe_year') (`se_fe_year') (`r2_fe_year') (`r2_within_fe_year') (`N_fe')
			postclose postfile_comp_short_full
			
			* format, save, and export postfile
			use "${DIR_TEMP}/comp_short_full.dta", clear
			label var b_level "Point estimate in levels (reg fy_fe_short fy_fe_full)"
			label var se_level "Standard error in levels (reg fy_fe_short fy_fe_full)"
			label var r2_level "R^2 in levels (reg fy_fe_short fy_fe_full)"
			label var b_level_year "Point estimate in levels w/ year FEs (reg fy_fe_short fy_fe_full)"
			label var se_level_year "Standard error in levels w/ year FEs (reg fy_fe_short fy_fe_full)"
			label var r2_level_year "R^2 in levels w/ year FEs (reg fy_fe_short fy_fe_full)"
			label var r2_within_level_year "Within R^2 in levels w/ year FEs (reg fy_fe_short fy_fe_full)"
			label var rho_level "Correlation in levels (corr fy_fe_short fy_fe_full)"
			label var rho_level_year "Correlation in levels after yearly demeaning (corr fy_fe_short fy_fe_full)"
			label var N_level "Number of worker-years in levels (reg fy_fe_short fy_fe_full)"
			label var b_diff "Point estimate in differences (reg D_fy_fe_short D_fy_fe_full)"
			label var se_diff "Standard error in differences (reg D_fy_fe_short D_fy_fe_full)"
			label var r2_diff "R^2 in differences (reg D_fy_fe_short D_fy_fe_full)"
			label var b_diff_year "Point estimate in differences w/ year FEs (reg D_fy_fe_short D_fy_fe_full)"
			label var se_diff_year "Standard error in differences w/ year FEs (reg D_fy_fe_short D_fy_fe_full)"
			label var r2_diff_year "R^2 in differences w/ year FEs (reg D_fy_fe_short D_fy_fe_full)"
			label var r2_within_diff_year "Within R^2 in differences w/ year FEs (reg D_fy_fe_short D_fy_fe_full)"
			label var rho_diff "Correlation in differences (corr D_fy_fe_short D_fy_fe_full)"
			label var rho_diff_year "Correlation in differences after yearly demeaning (corr D_fy_fe_short D_fy_fe_full)"
			label var N_diff "Number of worker-years in differences (reg D_fy_fe_short D_fy_fe_full)"
			label var b_fe "Point estimate with firm FE (reg fy_fe_short fy_fe_full, FE)"
			label var se_fe "Standard error with firm FE (reg fy_fe_short fy_fe_full, FE)"
			label var r2_fe "R^2 with firm FE (reg fy_fe_short fy_fe_full, FE)"
			label var r2_within_fe "Within R^2 with firm FE (reg fy_fe_short fy_fe_full, FE)"
			label var b_fe_year "Point estimate with firm FE w/ year FEs (reg fy_fe_short fy_fe_full, FE)"
			label var se_fe_year "Standard error with firm FE w/ year FEs (reg fy_fe_short fy_fe_full, FE)"
			label var r2_fe_year "R^2 with firm FE w/ year FEs (reg fy_fe_short fy_fe_full, FE)"
			label var r2_within_fe_year "Within R^2 with firm FE w/ year FEs (reg fy_fe_short fy_fe_full, FE)"
			label var N_fe "Number of worker-years with firm FE (reg fy_fe_short fy_fe_full, FE)"
			save "${DIR_TEMP}/comp_short_full.dta", replace
			export delim using "${DIR_EXPORT}/csv/7_comp_short_full.csv", replace
		}
	}
}
