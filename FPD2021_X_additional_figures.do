********************************************************************************
* additional figures for EMS (2021)
*
* Date: October 23, 2021
********************************************************************************


*** macros
global DIR_MAIN_ADDITIONAL_FIGURES = "/Users/cm3594/Dropbox (CBS)/AKM Conference"

global DIR_MAIN_ADDITIONAL_FIGURES = "/Users/cm3594/Dropbox (CBS)/AKM Conference"

// *** variance of (demeaned) log earnings and normalized percentiles of log earnings
// * load data
// import delim using "${DIR_MAIN_ADDITIONAL_FIGURES}/5_results/7_IFAU_20211016/csv/8_time_series_var_perc.csv", clear

// * define macros
// local inc_concept = "earnings"
// global year = "year"
// // local y_label = ".15(.05).4"
// local y_label = "0(.1).4"
// local norm_year = 1985

// * generate empty variable for plotting purposes
// gen byte empty = .

// * time series of var(demeaned log earnings)
// tw ///
// 	(connected `inc_concept'_norm_var ${year}, sort lcolor(blue none) lpattern(l blank) mcolor(blue none) msymbol(Oh none)) ///
// 	, graphregion(color(white)) plotregion(lcolor(black) margin(l=0 r=0 b=0 t=0)) ///
// 	xlabel(1985(5)2015, gmin gmax grid gstyle(dot)) ylabel(`y_label', gmin gmax grid gstyle(dot)) ///
// 	xtitle("") ytitle("Variance of log earnings") ///
// 	legend(region(lcolor(white))) ///
// 	name(ts_`inc_concept'_norm_var, replace)
// graph export "${DIR_MAIN_ADDITIONAL_FIGURES}/8_paper/1_inputs/8_time_series_`inc_concept'_var.eps", replace

// * comparison of time series of var(log earnings) and var(demeaned log earnings)
// // tw ///
// // 	(connected `inc_concept'_var `inc_concept'_norm_var ${year}, sort lcolor(blue red) lpattern(l _ -) mcolor(blue red) msymbol(Oh Dh)) ///
// // 	, graphregion(color(white)) plotregion(lcolor(black) margin(l=0 r=0 b=0 t=0)) ///
// // 	xlabel(, gmin gmax grid gstyle(dot)) ylabel(`y_label', gmin gmax grid gstyle(dot)) ///
// // 	xtitle("") ytitle("Variance of log earnings") ///
// // 	legend(order(1 "Earnings" 2 "Demeaned earnings") cols(2) region(lcolor(white))) ///
// // 	name(ts_`inc_concept'_var_comparison, replace)

// * time series of percentiles of log earnings
// tw ///
// 	(connected `inc_concept'_p5_norm `inc_concept'_p10_norm `inc_concept'_p25_norm `inc_concept'_p50_norm `inc_concept'_p75_norm `inc_concept'_p90_norm `inc_concept'_p95_norm ${year}, sort lcolor(blue red green orange cyan magenta lime gold purple) lpattern(solid longdash dash shortdash dot longdash_dot dash_dot shortdash_dot -.) mcolor(blue red green orange cyan magenta lime gold purple) msymbol(O D T S Oh Dh Th Sh X)) ///
// 	, graphregion(color(white)) plotregion(lcolor(black) margin(l=0 r=0 b=0 t=0)) ///
// 	xlabel(1985(5)2015, gmin gmax grid gstyle(dot)) ylabel(-.2(.2).8, gmin gmax grid gstyle(dot)) ///
// 	xtitle("") ytitle("Normalized percentiles of log earnings (`norm_year' = 0.0)") ///
// 	legend(order(1 "P5" 2 "P10" 3 "P25" 4 "P50" 5 "P75" 6 "P90" 7 "P95") cols(4) region(lcolor(white)))  ///
// 	name(ts_perc_norm_`inc_concept', replace)
// graph export "${DIR_MAIN_ADDITIONAL_FIGURES}/8_paper/1_inputs/8_time_series_perc_norm_`inc_concept'.eps", replace


*** comparison of firm pay measures
* define macros
local inc_concept = "earnings"

* loop through firm pay measures
*foreach var_x in `inc_concept'_norm $f_fe $fy_fe {
foreach var_x in `inc_concept'_demeaned $f_fe $fy_fe {
	
	* load data
	import delim using "${DIR_EXPORT}/csv/7_comp_var_x_`var_x'.csv", clear
	*import delim using "${DIR_MAIN_ADDITIONAL_FIGURES}/5_results/7_IFAU_20211016/csv/7_comp_var_x_`var_x'.csv", clear
	
	* prepare plots
	local dep_var = subinstr("`inc_concept'_norm ${f_fe} ${fy_fe}", "`var_x'", "", .)
	local counter = 1
	disp "dep_var = `dep_var'"
	foreach var2 of local dep_var {
		local var_`counter' = "`var2'"
		if "`var2'" == "`inc_concept'_norm" {
			local legend_`counter' = "Mean earnings"
			local legend_`counter'_small = "mean earnings"
		}
		else if "`var2'" == "${f_fe}" {
			local legend_`counter' = "Firm FEs"
			local legend_`counter'_small = "firm FEs"
		}
		else if "`var2'" == "${fy_fe}" {
			local legend_`counter' = "Firm-year FEs"
			local legend_`counter'_small = "firm-year FEs"
		}
		local ++counter
	}
	if "`var_x'" == "`inc_concept'_norm" local indep_var_title = "Firm-level mean earnings"
	else if "`var_x'" == "${f_fe}" local indep_var_title = "Firm FEs"
	else if "`var_x'" == "${fy_fe}" local indep_var_title = "Firm-year FEs"
	
	* plot mean of variables on y-axis vs. mean of variable on x-axis, without regression estimates
// 	tw ///
// 		(scatter `var_1'_mean `var_x'_mean, mcolor(blue) msymbol(O)) ///
// 		(lfit `var_1'_mean `var_x'_mean [aw=n], lcolor(blue)) ///
// 		(scatter `var_2'_mean `var_x'_mean, mcolor(red) msymbol(Oh)) ///
// 		(lfit `var_2'_mean `var_x'_mean [aw=n], lcolor(red)) ///
// 		(function y = x, range(-1 1) lcolor(black)) ///
// 		, graphregion(color(white)) plotregion(lcolor(black) margin(l=0 r=0 b=0 t=0)) ///
// 		xlabel(-1(.25)1, format(%3.2f) gmin gmax grid gstyle(dot)) ylabel(-1(.25)1, format(%3.2f) gmin gmax grid gstyle(dot)) ///
// 		xtitle("`indep_var_title'") ytitle("Mean firm pay measure") ///
// 		legend(order(1 "`legend_1'" 2 "Best fit for `legend_1_small'" 3 "`legend_2'" 4 "Best fit for `legend_2_small'") cols(2) region(lcolor(white))) ///
// 		name(comp_`var_x'_mean, replace)
// 	graph export "${DIR_MAIN_ADDITIONAL_FIGURES}/8_paper/1_inputs/7_comp_`var_x'_mean.eps", replace
	
	* plot mean of variables on y-axis vs. mean of variable on x-axis, with regression estimates
	if "`var_x'" == "`inc_concept'_norm" {
		reghdfe `var_1'_mean `var_x'_mean, noabsorb
		local b_1: di %4.3f `=_b[`var_x'_mean]'
		local se_1: di %4.3f `=_se[`var_x'_mean]'
		reghdfe `var_2'_mean `var_x'_mean, noabsorb
		local b_2: di %4.3f `=_b[`var_x'_mean]'
		local se_2: di %4.3f `=_se[`var_x'_mean]'
		tw ///
			(scatter `var_1'_mean `var_x'_mean, mcolor(blue) msymbol(O)) ///
			(lfit `var_1'_mean `var_x'_mean [aw=n], lcolor(blue)) ///
			(scatter `var_2'_mean `var_x'_mean, mcolor(red) msymbol(Oh)) ///
			(lfit `var_2'_mean `var_x'_mean [aw=n], lcolor(red)) ///
			(function y = x, range(-1 1) lcolor(black)) ///
			, text(-.3 -.03 "`legend_1' slope (s.e.) = `b_1' (`se_1')", color(blue) place(se)) ///
			text(-.4 -.03 "`legend_2' slope (s.e.) = `b_2' (`se_2')", color(red) place(se)) ///
			graphregion(color(white)) plotregion(lcolor(black) margin(l=0 r=0 b=0 t=0)) ///
			xlabel(-1(.25)1, format(%3.2f) gmin gmax grid gstyle(dot)) ylabel(-1(.25)1, format(%3.2f) gmin gmax grid gstyle(dot)) ///
			xtitle("`indep_var_title'") ytitle("Mean firm pay measure") ///
			legend(order(1 "`legend_1'" 2 "Best fit for `legend_1_small'" 3 "`legend_2'" 4 "Best fit for `legend_2_small'") cols(2) region(lcolor(white))) ///
			name(comp_`var_x'_mean_reg, replace)
		graph export "${DIR_EXPORT}/eps/7_comp_`var_x'_mean_reg.eps", replace
		*graph export "${DIR_MAIN_ADDITIONAL_FIGURES}/8_paper/1_inputs/7_comp_`var_x'_mean_reg.eps", replace
	}
	
	* plot mean of variables on y-axis vs. mean of variable on x-axis, with regression estimates (1)
	if "`var_x'" == "${f_fe}" {
		reghdfe ${fy_fe}_mean `var_x'_mean, noabsorb
		local b: di %4.3f `=_b[`var_x'_mean]'
		local se: di %4.3f `=_se[`var_x'_mean]'
		tw ///
			(function y = x, range(-.75 .75) lcolor(black)) ///
			(scatter ${fy_fe}_mean `var_x'_mean, mcolor(blue) msymbol(Oh)) ///
			(lfit ${fy_fe}_mean `var_x'_mean [aw=n], lcolor(blue)) ///
			, text(-.4 -.08 "Firm-year FEs slope (s.e.) = `b' (`se')", color(blue) place(se)) ///
			graphregion(color(white)) plotregion(lcolor(black) margin(l=0 r=0 b=0 t=0)) ///
			xlabel(-.75(.25).75, format(%3.2f) gmin gmax grid gstyle(dot)) ylabel(-.75(.25).75, format(%3.2f) gmin gmax grid gstyle(dot)) ///
			xtitle("`indep_var_title'") ytitle("Mean firm-year FE") ///
			legend(order(2 "Firm-year FEs" 3 "Best fit") cols(2) region(lcolor(white))) ///
			name(comp_`var_x'_m_r_short_1, replace)
		graph export "${DIR_EXPORT}/eps/7_comp_`var_x'_mean_reg_short_1.eps", replace
		*graph export "${DIR_MAIN_ADDITIONAL_FIGURES}/8_paper/1_inputs/7_comp_`var_x'_mean_reg_short_1.eps", replace
	}
	
	* plot mean of variables on y-axis vs. mean of variable on x-axis, with regression estimates (1)
	if "`var_x'" == "${fy_fe}" {
		reghdfe ${f_fe}_mean `var_x'_mean, noabsorb
		local b: di %4.3f `=_b[`var_x'_mean]'
		local se: di %4.3f `=_se[`var_x'_mean]'
		tw ///
			(function y = x, range(-.75 .75) lcolor(black)) ///
			(scatter ${f_fe}_mean `var_x'_mean, mcolor(red) msymbol(Dh)) ///
			(lfit ${f_fe}_mean `var_x'_mean [aw=n], lcolor(red)) ///
			, text(-.4 -.08 "Firm FEs slope (s.e.) = `b' (`se')", color(red) place(se)) ///
			graphregion(color(white)) plotregion(lcolor(black) margin(l=0 r=0 b=0 t=0)) ///
			xlabel(-.75(.25).75, format(%3.2f) gmin gmax grid gstyle(dot)) ylabel(-.75(.25).75, format(%3.2f) gmin gmax grid gstyle(dot)) ///
			xtitle("`indep_var_title'") ytitle("Mean firm FE") ///
			legend(order(2 "Firm FEs" 3 "Best fit") cols(2) region(lcolor(white))) ///
			name(comp_`var_x'_m_r_short_2, replace)
		graph export "${DIR_EXPORT}/eps/7_comp_`var_x'_mean_reg_short_2.eps", replace
		*graph export "${DIR_MAIN_ADDITIONAL_FIGURES}/8_paper/1_inputs/7_comp_`var_x'_mean_reg_short_2.eps", replace
	}
	
	* plot percentiles of variables on y-axis vs. quantiles of variable on x-axis
// 	tw ///
// 		(scatter `var_1'_p10 `var_1'_p50 `var_1'_p90 `var_x'_mean, mcolor(blue midblue eltblue) msymbol(O D T)) ///
// 		(scatter `var_2'_p10 `var_2'_p50 `var_2'_p90 `var_x'_mean, mcolor(red orange_red orange) msymbol(Oh Dh Th)) ///
// 		(function y = x, range(-1 1) lcolor(black)) ///
// 		, graphregion(color(white)) plotregion(lcolor(black) margin(l=0 r=0 b=0 t=0)) ///
// 		xlabel(-1(.25)1, format(%3.2f) gmin gmax grid gstyle(dot)) ylabel(-1(.25)1, format(%3.2f) gmin gmax grid gstyle(dot)) ///
// 		xtitle("`indep_var_title'") ytitle("Quantiles of firm pay measure") ///
// 		legend(order(1 "P10 `legend_1_small'" 2 "P50 `legend_1_small'" 3 "P90 `legend_1_small'" 4 "P10 `legend_2_small'" 5 "P50 `legend_2_small'" 6 "P90 `legend_2_small'") cols(3) region(lcolor(white))) ///
// 		name(comp_`var_x'_perc, replace)
// 	graph export "${DIR_MAIN_ADDITIONAL_FIGURES}/8_paper/1_inputs/7_comp_`var_x'_perc.eps", replace
}
