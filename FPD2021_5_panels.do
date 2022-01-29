********************************************************************************
* DESCRIPTION: Creates various panel datasets based on AKM estimates and firm
*              financials.
*
* CREDIT:      Please cite: Engbom, Niklas & Christian Moser & Jan Sauermann,
*              2021. "Firm Pay Dynamics," Working Paper.
*
* TIME STAMP:  October 23, 2021.
********************************************************************************


*** prepare firm financials data
* load data
local vars_list_fundamental = "${f_ind} ${f_size_fek}"
local vars_list_financial = "${f_rev} ${f_va} ${f_assets}"
if $test {
	use ${id_firm} ${year} `vars_list_fundamental' `vars_list_financial' using "${DIR_TEMP}/test_data_firm.dta", clear
	rename ${id_firm} ${id_firm}_original
}
else if "${user}" == "SWE" {
	use ${id_firm} ${year} `vars_list_fundamental' `vars_list_financial' using "${RawDataDir}/CombinedFirm.dta", clear
	rename ${id_firm} ${id_firm}_original
	label var ${id_firm}_original "Firm ID (string, original)"
	keep if !missing(${id_firm}_original)
	cap confirm str var ${id_firm}_original
	if !_rc keep if ${id_firm}_original != "00000000"
	${gtools}collapse (firstnm) ${f_ind} ${f_size_fek} ${f_rev} ${f_va} ${f_assets}, by(${id_firm}_original ${year})
}
else {
	disp as error `"USER ERROR: Cannot load firm financials data if \$test == 0 and "\${user}" != "SWE"."'
}

* label variables
label var ${f_ind} "Sector"
label var ${f_size_fek} "Firm size from FEK (number of employees)"
label var ${f_rev} "Sales (SEK)"
label var ${f_va} "Value added (SEK)"
// label var ${f_ebit} "EBIT (SEK)"
// label var ${f_profit} "Profits (SEK)"
label var ${f_assets} "Total assets (SEK)"
// label var ${f_assets_lt} "Long-term assets (SEK)"
// label var ${f_equip} "Equipment (SEK)"
// label var ${f_struct} "Structures (SEK)"
// label var ${f_cash} "Cash (SEK)"
// label var ${f_debt} "Total debt (SEK)"
// label var ${f_equity} "Total equity (SEK)"
// label var ${f_inv} "Investment (SEK)"

* ensure that firm size is >= 1
replace ${f_size_fek} = 1 if ${f_size_fek} < 1

* take logarithm of firm size and all financial variables
foreach var of varlist $f_size_fek `vars_list_financial' {
	local l_`var': var label `var'
	gen float ln_`var' = ln(`var')
	local l_`var'_sub = subinstr("`l_`var''", "(SEK", "(log SEK", .)
	local l_`var'_sub = subinstr("`l_`var'_sub'", "(number", "(log number", .)
	label var ln_`var' "`l_`var'_sub'"
	drop `var'
	rename ln_`var' `var'
}

* create log per-worker financial variables
local vars_list_financial_pw = ""
foreach var of varlist `vars_list_financial' {
	local l_`var': var label `var'
	gen float `var'_pw = `var' - ${f_size_fek}
	local l_`var'_sub = subinstr("`l_`var''", "(log SEK", "per worker (log SEK", .)
	label var `var'_pw "`l_`var'_sub'"
	local vars_list_financial_pw = "`vars_list_financial_pw' `var'_pw"
}

* save processed firm financials data
order ${id_firm}_original ${year} `vars_list_fundamental' `vars_list_financial' `vars_list_financial_pw'
sort ${id_firm}_original ${year}
prog_comp_desc_sum_save "${DIR_TEMP}/firm_financials.dta"


*** combine estimates from both AKM models: firm FEs and firm-year FEs
* loop through different individual income concepts, minimum firm size thresholds, and types of AKM model (firm FE and firm-year FE)
foreach inc_concept of global inc_concept_list {
	foreach thresh of global AKM_threshold_list {
		disp _newline(5)
		disp "--> inc_concept = `inc_concept', thresh = `thresh'"
		
		* load worker-level data containing firm FEs
		use ${id_firm} ${f_fe} using "${DIR_TEMP}/akm_estimates_f_`inc_concept'_`thresh'_${year_min}_${year_max}.dta", clear // list of all variables: ${id_pers} ${year} ${id_firm} ${gender} ${edu} ${age} `inc_concept' `inc_concept'${demeaned} pe ${f_fe} `var_xb_year' `var_xb_age' resid
		
		* collapse to firm level
		${gtools}collapse (firstnm) ${f_fe}, by(${id_firm}) fast
		recast float ${f_fe}, force
		label var ${f_fe} "Predicted AKM firm FE"
		
		* save temp data
		save "${DIR_TEMP}/temp_f_collapsed.dta", replace
		
		* load data with firm FE estimates
		use "${DIR_TEMP}/akm_estimates_fy_`inc_concept'_`thresh'_${year_min}_${year_max}.dta", clear
		
		* drop variables that are no longer needed
		drop pe resid
		cap drop *xb*
		
		* merge in data with firm FE estimates
		merge m:1 ${id_firm} using "${DIR_TEMP}/temp_f_collapsed.dta", keep(match) keepusing(${f_fe}) nogen
		rm "${DIR_TEMP}/temp_f_collapsed.dta"
		
		* generate firm year of incorporation (proxied by year of first observation)
		bys ${id_firm} (${year}): gen int ${f_yob} = ${year}[1]
		label var ${f_yob} "Firm year of incorporation"

		* generate firm age
		gen byte ${f_age} = ${year} - ${f_yob}
		label var ${f_age} "Firm age (years)"
		
		* generate firm size
		bys ${id_firm} ${year}: gen float ${f_size} = ln(_N)
		label var ${f_size} "Firm size from RAMS (log number of employees)"
		
		* save merged file
		sort ${id_pers} ${year}
		prog_comp_desc_sum_save "${DIR_TEMP}/akm_estimates_combined_`inc_concept'_`thresh'_${year_min}_${year_max}.dta"
	}
}


*** collapse combined data to the firm-year level
* loop through different individual income concepts and minimum firm size thresholds
foreach inc_concept of global inc_concept_list {
	foreach thresh of global AKM_threshold_list {
		disp _newline(5)
		disp "--> inc_concept = `inc_concept', thresh = `thresh'"
		
		* load data
		use ${id_firm} ${id_firm}_original ${year} `inc_concept' `inc_concept'${demeaned} ${f_fe} ${fy_fe} ${f_yob} ${f_age} ${f_size} using "${DIR_TEMP}/akm_estimates_combined_`inc_concept'_`thresh'_${year_min}_${year_max}.dta", clear
		
		* collapse to firm-year level
		foreach var of varlist * {
			local l_`var': var label `var'
		}
		${gtools}collapse (mean) `inc_concept'* (firstnm) ${id_firm} ${f_fe} ${fy_fe} ${f_yob} ${f_age} ${f_size} (count) N=${f_fe}, by(${id_firm}_original ${year}) fast
		recast float `inc_concept'* ${f_fe} ${fy_fe}, force
		label var ${id_firm} "Firm ID (numeric)"
		label var `inc_concept' "Mean `l_`inc_concept''"
		if "${demeaned}" != "" label var `inc_concept'${demeaned} "Mean `l_`inc_concept'${demeaned}'"
		foreach var of varlist $f_fe $fy_fe $f_yob $f_age $f_size {
			label var `var' "`l_`var''"
		}
		label var N "Weight (number of employees in RAMS)"
		
		* save collapsed file
		order ${id_firm} ${id_firm}_original ${year} `inc_concept' `inc_concept'${demeaned} ${f_fe} ${fy_fe} ${f_yob} ${f_age} ${f_size} N
		sort ${id_firm} ${year}
		prog_comp_desc_sum_save "${DIR_TEMP}/akm_estimates_combined_firm_`inc_concept'_`thresh'_${year_min}_${year_max}.dta"
	}
}


*** combine estimates from AKM models with firm financials
* loop through different individual income concepts and minimum firm size thresholds
foreach inc_concept of global inc_concept_list {
	foreach thresh of global AKM_threshold_list {
		disp _newline(5)
		disp "--> inc_concept = `inc_concept', thresh = `thresh'"
		
		* load data with merged AKM model results
		use "${DIR_TEMP}/akm_estimates_combined_firm_`inc_concept'_`thresh'_${year_min}_${year_max}.dta", clear
		
		* merge in processed firm financials data
		merge 1:1 ${id_firm}_original ${year} using "${DIR_TEMP}/firm_financials.dta", keep(match) keepusing(`vars_list_fundamental' `vars_list_financial' `vars_list_financial_pw') nogen
		
		* save merged file with firm financials data
		order ${id_firm} ${id_firm}_original ${year} `inc_concept' `inc_concept'${demeaned} ${f_fe} ${fy_fe} `vars_list_fundamental' ${f_yob} ${f_age} ${f_size} `vars_list_financial' `vars_list_financial_pw' N
		sort ${id_firm} ${year}
		prog_comp_desc_sum_save "${DIR_TEMP}/akm_estimates_combined_firm_financials_`inc_concept'_`thresh'_${year_min}_${year_max}.dta"
	}
}
