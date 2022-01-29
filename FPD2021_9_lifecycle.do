********************************************************************************
* DESCRIPTION: Investigates firm life cycle profiles.
*
* CREDIT:      Please cite: Engbom, Niklas & Christian Moser & Jan Sauermann,
*              2021. "Firm Pay Dynamics," Working Paper.
*
* TIME STAMP:  October 23, 2021.
********************************************************************************
* graphing details
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


* Firm life-cycle outcomes
global firms = 0

* Worker life-cycle outcomes
global workers = 0

* graph output
global graphs = 1


*** check estimated age profiles, especially their flat parts
foreach inc_concept of global inc_concept_list {
	foreach thresh of global AKM_threshold_list {
	
		* load data
		use "${DIR_TEMP}/akm_estimates_fy_`inc_concept'_`thresh'_${year_min}_${year_max}.dta", clear
		
		* detect interactions
		local vars_inter = ""
		if $akm_gender_inter {
			local vars_inter = "${gender}"
// 			if "${user}" == "SWE" {
// 				label define gen_l 1 "Male" 2 "Female", replace
// 				label val ${gender} gen_l
// 			}
		}
		if $akm_edu_inter {
			local vars_inter = "`vars_inter' ${edu}"
// 			if "${user}" == "SWE" {
// 				label define edu_l 2 "Less than High School" 3 "High School" 4 "Some College" 5 "College" 6 "Postgraduate", replace
// 				label val ${edu} edu_l
// 			}
		}
		
		* collapse estimated age effects to relevant subgroup and age
		${gtools}collapse (firstnm) xb_age, by(`vars_inter' ${age}) fast
		
		* create demographic (gender x education) variable
		if "`vars_inter'" == "" local vars_by = ""
		else local vars_by = "by(`vars_inter')"
		
		* export CSV file
		export delim using "${DIR_EXPORT}/csv/9_lifecycle_age_profiles.csv", replace
		
		* produce graphs
		tw connected xb_age age, sort `vars_by' name("age_profiles_joint", replace)
		graph_export_eps_pdf "9_lifecycle_age_profiles_joint"
	}
}


*** loop through different individual income concepts and minimum firm size thresholds
foreach inc_concept of global inc_concept_list {

	foreach thresh of global AKM_threshold_list {
		
		
		*** FIRM LIFECYCLE DYNAMICS
		if $firms {
		
			* make sure all export directories exist
// 			foreach folder in csv eps out pdf tex {
// 				local dir_file = "${DIR_EXPORT}/`folder'"
// 				${cmd_confirm_dir} "`dir_file'"
// 				if _rc {
// 					!mkdir "`dir_file'"
// 				}
// 			}
			
			use "${DIR_TEMP}/akm_estimates_combined_firm_financials_`inc_concept'_`thresh'_${year_min}_${year_max}.dta" if ${f_age} >= 0 & yob_firm > 1985, clear
			
			local moment_list = "mean Var"
			local varlist = "${f_va}_pw $f_size earnings_norm fe fye dfye"
			xtset ${id_firm} ${year}
			gen dfye = fye - l.fye
			qui sum yob_firm
			local y1 = r(min)
			local y2 = r(max)
			qui gen balanced = .
			forvalues yy = `y1' / `y2' {
				qui sum ${f_age} if yob_firm == `yy'
				local r1 = r(min)
				local r2 = r(max)
				local d = `r2'-`r1'+1
				bys ${id_firm}: gen temp = ( _N == `d' )
				replace balanced = temp if yob_firm == `yy'
				drop temp
			}
			if "${gtools}" == "" bys ${id_firm}: egen float weight = mean(N)
			else gegen float weight = mean(N), by(${id_firm})
			
			foreach moment of local moment_list {
				global vars_list_`moment' = ""
				foreach var of local varlist {	
					foreach version in "" "_balanced" "_weight" {
						global vars_list_`moment' = "${vars_list_`moment'} `var'`version'_`moment'"
					}
				}
			}
			
			foreach var of local varlist {
				foreach version in "" "_balanced" "_weight" {
					foreach moment of local moment_list {
						qui gen `var'`version'_`moment' = .
					}
					levelsof yob_firm, local(list1)
					foreach yy in `list1' {
						levelsof $f_age if yob_firm == `yy', local(list)
						foreach aa in `list' {	
							if "`version'" == "" local cond = "if ${f_age} == `aa' & yob_firm == `yy' [aw=N]"
							else if "`version'" == "_balanced" local cond = "if balanced & ${f_age} == `aa' & yob_firm == `yy' [aw=N]"
							else if "`version'" == "_weight" local cond = "if balanced & ${f_age} == `aa' & yob_firm == `yy' [aw=weight]"
							qui sum `var' `cond'
							foreach moment of local moment_list {
								qui replace `var'`version'_`moment' = r(`moment') if ${f_age} == `aa' & yob_firm == `yy'
							}
						}
					}
				}
			}
			
			gcollapse (firstnm) ${vars_list_Var} ${vars_list_mean}, by( $f_age yob_firm ) fast
			
			foreach var in $vars_list_Var $vars_list_mean {
				format `var' %9.6g
			}
			
			* save as csv file
			export delim using "${DIR_EXPORT}/csv/9_lifecycle_firm.csv", replace
			
		}
		
		
		
		
		
		
		
		
		
		*** WORKER LIFECYCLE DYNAMICS
		if $workers {
			
			* make sure all export directories exist
			foreach folder in csv eps out pdf tex {
				local dir_file = "${DIR_EXPORT}/`folder'"
				${cmd_confirm_dir} "`dir_file'"
				if _rc {
					!mkdir "`dir_file'"
				}
			}
			local list = "${earn}_norm fe fye fye_init fye_post"
			
			use "${DIR_TEMP}/akm_estimates_fy_`inc_concept'_`thresh'_${year_min}_${year_max}.dta", clear
			
			
			recode educ (1 2 3 4=1)(5 6 7=2)
			gen yob = year-age
			
			merge 1:1 $id_pers $year using "${DIR_TEMP}/akm_estimates_f_`inc_concept'_`thresh'_${year_min}_${year_max}.dta", keepusing(fe) keep(match) nogen
			
			* generate firm component at the start of an individual's job spell
			bys $id_pers $id_firm ( $year ): gen fye_init = fye[1]
			gen fye_post = fye-fye_init
			
			local collapselist = ""
			foreach var in `list' {
				local collapselist = "`collapselist' (mean) `var'_mean = `var' (sd) `var'_var = `var'"
			}
			
			* collapse to age education level
			gcollapse `collapselist', by(${age} ${edu} ${gender} yob) fast

			* generate firm component at the start of an individual's job spell
			foreach var in `list' {
				format `var'_mean %9.6g
				replace `var'_var = `var'_var^2
				format `var'_var %9.6g
			}
			
			* save as csv file
			export delim using "${DIR_EXPORT}/csv/9_lifecycle_worker.csv", replace
						
		}		
		
		
		
		
		
		
		
		
		
		
		

		
		if $graphs {
			
			set graph on

			local listcolor = "blue red green orange cyan magenta lime gold purple"
			local listpattern = "solid longdash dash shortdash dot longdash_dot dash_dot shortdash_dot -."
			local listmsymbol = "O D T S Oh Dh Th Sh X"
			*JS global DIR_EXPORT = "${DIR_MAIN}/5_results/7_IFAU_20211016"

			* DEC2021 local moment_list = "mean var"
			local moment_list = "var"
			* DEC2021 local varlist = "${f_va}_pw $f_size earnings_norm fe fye dfye"
			local varlist = "fye"

			* load data
			insheet using "${DIR_EXPORT}/csv/9_lifecycle_firm.csv", clear
			keep if yob_firm > 1985 
			drop if f_age == .
			gen year = yob_firm + ${f_age}
			
			* normalize means to initial firm age
// 			foreach var of local varlist {
// 				foreach version in "" "_balanced" "_weight" {
// 					qui sum `var'`version'_mean if ${f_age} == 0 & yob_firm == 1986
// 					replace `var'`version'_mean = `var'`version'_mean-r(mean)
// 				}
// 			}
			
			foreach moment of local moment_list {
				if "`moment'" == "mean" local ytitle = "Mean"
				else if "`moment'" == "var" local ytitle = "Variance"
				else disp "moment not defined"
				foreach var of local varlist {
					local list = "1986 1989 1992 1995 1998 2001 2005 2008"
					local plotlist = ""
					local leglist = ""
					local i = 1 
					foreach yy in `list' {
						local col: word `i' of `listcolor'
						local mcol: word `i' of `listcolor'
						local pat: word `i' of `listpattern'
						local sym: word `i' of `listmsymbol'
						local plotlist = "`plotlist' (${LineStyle} `var'_`moment' ${f_age} if yob_firm == `yy', lcolor(`col') lpattern(`pat') mcolor(`col') msymbol(`sym'))"
						local leglist = `"`leglist' `i' "`yy'""'
						local i=`i'+1
						di "Plotlist"
						di "`plotlist'"
					}
					tw `plotlist', ///
						   graphregion(color(white)) plotregion(lcolor(black) margin(l=0 r=0 b=0 t=0)) ///
						   xlabel(, gmin gmax grid gstyle(dot)) ylabel(, gmin gmax grid gstyle(dot)) ///
						   xtitle("Firm age") ytitle("`ytitle'") ///
						   legend(order(`leglist') cols(4) region(lcolor(white))) ///
						   name("`var'_`moment'_cohort", replace)
						   graph_export_eps_pdf_png "9_lifecycle_firm_`var'_`moment'_cohort", replace
					local plotlist = ""
					local leglist = ""
					local i = 1 
					foreach yy in `list' {
						local col: word `i' of `listcolor'
						local mcol: word `i' of `listcolor'
						local pat: word `i' of `listpattern'
						local sym: word `i' of `listmsymbol'					
						local plotlist = "`plotlist' (${LineStyle} `var'_`moment' year if yob_firm == `yy', lcolor(`col') lpattern(`pat') mcolor(`col') msymbol(`sym'))"
						local leglist = `"`leglist' `i' "`yy'""'
						local i=`i'+1
					}
					tw `plotlist', ///
						   graphregion(color(white)) plotregion(lcolor(black) margin(l=0 r=0 b=0 t=0)) ///
						   xlabel(, gmin gmax grid gstyle(dot)) ylabel(, gmin gmax grid gstyle(dot)) ///
						   xtitle("YEar") ytitle("`ytitle'") ///
						   legend(order(`leglist') cols(4) region(lcolor(white))) ///
						   name("`var'_`moment'_time", replace)
						   graph_export_eps_pdf_png "9_lifecycle_firm_`var'_`moment'_time", replace
				}
			}
			
			
			* dec2021: removed
			/*			
			
			* load data
			local list = "${earn}_norm fe fye fye_init fye_post"
			local varlist = "fye"

			insheet using "${DIR_EXPORT}/csv/9_lifecycle_worker.csv", clear
			gen year = yob+${age}
// 			* normalize pay at initial age by education-gender-cohort groups
// 			foreach var in `list' {
// 				bys ${edu} ${gender} ( ${age} ): gen temp = `var'_mean[1]
// 				replace `var'_mean = `var'_mean-temp
// 				drop temp
// 			}			
			
			* plot

			foreach gender in "Male" "Female" {
				foreach moment in mean var {
					if "`moment'" == "mean" local ytitle = "Mean"
					else local ytitle = "Variance"
					forvalue educ = 1/2 {
						local list = "1950 1955 1960 1965 1970 1975 1980 1985"
						local plotlist = ""
						local leglist = ""
						local i = 1 
						foreach yy in `list' {
							local col: word `i' of `listcolor'
							local mcol: word `i' of `listcolor'
							local pat: word `i' of `listpattern'
							local sym: word `i' of `listmsymbol'						
							local plotlist = `"`plotlist' (${LineStyle} fye_`moment' ${age} if ${gender} == "`gender'" & ${edu} == `educ' & yob == `yy', lcolor(`col') lpattern(`pat') mcolor(`col') msymbol(`sym'))"'
							local leglist = `"`leglist' `i' "`yy'""'
							local i=`i'+1
						}
						tw `plotlist', ///
							graphregion(color(white)) plotregion(lcolor(black) margin(l=0 r=0 b=0 t=0)) ///
							xlabel(, gmin gmax grid gstyle(dot)) ylabel(, gmin gmax grid gstyle(dot)) ///
							xtitle("Worker age") ytitle("`ytitle'") ///
							legend(order(`leglist') cols(4) region(lcolor(none) fcolor(none)) size(${LegendSize}) margin(0 0 0 0)) ///
							name("`moment'_`gender'`educ'", replace)
							graph_export_eps_pdf_png "9_lifecycle_worker_`moment'_`gender'_educ`educ'_cohort", replace
						local plotlist = ""
						local leglist = ""
						local i = 1 
						foreach yy in `list' {
							local col: word `i' of `listcolor'
							local mcol: word `i' of `listcolor'
							local pat: word `i' of `listpattern'
							local sym: word `i' of `listmsymbol'						
							local plotlist = `"`plotlist' (${LineStyle} fye_`moment' year if ${gender} == "`gender'" & ${edu} == `educ' & yob == `yy', lcolor(`col') lpattern(`pat') mcolor(`col') msymbol(`sym'))"'
							local leglist = `"`leglist' `i' "`yy'""'
							local i=`i'+1
						}
						tw `plotlist', ///
							graphregion(color(white)) plotregion(lcolor(black) margin(l=0 r=0 b=0 t=0)) ///
							xlabel(, gmin gmax grid gstyle(dot)) ylabel(, gmin gmax grid gstyle(dot)) ///
							xtitle("Year") ytitle("`ytitle'") ///
							legend(order(`leglist') cols(4) region(lcolor(none) fcolor(none)) size(${LegendSize}) margin(0 0 0 0)) ///
							name("`moment'_`gender'`educ'", replace)
							graph_export_eps_pdf_png "9_lifecycle_worker_`moment'_`gender'_educ`educ'_time", replace
							
					}
					
				}
				
			}
			*/
			
			
			
		}
	
	}
	
}
			
