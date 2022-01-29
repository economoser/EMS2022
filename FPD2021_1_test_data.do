********************************************************************************
* DESCRIPTION: Creates test dataset with linked employer-employee structure.
*
* CREDIT:      Please cite: Engbom, Niklas & Christian Moser & Jan Sauermann,
*              2021. "Firm Pay Dynamics," Working Paper.
*
* TIME STAMP:  October 27, 2021.
********************************************************************************


* clear memory
clear

* compute number of years
local n_years = ${year_max_data} - ${year_min_data} + 1

* generate year
gen int ${year} = .
forval y = $year_min_data/$year_max_data {
	count
	local n_obs_new = r(N) + ${n_obs_year}
	set obs `n_obs_new'
	replace ${year} = `y' if ${year} == .
}
label var ${year} "Year"

* generate observation number used to achieve a complete sort -- Note: This is needed so that seed produces a deterministic result.
gen long n_unique = _n
label var n_unique "Unique observation number"

* generate worker IDs
bys ${year} (n_unique): gen double ${id_pers} = _n
tostring ${id_pers}, replace
label var ${id_pers} "Worker ID (string)"

* generate firm IDs
// gen double ${id_firm} = runiform()^3 // Note: Taking the third power makes firm size distribution non-uniform.
gen double ${id_firm} = 1/exp(-ln(runiform())) // Note: Make firm size distribution Pareto with tail index 1
replace ${id_firm} = ceil(${id_firm}*${n_firms})

* have firms exit once in a while
forvalues yy = $year_min_data / $year_max_data {
	qui sum ${id_firm} if ${year} == `yy'
	local r = r(max)
	replace ${id_firm} = floor(`r'+ 1/exp(-ln(runiform()))) if year == `yy' & mod(${id_firm},89) == 1
}
tostring ${id_firm}, replace
label var ${id_firm} "Firm ID (string)"

* ensure that only a small fraction of workers switch firms across periods
bys ${id_pers} n_unique: replace ${id_firm} = ${id_firm}[1] if mod(real(${id_pers}), round(1/${test_switcher_share}))

* generate establishment IDs
// 	gen double ${id_est} = ${id_firm}
// 	label ${id_est} "Establishment ID"

* generate worker component of earnings
gen float pe = rnormal(${pe_mean}, ${pe_sd})
bys ${id_pers} (n_unique): replace pe = pe[1]

* generate firm component of earnings
gen float f_fe = rnormal(${f_fe_mean}, ${f_fe_sd})
bys ${id_firm} (n_unique): replace f_fe = f_fe[1]
label var f_fe "Firm fixed effect"

* generate firm-year component of earnings
gen float fy_fe = rnormal(${fy_fe_mean}, ${fy_fe_sd})
bys ${id_firm} ${year} (n_unique): replace fy_fe = fy_fe[1]
if "${gtools}" == "g" gegen float fy_fe_constant = mean(fy_fe), by(${year})
else bys ${year}: egen float fy_fe_constant = mean(fy_fe)
replace fy_fe = fy_fe - fy_fe_constant
drop fy_fe_constant
label var fy_fe "Firm-year fixed effect"

* generate time component of earnings
gen float ye = 0.1*sin((${year} - ${year_min_data})*_pi/2.5) + 0.01*`n_years'*(${year} - ${year_min_data})/(${year_max_data} - ${year_min_data}) // business cycle of period 5 years and amplitude 0.1 + average linear earnings growth of 1% p.a.

* generate residual component of earnings
gen float eps = rnormal(${eps_mean}, ${eps_sd})
if "${gtools}" == "g" gegen float eps_constant = mean(eps), by(${year})
else bys ${year}: egen float eps_constant = mean(eps)
replace eps = eps - eps_constant
drop eps_constant

* generate earnings (nominal LCU)
gen float ${earn} = exp(pe + f_fe + fy_fe + ye + eps)
label var ${earn} "Monthly earnings (nominal LCU)"

* generate weekly hours worked
gen int ${hours} = 20 + floor(runiform()*21)
label var hours "Weekly hours worked"

* generate wage (real LCU)
gen float ${wage} = ${earn}/(${hours}*365/12/7) // Note: 365/12/7 = mean number of weeks per month.
label var ${wage} "Hourly wage (nominal LCU)"
	
* summarize income measure and earnings components
sum ${earn} ${wage} pe f_fe fy_fe ye eps, d

* drop earnings components
drop pe f_fe fy_fe ye eps

* generate gender
gen byte ${gender} = .
replace ${gender} = 0 if mod(_n, 2)
replace ${gender} = 1 if !mod(_n, 2)
bys ${id_pers} (n_unique): replace ${gender} = ${gender}[1]
label var ${gender} "Gender"
label define gen_l 0 "male" 1 "female"
label val gender gen_l

* generate education
gen byte ${edu} = ceil(runiform()*4)
bys ${id_pers} (n_unique): replace ${edu} = ${edu}[1]
label var ${edu} "Education"

* generate age
gen int yob = ${year_min_data} - 23 + ceil(runiform()*5)
bys ${id_pers} (n_unique): replace yob = yob[1]
gen int ${age} = ${year} - yob
label var ${age} "Age (years)"
drop yob

* generate employment duration during current year
// 	gen int ${emp_dur} = ceil(runiform()*12)
// 	label var ${emp_dur} "Employment duration during current year (months)"

* generate occupational category
gen byte ${occ} = ceil(runiform()*25)
label var ${occ} "Occupational category"

* keep only relevant variables
keep ${id_pers} ${year} n_unique ${id_firm} ${earn} ${wage} ${gender} ${edu} ${age} ${occ}
order ${id_pers} ${year} n_unique ${id_firm} ${earn} ${wage} ${gender} ${edu} ${age} ${occ}
sort ${id_pers} ${year} n_unique
drop n_unique

* save
prog_comp_desc_sum_save "${DIR_TEMP}/test_data.dta"


*** firm-level covariates
* load worker-year-level data
use "${DIR_TEMP}/test_data.dta", clear

* collapse to firm-year level
gen long n = _n
label var n "Observation number"
${gtools}collapse (count) ${f_size_fek}=n, by(${id_firm} ${year}) fast

* generate observation number used to achieve a complete sort -- Note: This is needed so that seed produces a deterministic result.
gen long n_unique = _n
label var n_unique "Unique observation number"

* generate sector
gen float temp = runiform()
bys ${id_firm} (n_unique): replace temp = temp[1]
gen byte ${f_ind} = 1 + (temp > .2) + (temp > .4) + (temp > .6) + (temp > .8)
label var ${f_ind} "Sector"
drop temp

* make panel unbalanced by removing some random firms (especially in first few years) so as to generate variation in firm year of incorporation
gen float temp = runiform()*min(5, ${year_max_data} - ${year_min_data})
bys ${id_firm}: replace temp = temp[1]
count
bys ${id_firm}: drop if _n < temp
drop temp

* firm size
label var ${f_size_fek} "Firm size from FEK (number of employees)"

* firm financials in levels and per worker
foreach var in f_rev f_va f_ebit f_profit f_assets f_assets_lt f_equip f_struct f_cash f_debt f_equity f_inv {
	gen float ${`var'} = runiform()*${f_size_fek}
}
label var ${f_rev} "Sales (SEK)"
label var ${f_va} "Value added (SEK)"
label var ${f_ebit} "EBIT (SEK)"
label var ${f_profit} "Profits (SEK)"
label var ${f_assets} "Total assets (SEK)"
label var ${f_assets_lt} "Long-term assets (SEK)"
label var ${f_equip} "Equipment (SEK)"
label var ${f_struct} "Structures (SEK)"
label var ${f_cash} "Cash (SEK)"
label var ${f_debt} "Total debt (SEK)"
label var ${f_equity} "Total equity (SEK)"
label var ${f_inv} "Investment (SEK)"

* save
keep ///
	${id_firm} ${year} n_unique ///
	${f_ind} ${f_size_fek} ///
	${f_rev} ${f_va} ${f_ebit} ${f_profit} ${f_assets} ${f_assets_lt} ${f_equip} ${f_struct} ${f_cash} ${f_debt} ${f_equity} ${f_inv}
order ///
	${id_firm} ${year} n_unique ///
	${f_ind} ${f_size_fek} ///
	${f_rev} ${f_va} ${f_ebit} ${f_profit} ${f_assets} ${f_assets_lt} ${f_equip} ${f_struct} ${f_cash} ${f_debt} ${f_equity} ${f_inv}
sort ${id_firm} ${year} n_unique
drop n_unique
prog_comp_desc_sum_save "${DIR_TEMP}/test_data_firm.dta"
