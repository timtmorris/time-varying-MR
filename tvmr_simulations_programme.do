********************************************************************************
* Time-varying MR simulations
* Tim Morris
* 18/11/2021
********************************************************************************

cd "<PATHWAY>"
capture mkdir datasets
capture mkdir figures

********************************************************************************
* Creating bootstrap programme:
capture program drop mydiff
qui program mydiff, rclass
	  capture summarize _b_x1
	  local beta = r(mean)
	  capture summarize _b_x2
	  local beta = r(mean)
	  summarize _eq2_sim_1
	  local eq = r(mean)
	  return scalar diff = `eq'-`beta'
end

********************************************************************************
* Creating programmes: 
capture program drop ols_labreque_model
qui program define ols_labreque_model
args obs gx1 gx2 x1x2 x1y1 x1y2 y1x2 y1y2 x2y2 u yvar xvar
drop _all
set obs `obs'
gen e1=rnormal()
gen e2=rnormal()
gen e3=rnormal()
gen e4=rnormal()
gen temp1=0
replace temp1 = 1 if runiform() < 0.2
gen temp2=0
replace temp2 = 1 if runiform() < 0.2
gen snp=temp1+temp2
drop temp1 temp2
gen age1=0
gen age2=1
gen u=rnormal()
gen x1 = snp*`gx1' + age1 + u*`u' + e1
gen x2 = snp*`gx2' + age2 + x1*`x1x2' + u*`u' + e3
gen y2 = x1*`x1y2' + x2*`x2y2' + u*`u' + e4
local gx1 `gx1'
local gx2 `gx2'
local x1x2 `x1x2'
local x1y2 `x1y2'
local x2y2 `x2y2'
gen effy2x2=`x2y2'
local x1x2 `x1x2'
summ snp
local varg r(Var)
gen biasn=`gx1'*`gx2'*`varg'
summ x1
local varx1 r(Var)
gen biasd=`varx1'
gen bias=biasn/biasd
gen newx1x2=`x1x2'+bias
gen effy2x1=`x1y2'+newx1x2*`x2y2'
gen eff=eff`yvar'`xvar'
capture reg `yvar' `xvar'
capture reg `yvar' `xvar' x1
end

capture program drop iv_model
qui program define iv_model
args obs gx1 gx2 x1x2 x1y1 x1y2 y1x2 y1y2 x2y2 u yvar xvar
drop _all
set obs `obs'
gen e1=rnormal()
gen e2=rnormal()
gen e3=rnormal()
gen e4=rnormal()
gen temp1=0
replace temp1 = 1 if runiform() < 0.2
gen temp2=0
replace temp2 = 1 if runiform() < 0.2
gen snp=temp1+temp2
drop temp1 temp2
gen age1=0
gen age2=1
gen u=rnormal()
gen x1 = snp*`gx1' + age1 + u*`u' + e1
gen y1 = x1*`x1y1' + u*`u' + e2
gen x2 = snp*`gx2' + age2 + x1*`x1x2' + y1*`y1x2' + u*`u' + e3
gen y2 = x1*`x1y2' + x2*`x2y2' + y1*`y1y2' + u*`u' + e4
ivregress 2sls `yvar' (`xvar' = snp)
local gx1 `gx1'
local gx2 `gx2'
local x1x2 `x1x2'
local x1y1 `x1y1'
local x1y2 `x1y2'
local y1x2 `y1x2'
local x2y2 `x2y2'
local y1y2 `y1y2'
gen effy2x1=`x1y1'*`y1x2'*`x2y2'+`x1y1'*`y1y2'+`x1y2'+`x1x2'*`x2y2'+`x2y2'*`gx2'/`gx1'
local multipx2=1/(`gx1'*`x1x2'+`gx2'+`gx1'*`x1y1'*`y1x2')
gen effy1x1=`x1y1'
gen effy2x2=`gx1'*`multipx2'*(`x1y1'*`y1x2'*`x2y2'+`x1y1'*`y1y2'+`x1y2'+`x1x2'*`x2y2')+`gx2'*`multipx2'*`x2y2'
gen effy1x2=0
gen eff=eff`yvar'`xvar'
end


********************************************************************************
* Labreque & Swanson model:
foreach model in ols_labreque iv {
foreach target in "x1" "x2" {
foreach number in 0 1 {
	simulate _b _se eff, reps(1000): `model'_model 10000 0.5*`number' 0.5 0.3 0 0.4 0 0 0.4 0.3 y2 `target'
	summ *
	bootstrap r(diff), reps(1000): mydiff
	gen bias=_b[_bs_1]
	gen bias_se=_se[_bs_1]
	collapse (mean) _b_`target' _se_`target' bias bias_se _eq
	rename _eq2_sim_1 expected
	gen model="`model'"
	gen parameter="gx1"
	gen target="`target'"
	gen number=`number'
	save "datasets/`model'_`target'_gx1_`number'.dta", replace
	mkmat _b _se bias bias_se expected, matrix(`model'_`target'_gx1_`number')
	
	simulate _b _se eff, reps(1000): `model'_model 10000 0.5 0.5*`number' 0.3 0 0.4 0 0 0.4 0.3 y2 `target'
	summ *
	bootstrap r(diff), reps(1000): mydiff
	gen bias=_b[_bs_1]
	gen bias_se=_se[_bs_1]
	collapse (mean) _b_`target' _se_`target' bias bias_se _eq
	rename _eq2_sim_1 expected
	gen model="`model'"
	gen parameter="gx2"
	gen target="`target'"
	gen number=`number'
	save "datasets/`model'_`target'_gx2_`number'.dta", replace
	mkmat _b _se bias bias_se expected, matrix(`model'_`target'_gx2_`number')

	simulate _b _se eff, reps(1000): `model'_model 10000 0.5 0.5 0.3*`number' 0 0.4 0 0 0.4 0.3 y2 `target'
	summ *
	bootstrap r(diff), reps(1000): mydiff
	gen bias=_b[_bs_1]
	gen bias_se=_se[_bs_1]
	collapse (mean) _b_`target' _se_`target' bias bias_se _eq
	rename _eq2_sim_1 expected
	gen model="`model'"
	gen parameter="x1x2"
	gen target="`target'"
	gen number=`number'
	save "datasets/`model'_`target'_x1x2_`number'.dta", replace
	mkmat _b _se bias bias_se expected, matrix(`model'_`target'_x1x2_`number')

	simulate _b _se eff, reps(1000): `model'_model 10000 0.5 0.5 0.3 0 0.4*`number' 0 0 0.4 0.3 y2 `target'
	summ *
	bootstrap r(diff), reps(1000): mydiff
	gen bias=_b[_bs_1]
	gen bias_se=_se[_bs_1]
	collapse (mean) _b_`target' _se_`target' bias bias_se _eq
	rename _eq2_sim_1 expected
	gen model="`model'"
	gen parameter="x1y2"
	gen target="`target'"
	gen number=`number'
	save "datasets/`model'_`target'_x1y2_`number'.dta", replace
	mkmat _b _se bias bias_se expected, matrix(`model'_`target'_x1y2_`number')

	simulate _b _se eff, reps(1000): `model'_model 10000 0.5 0.5 0.3 0 0.4 0 0 0.4*`number' 0.3 y2 `target'
	summ *
	bootstrap r(diff), reps(1000): mydiff
	gen bias=_b[_bs_1]
	gen bias_se=_se[_bs_1]
	collapse (mean) _b_`target' _se_`target' bias bias_se _eq
	rename _eq2_sim_1 expected
	gen model="`model'"
	gen parameter="x2y2"
	gen target="`target'"
	gen number=`number'
	save "datasets/`model'_`target'_x2y2_`number'.dta", replace
	mkmat _b _se bias bias_se expected, matrix(`model'_`target'_x2y2_`number')

	simulate _b _se eff, reps(1000): `model'_model 10000 0.5 0.5 0.3 0 0.4 0 0 0.4 0.3*`number' y2 `target'
	summ *
	bootstrap r(diff), reps(1000): mydiff
	gen bias=_b[_bs_1]
	gen bias_se=_se[_bs_1]
	collapse (mean) _b_`target' _se_`target' bias bias_se _eq
	rename _eq2_sim_1 expected
	gen model="`model'"
	gen parameter="u"
	gen target="`target'"
	gen number=`number'
	save "datasets/`model'_`target'_u_`number'.dta", replace
	mkmat _b _se bias bias_se expected, matrix(`model'_`target'_u_`number')
}
}
}


* Creating tables:
clear
foreach model in ols_labreque iv {
foreach target in "x1" "x2" {
foreach parameter in gx1 gx2 x1x2 x1y2 x2y2 u {	
foreach number in 0 1 {
	append using "datasets/`model'_`target'_`parameter'_`number'.dta"
}
}
}
}
table parameter model number if target=="x1", c(mean _b_x1 mean _se_x1 mean bias) format(%9.3f)
table parameter model number if target=="x2", c(mean _b_x2 mean _se_x2 mean bias) format(%9.3f)


********************************************************************************
* Code to calcualte results for Figure 1 from Hardy et al (2010) and Timpson et 
* 		al (2009). Standard errors are computed using eq4 in Burgesset al (2017)
* 		ignoring covariance between SNP-exposure associations at different ages.

import delimited using https://github.com/timtmorris/time-varying-MR/blob/293801982fe792a496165b5dd5afe12dbad912d7/hardy_datapoints.txt, clear

* Calculating MR results:
gen timpson_sbp_beta=0.63
gen timpson_sbp_se=(0.93-0.33)/3.92
gen mr_sbp=timpson_sbp_beta/bmisds_mean
gen mr_sbp_se=sqrt(timpson_sbp_se^2/bmisds_mean^2 + timpson_sbp_beta^2*se^2 / bmisds_mean^4)

* Calculating genetically induced exposure differences:
local reference_beta=.1411684
local reference_se=.0311508
gen lchange=1/bmisds_mean
gen lchange_se=se/bmisds_mean^2
gen gied=bmisds_mean/`reference_beta'
gen gied_se=sqrt(se^2/`reference_beta'^2 + bmisds_mean^2*`reference_se'^2/`reference_beta'^4)

table age, c(mean bmisds_mean mean se) format(%9.3f) center
table age, c(mean mr_sbp mean mr_sbp_se mean lchange mean lchange_se) format(%9.2f) center
table age, c(mean gied mean gied_se) format(%9.2f) center


********************************************************************************
* Code to create Figure 4 (Hardy et al, Table 2): 

import delimited using https://github.com/timtmorris/time-varying-MR/blob/293801982fe792a496165b5dd5afe12dbad912d7/hardy_genotype.txt, clear

twoway (scatter tt age, color(navy%30) msymbol(T)) ///
	(scatter ta age, color(red%30) msymbol(X)) ///
	(scatter aa age, color(green%30) ///
	ytitle("BMI") xtitle("Age") legend(order (1 "TT" 2 "TA" 3 "AA") rows(1)) ///
	ylab(, nogrid) ///
	graphregion(color(white)))
graph export "figures/hardy.tif", replace width(2000)


********************************************************************************
* MR in the presence of reverse causation/feedback effects (Appendix 2):
capture program drop ols_feedback_model
qui program define ols_feedback_model
args obs gx1 gx2 x1x2 x1y1 x1y2 y1x2 y1y2 x2y2 u yvar xvar
drop _all
set obs `obs'
gen e1=rnormal()
gen e2=rnormal()
gen e3=rnormal()
gen e4=rnormal()
gen temp1=0
replace temp1 = 1 if runiform() < 0.2
gen temp2=0
replace temp2 = 1 if runiform() < 0.2
gen snp=temp1+temp2
drop temp1 temp2
gen age1=0
gen age2=1
gen u=rnormal()
gen x1 = snp*`gx1' + age1 + u*`u' + e1
gen y1 = x1*`x1y1' + u*`u' + e2
gen x2 = snp*`gx2' + age2 + x1*`x1x2' + y1*`y1x2' + u*`u' + e3
gen y2 = x1*`x1y2' + x2*`x2y2' + y1*`y1y2' + u*`u' + e4
capture reg `yvar' `xvar'
capture reg `yvar' `xvar' x1
local gx1 `gx1'
local gx2 `gx2'
local x1x2 `x1x2'
summ snp
local varg r(Var)
gen biasn=`gx1'*`gx2'*`varg'
summ x1
local varx1 r(Var)
gen biasd=`varx1'
gen bias=biasn/biasd
gen newx1x2=`x1x2'+bias
gen effy2x1=`x1y1'*`y1y2'+`x1y2'+(newx1x2*`x2y2')+`x1y1'*`y1x2'*`x2y2'
local x1y1 `x1y1'
local x1y2 `x1y2'
local y1x2 `y1x2'
local x2y2 `x2y2'
local y1y2 `y1y2'
gen effy2x2=`x2y2'
gen effy1x1=`x1y1'
gen effy1x2=0
gen eff=eff`yvar'`xvar'
end

foreach model in ols_feedback iv {
foreach target in "x1" "x2" {
foreach outcome in "y1" "y2" {
foreach number in 0 1 {
	simulate _b _se eff, reps(1000): `model'_model 10000 0.5*`number' 0.5 0.3 0.4 0.4 0.2 0.2 0.4 0.3 `outcome' `target'
	summ *
	bootstrap r(diff), reps(1000): mydiff
	gen bias=_b[_bs_1]
	gen bias_se=_se[_bs_1]
	collapse (mean) _b_`target' _se_`target' bias bias_se _eq
	rename _eq2_sim_1 expected
	gen model="`model'"
	gen parameter="gx1"
	gen target="`target'"
	gen number=`number'
	gen outcome="`outcome'"
	save "datasets/`model'_`outcome'_`target'_gx1_`number'.dta", replace
	mkmat _b _se bias bias_se expected, matrix(`model'_`outcome'_`target'_gx1_`number')

	simulate _b _se eff, reps(1000): `model'_model 10000 0.5 0.5*`number' 0.3 0.4 0.4 0.2 0.2 0.4 0.3 `outcome' `target'
	summ *
	bootstrap r(diff), reps(1000): mydiff
	gen bias=_b[_bs_1]
	gen bias_se=_se[_bs_1]
	collapse (mean) _b_`target' _se_`target' bias bias_se _eq
	rename _eq2_sim_1 expected
	gen model="`model'"
	gen parameter="gx2"
	gen target="`target'"
	gen number=`number'
	gen outcome="`outcome'"
	save "datasets/`model'_`outcome'_`target'_gx2_`number'.dta", replace
	mkmat _b _se bias bias_se, matrix(`model'_`outcome'_`target'_gx2_`number')

	simulate _b _se eff, reps(1000): `model'_model 10000 0.5 0.5 0.3*`number' 0.4 0.4 0.2 0.2 0.4 0.3 `outcome' `target'
	summ *
	bootstrap r(diff), reps(1000): mydiff
	gen bias=_b[_bs_1]
	gen bias_se=_se[_bs_1]
	collapse (mean) _b_`target' _se_`target' bias bias_se _eq
	rename _eq2_sim_1 expected
	gen model="`model'"
	gen parameter="x1x2"
	gen target="`target'"
	gen number=`number'
	gen outcome="`outcome'"
	save "datasets/`model'_`outcome'_`target'_x1x2_`number'.dta", replace
	mkmat _b _se bias bias_se, matrix(`model'_`outcome'_`target'_x1x2_`number')
	
	simulate _b _se eff, reps(1000): `model'_model 10000 0.5 0.5 0.3 0.4*`number' 0.4 0.2 0.2 0.4 0.3 `outcome' `target'
	summ *
	bootstrap r(diff), reps(1000): mydiff
	gen bias=_b[_bs_1]
	gen bias_se=_se[_bs_1]
	collapse (mean) _b_`target' _se_`target' bias bias_se _eq
	rename _eq2_sim_1 expected
	gen model="`model'"
	gen parameter="x1y1"
	gen target="`target'"
	gen number=`number'
	gen outcome="`outcome'"
	save "datasets/`model'_`outcome'_`target'_x1y1_`number'.dta", replace
	mkmat _b _se bias bias_se, matrix(`model'_`outcome'_`target'_x1y1_`number')

	simulate _b _se eff, reps(1000): `model'_model 10000 0.5 0.5 0.3 0.4 0.4*`number' 0.2 0.2 0.4 0.3 `outcome' `target'
	summ *
	bootstrap r(diff), reps(1000): mydiff
	gen bias=_b[_bs_1]
	gen bias_se=_se[_bs_1]
	collapse (mean) _b_`target' _se_`target' bias bias_se _eq
	rename _eq2_sim_1 expected
	gen model="`model'"
	gen parameter="x1y2"
	gen target="`target'"
	gen number=`number'
	gen outcome="`outcome'"
	save "datasets/`model'_`outcome'_`target'_x1y2_`number'.dta", replace
	mkmat _b _se bias bias_se, matrix(`model'_`outcome'_`target'_x1y2_`number')
	
	simulate _b _se eff, reps(1000): `model'_model 10000 0.5 0.5 0.3 0.4 0.4 0.2*`number' 0.2 0.4 0.3 `outcome' `target'
	summ *
	bootstrap r(diff), reps(1000): mydiff
	gen bias=_b[_bs_1]
	gen bias_se=_se[_bs_1]
	collapse (mean) _b_`target' _se_`target' bias bias_se _eq
	rename _eq2_sim_1 expected
	gen model="`model'"
	gen parameter="y1x2"
	gen target="`target'"
	gen number=`number'
	gen outcome="`outcome'"
	save "datasets/`model'_`outcome'_`target'_y1x2_`number'.dta", replace
	mkmat _b _se bias bias_se, matrix(`model'_`outcome'_`target'_y1x2_`number')
	
	simulate _b _se eff, reps(1000): `model'_model 10000 0.5 0.5 0.3 0.4 0.4 0.2 0.2*`number' 0.4 0.3 `outcome' `target'
	summ *
	bootstrap r(diff), reps(1000): mydiff
	gen bias=_b[_bs_1]
	gen bias_se=_se[_bs_1]
	collapse (mean) _b_`target' _se_`target' bias bias_se _eq
	rename _eq2_sim_1 expected
	gen model="`model'"
	gen parameter="y1y2"
	gen target="`target'"
	gen number=`number'
	gen outcome="`outcome'"
	save "datasets/`model'_`outcome'_`target'_y1y2_`number'.dta", replace
	mkmat _b _se bias bias_se, matrix(`model'_`outcome'_`target'_y1y2_`number')

	simulate _b _se eff, reps(1000): `model'_model 10000 0.5 0.5 0.3 0.4 0.4 0.2 0.2 0.4*`number' 0.3 `outcome' `target'
	summ *
	bootstrap r(diff), reps(1000): mydiff
	gen bias=_b[_bs_1]
	gen bias_se=_se[_bs_1]
	collapse (mean) _b_`target' _se_`target' bias bias_se _eq
	rename _eq2_sim_1 expected
	gen model="`model'"
	gen parameter="x2y2"
	gen target="`target'"
	gen number=`number'
	gen outcome="`outcome'"
	save "datasets/`model'_`outcome'_`target'_x2y2_`number'.dta", replace
	mkmat _b _se bias bias_se, matrix(`model'_`outcome'_`target'_x2y2_`number')

	simulate _b _se eff, reps(1000): `model'_model 10000 0.5 0.5 0.3 0.4 0.4 0.2 0.2 0.4 0.3*`number' `outcome' `target'
	summ *
	bootstrap r(diff), reps(1000): mydiff
	gen bias=_b[_bs_1]
	gen bias_se=_se[_bs_1]
	collapse (mean) _b_`target' _se_`target' bias bias_se _eq
	rename _eq2_sim_1 expected
	gen model="`model'"
	gen parameter="u"
	gen target="`target'"
	gen number=`number'
	gen outcome="`outcome'"
	save "datasets/`model'_`outcome'_`target'_u_`number'.dta", replace
	mkmat _b _se bias bias_se, matrix(`model'_`outcome'_`target'_u_`number')
}
}
}
}

* Creating tables:
clear
foreach model in ols_feedback iv {
	foreach target in "x1" "x2" {
		foreach outcome in "y1" "y2" {
			foreach parameter in gx1 gx2 x1x2 x1y1 x1y2 y1x2 y1y2 x2y2 u {	
				foreach number in 0 1 {
					append using "datasets/`model'_`outcome'_`target'_`parameter'_`number'.dta"
				}
			}
		}
	}
}
table parameter model number if target=="x1" & outcome=="y1", c(mean _b_x1 mean _se_x1 mean bias mean bias_se) format(%9.3f)
table parameter model number if target=="x1" & outcome=="y2", c(mean _b_x1 mean _se_x1 mean bias mean bias_se) format(%9.3f)
table parameter model number if target=="x2" & outcome=="y2", c(mean _b_x2 mean _se_x2 mean bias mean bias_se) format(%9.3f)
table parameter model number if target=="x2" & outcome=="y1", c(mean _b_x2 mean _se_x2 mean bias mean bias_se) format(%9.3f)


********************************************************************************
* Cross-sectional OLS confounding by genotype (Appendix 3):
capture program drop ols_naive
qui program define ols_naive
args obs gx1 gx2 x1x2 x1y2 x2y2 u yvar xvar cvar
drop _all
set obs `obs'
gen e1=rnormal()
gen e2=rnormal()
gen e3=rnormal()
gen e4=rnormal()
gen temp1=0
replace temp1 = 1 if runiform() < 0.2
gen temp2=0
replace temp2 = 1 if runiform() < 0.2
gen snp=temp1+temp2
drop temp1 temp2
gen age1=0
gen age2=1
gen u=rnormal()
gen x1 = snp*`gx1' + age1 + u*`u' + e1
gen x2 = snp*`gx2' + age2 + x1*`x1x2' + u*`u' + e3
gen y2 = x1*`x1y2' + x2*`x2y2' + u*`u' + e4
local gx1 `gx1'
local gx2 `gx2'
local x1x2 `x1x2'
local x1y2 `x1y2'
local x2y2 `x2y2'
gen effy2x2=`x2y2'
reg x2 x1
local newx1x2 _b[x1]
gen efftrue=`x1y2'+`newx1x2'*`x2y2'
gen effnaive=`x1y2'+`x1x2'*`x2y2'
capture reg `yvar' `xvar' `cvar'
end

foreach number in 0 1 {
	simulate _b _se effnaive efftrue, reps(1000): ols_naive 10000 0.5 0 0.3 0.4 0.4 0.3*`number' y2 x1 
	summ *
	collapse (mean) _b_x1 _se_x1 _eq2_sim_1 _eq2_sim_2
	gen bias_naive=_b_x-_eq2_sim_1
	gen bias_true=_b_x-_eq2_sim_2
	gen time="invariant"
	gen target="x1"
	gen pathway="gx1"
	gen number=`number'
	gen snp="No"
	save "datasets/ols_confounding_gx1_nosnp_`number'.dta", replace
	
	simulate _b _se effnaive efftrue, reps(1000): ols_naive 10000 0.5 0 0.3 0.4 0.4 0.3*`number' y2 x1 snp
	summ *
	collapse (mean) _b_x1 _se_x1 _eq2_sim_1 _eq2_sim_2
	gen bias_naive=_b_x-_eq2_sim_1
	gen bias_true=_b_x-_eq2_sim_2
	gen time="invariant"
	gen target="x1"
	gen pathway="gx1"
	gen number=`number'
	gen snp="Yes"
	save "datasets/ols_confounding_gx1_snp_`number'.dta", replace
	
	simulate _b _se effnaive efftrue, reps(1000): ols_naive 10000 0.5 0.5 0 0.4 0.4 0.3*`number' y2 x1 
	summ *
	collapse (mean) _b_x1 _se_x1 _eq2_sim_1 _eq2_sim_2
	gen bias_naive=_b_x-_eq2_sim_1
	gen bias_true=_b_x-_eq2_sim_2
	gen time="varying"
	gen pathway="x0x1"
	gen target="x1"
	gen number=`number'
	gen snp="No"
	save "datasets/ols_confounding_x0x1_nosnp_`number'.dta", replace
	
	simulate _b _se effnaive efftrue, reps(1000): ols_naive 10000 0.5 0.5 0 0.4 0.4 0.3*`number' y2 x1 snp
	summ *
	collapse (mean) _b_x1 _se_x1 _eq2_sim_1 _eq2_sim_2
	gen bias_naive=_b_x-_eq2_sim_1
	gen bias_true=_b_x-_eq2_sim_2
	gen time="varying"
	gen pathway="x0x1"
	gen target="x1"
	gen number=`number'
	gen snp="Yes"
	save "datasets/ols_confounding_x0x1_snp_`number'.dta", replace
	
	simulate _b _se effnaive efftrue, reps(1000): ols_naive 10000 0.5 0.5 0.3 0.4 0.4 0.3*`number' y2 x1 
	summ *
	collapse (mean) _b_x1 _se_x1 _eq2_sim_1 _eq2_sim_2
	gen bias_naive=_b_x-_eq2_sim_1
	gen bias_true=_b_x-_eq2_sim_2
	gen time="varying"
	gen target="x1"
	gen pathway="none"
	gen number=`number'
	gen snp="No"
	save "datasets/ols_confounding_full_nosnp_`number'.dta", replace
	
	simulate _b _se effnaive efftrue, reps(1000): ols_naive 10000 0.5 0.5 0.3 0.4 0.4 0.3*`number' y2 x1 snp
	summ *
	collapse (mean) _b_x1 _se_x1 _eq2_sim_1 _eq2_sim_2
	gen bias_naive=_b_x-_eq2_sim_1
	gen bias_true=_b_x-_eq2_sim_2
	gen time="varying"
	gen target="x1"
	gen pathway="none"
	gen number=`number'
	gen snp="Yes"
	save 	"datasets/ols_confounding_full_snp_`number'.dta", replace
}
clear
foreach number in 0 1 {
	append using "datasets/ols_confounding_gx1_nosnp_`number'.dta"
	append using "datasets/ols_confounding_gx1_snp_`number'.dta"
	append using "datasets/ols_confounding_x0x1_nosnp_`number'.dta"
	append using "datasets/ols_confounding_x0x1_snp_`number'.dta"
	append using "datasets/ols_confounding_full_nosnp_`number'.dta"
	append using "datasets/ols_confounding_full_snp_`number'.dta"
}
table number if snp=="No" & pathway=="gx1", c(mean _b_x1 mean _se_x1 mean _eq2_sim_1 mean _eq2_sim_2) format(%9.3f)
table number if snp=="Yes" & pathway=="gx1", c(mean _b_x1 mean _se_x1 mean _eq2_sim_1 mean _eq2_sim_2) format(%9.3f)
table number if snp=="No" & pathway=="x0x1", c(mean _b_x1 mean _se_x1 mean _eq2_sim_1 mean _eq2_sim_2) format(%9.3f)
table number if snp=="Yes" & pathway=="x0x1", c(mean _b_x1 mean _se_x1 mean _eq2_sim_1 mean _eq2_sim_2) format(%9.3f)
table number if snp=="No" & pathway=="none", c(mean _b_x1 mean _se_x1 mean _eq2_sim_1 mean _eq2_sim_2) format(%9.3f)
table number if snp=="Yes" & pathway=="none", c(mean _b_x1 mean _se_x1 mean _eq2_sim_1 mean _eq2_sim_2) format(%9.3f)


* END
