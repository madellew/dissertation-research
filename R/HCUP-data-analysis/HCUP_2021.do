/*Create indicator variables for diagnoses (DX) and procedure (PR) codes of interest*/

/*Created on: 11 Sep 2025 */
/*Modified on: 11 Nov 2025 */

/*Exclude if primary diagonis (i.e., ICD_DX1)*/

/*Create indicator variable for non-reccurent CDI in ICD DX2-40*/
foreach var of varlist I10_DX2-I10_DX40 {
generate CDI1`var' = `var' == "A0472"
}

/*Create indicator variable if CDI present in any ICD DX*/
egen CDI_np = rowtotal(CDI1I10_DX2-CDI1I10_DX40)

/*Get rid of unnecessary indicator variables*/
drop CDI1I10_DX2-CDI1I10_DX40

/*Create indicator variable for Candidal UTI in secondary ICD DX*/
foreach var of varlist I10_DX2-I10_DX40 {
generate CAN_UTI`var' = `var' == "B3741"
}
egen CAN_UTI_np = rowtotal(CAN_UTII10_DX2-CAN_UTII10_DX40)
drop CAN_UTII10_DX2-CAN_UTII10_DX40

/*Create indicator variable for Candidiasis (pulmonary, UTI, meningitis, 
endocarditis, sepsis) in secondary ICD DX*/
foreach var of varlist I10_DX2-I10_DX40 {
generate CAN`var' = `var' == "B371" | `var' == "B3741" | `var' == "B375" | `var' == "B376" | `var' == "B377"
}
egen CAN_np = rowtotal(CANI10_DX2-CANI10_DX40)
drop CANI10_DX2-CANI10_DX40
replace CAN_np = 1 if CAN_np >= 1

/*Create indicator variable for CLABSI in secondary ICD DX*/
foreach var of varlist I10_DX2-I10_DX40 {
generate CLABSI`var' = `var' == "T80211A"
}
egen CLABSI_np = rowtotal(CLABSII10_DX2-CLABSII10_DX40)
drop CLABSII10_DX2-CLABSII10_DX40

/*Create indicator variable for CLABSI in any ICD DX, fixed (remove subsequent encounter)*/
foreach var of varlist I10_DX1-I10_DX40 {
generate CLABSI`var' = `var' == "T80211A"
}
egen CLABSI = rowtotal(CLABSII10_DX1-CLABSII10_DX40)
drop CLABSII10_DX1-CLABSII10_DX40

/*Create indicator variable for CAUTI in secondary ICD DX*/
foreach var of varlist I10_DX2-I10_DX40 {
generate CAUTI`var' = `var' == "T83511A"
}
egen CAUTI_np = rowtotal(CAUTII10_DX2-CAUTII10_DX40)
drop CAUTII10_DX2-CAUTII10_DX40

/*Create indicator variable for CAUTI in any ICD DX, fixed (remove subsequent encounter*/
foreach var of varlist I10_DX1-I10_DX40 {
generate CAUTI`var' = `var' == "T83511A"
}
egen CAUTI = rowtotal(CAUTII10_DX1-CAUTII10_DX40)
drop CAUTII10_DX1-CAUTII10_DX40
replace CAUTI = 1 if CAUTI == 2

/*Create indicator variable for HAI in any ICD DX*/
foreach var of varlist I10_DX1-I10_DX40 {
generate HAI`var' = `var' == "Y95"
}
egen HAI = rowtotal(HAII10_DX1-HAII10_DX40)
drop HAII10_DX1-HAII10_DX40

/*Create indicator variable for antibiotic exposure in any ICD DX*/
foreach var of varlist I10_DX1-I10_DX40 {
generate ABX`var' = `var' == "Z792"
}
egen ABX = rowtotal(ABXI10_DX1-ABXI10_DX40)
drop ABXI10_DX1-ABXI10_DX40

/*Create indicator variable for MRSA carrier in any ICD DX*/
foreach var of varlist I10_DX1-I10_DX40 {
generate MRSA`var' = `var' == "Z22322"
}
egen MRSA = rowtotal(MRSAI10_DX1-MRSAI10_DX40)
drop MRSAI10_DX1-MRSAI10_DX40

/*Create indicator variable for urinary catheter placement in any ICD PR*/
foreach var of varlist I10_PR1-I10_PR25 {
generate URCATH`var' = `var' == "0T9B70Z"
}
egen URCATH = rowtotal(URCATHI10_PR1-URCATHI10_PR25)
drop URCATHI10_PR1-URCATHI10_PR25

/*Create indicator variable for CVC placement in any ICD PR*/
foreach var of varlist I10_PR1-I10_PR25 {
generate CVC`var' = `var' == "05HM33Z" | `var' == "05HN33Z" | `var' == "05H533Z" | `var' == "05H633Z" | `var' == "05H733Z" | `var' == "05H833Z" | `var' == "06HM33Z" | `var' == "06HN33Z"
}
egen CVC = rowtotal(CVCI10_PR1-CVCI10_PR25)
drop CVCI10_PR1-CVCI10_PR25
