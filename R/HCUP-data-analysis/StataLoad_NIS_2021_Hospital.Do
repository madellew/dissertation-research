/*****************************************************************************
* Stataload_NIS_2021_Hospital.Do
* This program will load the NIS_2021_Hospital ASCII File into Stata.
* Because Stata loads the entire file into memory, it may not be possible
* to load every data element for large files.  If necessary, edit this
* program to change the memory size or to load only selected data elements.
* Created on 12/22/2023.
*****************************************************************************/

*** Read data elements from the ASCII file ***
infix double DISCWT                       1-  11   ///
      byte   HOSP_BEDSIZE                12-  13   ///
      byte   HOSP_DIVISION               14-  15   ///
      byte   HOSP_LOCTEACH               16-  17   ///
      long   HOSP_NIS                    18-  22   ///
      byte   HOSP_REGION                 23-  24   ///
      byte   H_CONTRL                    25-  26   ///
      int    NIS_STRATUM                 27-  30   ///
      long   N_DISC_U                    31-  38   ///
      int    N_HOSP_U                    39-  42   ///
      long   S_DISC_U                    43-  50   ///
      int    S_HOSP_U                    51-  54   ///
      long   TOTAL_DISC                  55-  60   ///
      int    YEAR                        61-  64   ///
      using  "NIS_2021_Hospital.ASC"

***  Assign labels to the data elements ***
label var DISCWT                   "NIS discharge weight"
label var HOSP_BEDSIZE             "Relative bed size category of hospital (STRATA)"
label var HOSP_DIVISION            "Census Division of hospital (STRATA)"
label var HOSP_LOCTEACH            "Location/teaching status of hospital (STRATA)"
label var HOSP_NIS                 "NIS hospital number"
label var HOSP_REGION              "Region of hospital"
label var H_CONTRL                 "Control/ownership of hospital (STRATA)"
label var NIS_STRATUM              "NIS hospital stratum"
label var N_DISC_U                 "Number of universe discharges in the stratum"
label var N_HOSP_U                 "Number of universe hospitals in the stratum"
label var S_DISC_U                 "Number of sample discharges in the stratum"
label var S_HOSP_U                 "Number of sample hospitals in the stratum"
label var TOTAL_DISC               "Total number of discharges from this hospital in the NIS"
label var YEAR                     "Calendar year"

*** Convert special values to missing values ***
recode DISCWT                    (-99.9999999 -88.8888888 -66.6666666=.)
recode HOSP_BEDSIZE              (-9 -8 -6 -5=.)
recode HOSP_DIVISION             (-9 -8 -6 -5=.)
recode HOSP_LOCTEACH             (-9 -8 -6 -5=.)
recode HOSP_NIS                  (-9999 -8888 -6666=.)
recode HOSP_REGION               (-9 -8 -6 -5=.)
recode H_CONTRL                  (-9 -8 -6 -5=.)
recode NIS_STRATUM               (-999 -888 -666=.)
recode N_DISC_U                  (-9999999 -8888888 -6666666=.)
recode N_HOSP_U                  (-999 -888 -666=.)
recode S_DISC_U                  (-9999999 -8888888 -6666666=.)
recode S_HOSP_U                  (-999 -888 -666=.)
recode TOTAL_DISC                (-99999 -88888 -66666=.)
recode YEAR                      (-999 -888 -666=.)

save "NIS_2021_Hospital.dta", replace
