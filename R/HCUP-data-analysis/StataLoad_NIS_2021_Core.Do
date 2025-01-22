/*****************************************************************************
* Stataload_NIS_2021_Core.Do
* This program will load the NIS_2021_Core ASCII File into Stata.
* Because Stata loads the entire file into memory, it may not be possible
* to load every data element for large files.  If necessary, edit this
* program to change the memory size or to load only selected data elements.
* Created on 12/22/2023.
*****************************************************************************/

*** Read data elements from the ASCII file ***
infix int    AGE                          1-   3   ///
      byte   AGE_NEONATE                  4-   5   ///
      byte   AMONTH                       6-   7   ///
      byte   AWEEKEND                     8-   9   ///
      byte   DIED                        10-  11   ///
      double DISCWT                      12-  22   ///
      byte   DISPUNIFORM                 23-  24   ///
      byte   DQTR                        25-  26   ///
      int    DRG                         27-  29   ///
      byte   DRGVER                      30-  31   ///
      int    DRG_NoPOA                   32-  34   ///
      byte   ELECTIVE                    35-  36   ///
      byte   FEMALE                      37-  38   ///
      int    HCUP_ED                     39-  41   ///
      byte   HOSP_DIVISION               42-  43   ///
      long   HOSP_NIS                    44-  48   ///
      int    I10_BIRTH                   49-  51   ///
      int    I10_DELIVERY                52-  54   ///
      str    I10_DX1                     55-  61   ///
      str    I10_DX2                     62-  68   ///
      str    I10_DX3                     69-  75   ///
      str    I10_DX4                     76-  82   ///
      str    I10_DX5                     83-  89   ///
      str    I10_DX6                     90-  96   ///
      str    I10_DX7                     97- 103   ///
      str    I10_DX8                    104- 110   ///
      str    I10_DX9                    111- 117   ///
      str    I10_DX10                   118- 124   ///
      str    I10_DX11                   125- 131   ///
      str    I10_DX12                   132- 138   ///
      str    I10_DX13                   139- 145   ///
      str    I10_DX14                   146- 152   ///
      str    I10_DX15                   153- 159   ///
      str    I10_DX16                   160- 166   ///
      str    I10_DX17                   167- 173   ///
      str    I10_DX18                   174- 180   ///
      str    I10_DX19                   181- 187   ///
      str    I10_DX20                   188- 194   ///
      str    I10_DX21                   195- 201   ///
      str    I10_DX22                   202- 208   ///
      str    I10_DX23                   209- 215   ///
      str    I10_DX24                   216- 222   ///
      str    I10_DX25                   223- 229   ///
      str    I10_DX26                   230- 236   ///
      str    I10_DX27                   237- 243   ///
      str    I10_DX28                   244- 250   ///
      str    I10_DX29                   251- 257   ///
      str    I10_DX30                   258- 264   ///
      str    I10_DX31                   265- 271   ///
      str    I10_DX32                   272- 278   ///
      str    I10_DX33                   279- 285   ///
      str    I10_DX34                   286- 292   ///
      str    I10_DX35                   293- 299   ///
      str    I10_DX36                   300- 306   ///
      str    I10_DX37                   307- 313   ///
      str    I10_DX38                   314- 320   ///
      str    I10_DX39                   321- 327   ///
      str    I10_DX40                   328- 334   ///
      byte   I10_INJURY                 335- 336   ///
      byte   I10_MULTINJURY             337- 338   ///
      byte   I10_NDX                    339- 340   ///
      byte   I10_NPR                    341- 342   ///
      str    I10_PR1                    343- 349   ///
      str    I10_PR2                    350- 356   ///
      str    I10_PR3                    357- 363   ///
      str    I10_PR4                    364- 370   ///
      str    I10_PR5                    371- 377   ///
      str    I10_PR6                    378- 384   ///
      str    I10_PR7                    385- 391   ///
      str    I10_PR8                    392- 398   ///
      str    I10_PR9                    399- 405   ///
      str    I10_PR10                   406- 412   ///
      str    I10_PR11                   413- 419   ///
      str    I10_PR12                   420- 426   ///
      str    I10_PR13                   427- 433   ///
      str    I10_PR14                   434- 440   ///
      str    I10_PR15                   441- 447   ///
      str    I10_PR16                   448- 454   ///
      str    I10_PR17                   455- 461   ///
      str    I10_PR18                   462- 468   ///
      str    I10_PR19                   469- 475   ///
      str    I10_PR20                   476- 482   ///
      str    I10_PR21                   483- 489   ///
      str    I10_PR22                   490- 496   ///
      str    I10_PR23                   497- 503   ///
      str    I10_PR24                   504- 510   ///
      str    I10_PR25                   511- 517   ///
      int    I10_SERVICELINE            518- 520   ///
      double KEY_NIS                    521- 530   ///
      long   LOS                        531- 535   ///
      byte   MDC                        536- 537   ///
      byte   MDC_NoPOA                  538- 539   ///
      int    NIS_STRATUM                540- 543   ///
      byte   PAY1                       544- 545   ///
      byte   PCLASS_ORPROC              546- 547   ///
      int    PL_NCHS                    548- 550   ///
      int    PRDAY1                     551- 553   ///
      int    PRDAY2                     554- 556   ///
      int    PRDAY3                     557- 559   ///
      int    PRDAY4                     560- 562   ///
      int    PRDAY5                     563- 565   ///
      int    PRDAY6                     566- 568   ///
      int    PRDAY7                     569- 571   ///
      int    PRDAY8                     572- 574   ///
      int    PRDAY9                     575- 577   ///
      int    PRDAY10                    578- 580   ///
      int    PRDAY11                    581- 583   ///
      int    PRDAY12                    584- 586   ///
      int    PRDAY13                    587- 589   ///
      int    PRDAY14                    590- 592   ///
      int    PRDAY15                    593- 595   ///
      int    PRDAY16                    596- 598   ///
      int    PRDAY17                    599- 601   ///
      int    PRDAY18                    602- 604   ///
      int    PRDAY19                    605- 607   ///
      int    PRDAY20                    608- 610   ///
      int    PRDAY21                    611- 613   ///
      int    PRDAY22                    614- 616   ///
      int    PRDAY23                    617- 619   ///
      int    PRDAY24                    620- 622   ///
      int    PRDAY25                    623- 625   ///
      byte   RACE                       626- 627   ///
      double TOTCHG                     628- 637   ///
      byte   TRAN_IN                    638- 639   ///
      byte   TRAN_OUT                   640- 641   ///
      int    YEAR                       642- 645   ///
      byte   ZIPINC_QRTL                646- 647   ///
      using  "NIS_2021_Core.ASC"

***  Assign labels to the data elements ***
label var AGE                      "Age in years at admission"
label var AGE_NEONATE              "Neonatal age (first 28 days after birth) indicator"
label var AMONTH                   "Admission month"
label var AWEEKEND                 "Admission day is a weekend"
label var DIED                     "Died during hospitalization"
label var DISCWT                   "NIS discharge weight"
label var DISPUNIFORM              "Disposition of patient (uniform)"
label var DQTR                     "Discharge quarter"
label var DRG                      "DRG in effect on discharge date"
label var DRGVER                   "DRG grouper version used on discharge date"
label var DRG_NoPOA                "DRG in use on discharge date, calculated without POA"
label var ELECTIVE                 "Elective versus non-elective admission"
label var FEMALE                   "Indicator of sex"
label var HCUP_ED                  "HCUP Emergency Department service indicator"
label var HOSP_DIVISION            "Census Division of hospital (STRATA)"
label var HOSP_NIS                 "NIS hospital number"
label var I10_BIRTH                "ICD-10-CM Birth Indicator"
label var I10_DELIVERY             "ICD-10-CM Delivery Indicator"
label var I10_DX1                  "ICD-10-CM Diagnosis 1"
label var I10_DX2                  "ICD-10-CM Diagnosis 2"
label var I10_DX3                  "ICD-10-CM Diagnosis 3"
label var I10_DX4                  "ICD-10-CM Diagnosis 4"
label var I10_DX5                  "ICD-10-CM Diagnosis 5"
label var I10_DX6                  "ICD-10-CM Diagnosis 6"
label var I10_DX7                  "ICD-10-CM Diagnosis 7"
label var I10_DX8                  "ICD-10-CM Diagnosis 8"
label var I10_DX9                  "ICD-10-CM Diagnosis 9"
label var I10_DX10                 "ICD-10-CM Diagnosis 10"
label var I10_DX11                 "ICD-10-CM Diagnosis 11"
label var I10_DX12                 "ICD-10-CM Diagnosis 12"
label var I10_DX13                 "ICD-10-CM Diagnosis 13"
label var I10_DX14                 "ICD-10-CM Diagnosis 14"
label var I10_DX15                 "ICD-10-CM Diagnosis 15"
label var I10_DX16                 "ICD-10-CM Diagnosis 16"
label var I10_DX17                 "ICD-10-CM Diagnosis 17"
label var I10_DX18                 "ICD-10-CM Diagnosis 18"
label var I10_DX19                 "ICD-10-CM Diagnosis 19"
label var I10_DX20                 "ICD-10-CM Diagnosis 20"
label var I10_DX21                 "ICD-10-CM Diagnosis 21"
label var I10_DX22                 "ICD-10-CM Diagnosis 22"
label var I10_DX23                 "ICD-10-CM Diagnosis 23"
label var I10_DX24                 "ICD-10-CM Diagnosis 24"
label var I10_DX25                 "ICD-10-CM Diagnosis 25"
label var I10_DX26                 "ICD-10-CM Diagnosis 26"
label var I10_DX27                 "ICD-10-CM Diagnosis 27"
label var I10_DX28                 "ICD-10-CM Diagnosis 28"
label var I10_DX29                 "ICD-10-CM Diagnosis 29"
label var I10_DX30                 "ICD-10-CM Diagnosis 30"
label var I10_DX31                 "ICD-10-CM Diagnosis 31"
label var I10_DX32                 "ICD-10-CM Diagnosis 32"
label var I10_DX33                 "ICD-10-CM Diagnosis 33"
label var I10_DX34                 "ICD-10-CM Diagnosis 34"
label var I10_DX35                 "ICD-10-CM Diagnosis 35"
label var I10_DX36                 "ICD-10-CM Diagnosis 36"
label var I10_DX37                 "ICD-10-CM Diagnosis 37"
label var I10_DX38                 "ICD-10-CM Diagnosis 38"
label var I10_DX39                 "ICD-10-CM Diagnosis 39"
label var I10_DX40                 "ICD-10-CM Diagnosis 40"
label var I10_INJURY               "Injury ICD-10-CM diagnosis reported on record (1: First-listed injury; 2: Other than first-listed injury; 0: No injury)"
label var I10_MULTINJURY           "Multiple ICD-10-CM injuries reported on record"
label var I10_NDX                  "ICD-10-CM Number of diagnoses on this record"
label var I10_NPR                  "ICD-10-PCS Number of procedures on this record"
label var I10_PR1                  "ICD-10-PCS Procedure 1"
label var I10_PR2                  "ICD-10-PCS Procedure 2"
label var I10_PR3                  "ICD-10-PCS Procedure 3"
label var I10_PR4                  "ICD-10-PCS Procedure 4"
label var I10_PR5                  "ICD-10-PCS Procedure 5"
label var I10_PR6                  "ICD-10-PCS Procedure 6"
label var I10_PR7                  "ICD-10-PCS Procedure 7"
label var I10_PR8                  "ICD-10-PCS Procedure 8"
label var I10_PR9                  "ICD-10-PCS Procedure 9"
label var I10_PR10                 "ICD-10-PCS Procedure 10"
label var I10_PR11                 "ICD-10-PCS Procedure 11"
label var I10_PR12                 "ICD-10-PCS Procedure 12"
label var I10_PR13                 "ICD-10-PCS Procedure 13"
label var I10_PR14                 "ICD-10-PCS Procedure 14"
label var I10_PR15                 "ICD-10-PCS Procedure 15"
label var I10_PR16                 "ICD-10-PCS Procedure 16"
label var I10_PR17                 "ICD-10-PCS Procedure 17"
label var I10_PR18                 "ICD-10-PCS Procedure 18"
label var I10_PR19                 "ICD-10-PCS Procedure 19"
label var I10_PR20                 "ICD-10-PCS Procedure 20"
label var I10_PR21                 "ICD-10-PCS Procedure 21"
label var I10_PR22                 "ICD-10-PCS Procedure 22"
label var I10_PR23                 "ICD-10-PCS Procedure 23"
label var I10_PR24                 "ICD-10-PCS Procedure 24"
label var I10_PR25                 "ICD-10-PCS Procedure 25"
label var I10_SERVICELINE          "ICD-10-CM/PCS Hospital Service Line"
label var KEY_NIS                  "NIS record number"
label var LOS                      "Length of stay (cleaned)"
label var MDC                      "MDC in effect on discharge date"
label var MDC_NoPOA                "MDC in use on discharge date, calculated without POA"
label var NIS_STRATUM              "NIS hospital stratum"
label var PAY1                     "Primary expected payer (uniform)"
label var PCLASS_ORPROC            "Indicates operating room (major diagnostic or therapeutic) procedure on the record"
label var PL_NCHS                  "Patient Location: NCHS Urban-Rural Code"
label var PRDAY1                   "Number of days from admission to I10_PR1"
label var PRDAY2                   "Number of days from admission to I10_PR2"
label var PRDAY3                   "Number of days from admission to I10_PR3"
label var PRDAY4                   "Number of days from admission to I10_PR4"
label var PRDAY5                   "Number of days from admission to I10_PR5"
label var PRDAY6                   "Number of days from admission to I10_PR6"
label var PRDAY7                   "Number of days from admission to I10_PR7"
label var PRDAY8                   "Number of days from admission to I10_PR8"
label var PRDAY9                   "Number of days from admission to I10_PR9"
label var PRDAY10                  "Number of days from admission to I10_PR10"
label var PRDAY11                  "Number of days from admission to I10_PR11"
label var PRDAY12                  "Number of days from admission to I10_PR12"
label var PRDAY13                  "Number of days from admission to I10_PR13"
label var PRDAY14                  "Number of days from admission to I10_PR14"
label var PRDAY15                  "Number of days from admission to I10_PR15"
label var PRDAY16                  "Number of days from admission to I10_PR16"
label var PRDAY17                  "Number of days from admission to I10_PR17"
label var PRDAY18                  "Number of days from admission to I10_PR18"
label var PRDAY19                  "Number of days from admission to I10_PR19"
label var PRDAY20                  "Number of days from admission to I10_PR20"
label var PRDAY21                  "Number of days from admission to I10_PR21"
label var PRDAY22                  "Number of days from admission to I10_PR22"
label var PRDAY23                  "Number of days from admission to I10_PR23"
label var PRDAY24                  "Number of days from admission to I10_PR24"
label var PRDAY25                  "Number of days from admission to I10_PR25"
label var RACE                     "Race (uniform)"
label var TOTCHG                   "Total charges (cleaned)"
label var TRAN_IN                  "Transfer in indicator"
label var TRAN_OUT                 "Transfer out indicator"
label var YEAR                     "Calendar year"
label var ZIPINC_QRTL              "Median household income national quartile for patient ZIP Code"

*** Convert special values to missing values ***
recode AGE                       (-99 -88 -66=.)
recode AGE_NEONATE               (-9 -8 -6 -5=.)
recode AMONTH                    (-9 -8 -6 -5=.)
recode AWEEKEND                  (-9 -8 -6 -5=.)
recode DIED                      (-9 -8 -6 -5=.)
recode DISCWT                    (-99.9999999 -88.8888888 -66.6666666=.)
recode DISPUNIFORM               (-9 -8 -6 -5=.)
recode DQTR                      (-9 -8 -6 -5=.)
recode DRG                       (-99 -88 -66=.)
recode DRGVER                    (-9 -8 -6 -5=.)
recode DRG_NoPOA                 (-99 -88 -66=.)
recode ELECTIVE                  (-9 -8 -6 -5=.)
recode FEMALE                    (-9 -8 -6 -5=.)
recode HCUP_ED                   (-99 -88 -66=.)
recode HOSP_DIVISION             (-9 -8 -6 -5=.)
recode HOSP_NIS                  (-9999 -8888 -6666=.)
recode I10_BIRTH                 (-99 -88 -66=.)
recode I10_DELIVERY              (-99 -88 -66=.)
recode I10_INJURY                (-9 -8 -6 -5=.)
recode I10_MULTINJURY            (-9 -8 -6 -5=.)
recode I10_NDX                   (-9 -8 -6 -5=.)
recode I10_NPR                   (-9 -8 -6 -5=.)
recode I10_SERVICELINE           (-99 -88 -66=.)
recode KEY_NIS                   (-999999999 -888888888 -666666666=.)
recode LOS                       (-9999 -8888 -6666=.)
recode MDC                       (-9 -8 -6 -5=.)
recode MDC_NoPOA                 (-9 -8 -6 -5=.)
recode NIS_STRATUM               (-999 -888 -666=.)
recode PAY1                      (-9 -8 -6 -5=.)
recode PCLASS_ORPROC             (-9 -8 -6 -5=.)
recode PL_NCHS                   (-99 -88 -66=.)
recode PRDAY1                    (-99 -88 -66=.)
recode PRDAY2                    (-99 -88 -66=.)
recode PRDAY3                    (-99 -88 -66=.)
recode PRDAY4                    (-99 -88 -66=.)
recode PRDAY5                    (-99 -88 -66=.)
recode PRDAY6                    (-99 -88 -66=.)
recode PRDAY7                    (-99 -88 -66=.)
recode PRDAY8                    (-99 -88 -66=.)
recode PRDAY9                    (-99 -88 -66=.)
recode PRDAY10                   (-99 -88 -66=.)
recode PRDAY11                   (-99 -88 -66=.)
recode PRDAY12                   (-99 -88 -66=.)
recode PRDAY13                   (-99 -88 -66=.)
recode PRDAY14                   (-99 -88 -66=.)
recode PRDAY15                   (-99 -88 -66=.)
recode PRDAY16                   (-99 -88 -66=.)
recode PRDAY17                   (-99 -88 -66=.)
recode PRDAY18                   (-99 -88 -66=.)
recode PRDAY19                   (-99 -88 -66=.)
recode PRDAY20                   (-99 -88 -66=.)
recode PRDAY21                   (-99 -88 -66=.)
recode PRDAY22                   (-99 -88 -66=.)
recode PRDAY23                   (-99 -88 -66=.)
recode PRDAY24                   (-99 -88 -66=.)
recode PRDAY25                   (-99 -88 -66=.)
recode RACE                      (-9 -8 -6 -5=.)
recode TOTCHG                    (-999999999 -888888888 -666666666=.)
recode TRAN_IN                   (-9 -8 -6 -5=.)
recode TRAN_OUT                  (-9 -8 -6 -5=.)
recode YEAR                      (-999 -888 -666=.)
recode ZIPINC_QRTL               (-9 -8 -6 -5=.)

save "NIS_2021_Core.dta", replace
