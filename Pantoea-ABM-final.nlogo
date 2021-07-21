;Application of Machine Learning Techniques to an Agent-Based Model of Pantoea
;
;Authors: Serena H. Chen, Pablo Londoño-Larrea, Andrew Stephen McGough, Amber N. Bible, Chathika Gunaratne,
;         Pablo A. Araujo-Granda, Jennifer L. Morrell-Falvey, Debsindhu Bhowmik, Miguel Fuentes-Cabrera
;..................................................................................................................................................

; Model variables
extensions [ csv matrix ] ; to export data and calculation
globals ; global variables used in model
[
  world                                                                      ; culture medium petri dish size
  steptime                                                                   ; step time [h]
  time_now                                                                   ; simulated time [hours]
  DST                                                                        ;standard deviation factor
  density_pa width_big_pa width_small_pa length_big_pa length_small_pa mmol_small_pa mmol_big_pa ini_biomass_pa  ; size parameters
  carbon hydrogen oxygen b-nitrogen molecular-weight e-eeq/mol dGpc gamma_cells dGp n dGr ;variables for thermodynamics (microorganims) pantoea
  column_matrix row_matrix rd ra rc energy_reaction synthesis_reaction global_reaction A fs fe;variables for matrix of helf-reactions pantoea
  uglucose uoxygen uamonium ubicarbonate                                     ; maximun uptake rates pantoea
  nglucose noxygen namonium nbicarbonate   Glu Amo Oxy Bic                   ; Culture medium nutrients availability
  gbacteria_ec gbacteria_pa gglucose gamonium goxygen gbicarbonate gdioxide gwater  ; global variables for graphical and numerical outputs
  glc_pa ox_pa                                                               ; global variables for cellular maintenance

  ; for mass balance
  ib-aceptors ib-donors iglucose iamonium ioxygen iwater idioxide ibicarbonate
  fb-aceptors fb-donors fglucose famonium foxygen fwater fdioxide fbicarbonate
  icb-aceptors inb-aceptors fcb-aceptors fnb-aceptors c-b-aceptors n-b-aceptors c-glucose c-co2 c-hco3 n-amonium n-obtained n-delivered
 icb-donors inb-donors fcb-donors fnb-donors c-b-donors n-b-donors
  c-obtained c-delivered
  %_e_c %_e_n yc/c

  name
  output_file
]

breed [b-aceptors b-aceptor] ; pantoea

b-aceptors-own
[
  pa-biomass                                             ; pantoea biomass [ mmol ]
  pa-biomass-reproduction                                ; reproduction treshold biomass [ mmol ]
  Rec_pa                                                    ; reproductions counter
  glucose-useful-pa oxygen-useful-pa amonium-useful-pa bicarbonate-useful-pa         ; usefuls nutrients variables
  time-viability_pa
  repr_time_pa
  inactive_pa
]

patches-own
  [
  glucose-medium         ; culture medium variable for the electron aceptor and C-Source in aerobic phase [ mmol ]
  oxygen-medium          ; culture medium variable for the electron aceptor in aerobic phase [ mmol ]
  amonium-medium         ; culture medium variable for the N-Source in aerobic phase [ mmol ]
  bicarbonate-medium     ; culture medium variable for the bicarbonate in aerobic phase [ mmol ]
  dioxide
  water
  num_microorg_pa
  ]
;....................................................................................................
; MODEL PROCEDURES
;
to setup
  clear-all                      ; clear all variables, graphical and numerical outputs
  reset-ticks                    ; reset tick counter
  set output_file (word "simul_"pmax"-"Microrganism"-"depth"-"umax_pa"-"diffusion-coefficient"-"
    Glucose"-"efficiency"-"Min/steptime"-"energy_maintenance_pa"-"
    rep_pa"-"max-time-viability_pa"-"Amonium"-"total-length-world"-"file-id"-"file-ir".dat")
  carefully
  [file-delete output_file] [ ]
  file-open output_file
  file-type "# pmax " file-print pmax
  file-type "# Microrganism " file-print Microrganism
  file-type "# depth " file-print depth
  file-type "# umax_pa " file-print umax_pa
  file-type "# diffusion-coefficient " file-print diffusion-coefficient
  file-type "# Glucose " file-print Glucose
  file-type "# efficiency " file-print efficiency
  file-type "# Min/steptime " file-print Min/steptime
  file-type "# energy_maintenance_pa " file-print energy_maintenance_pa
  file-type "# rep_pa " file-print rep_pa
  file-type "# max-time-viability_pa " file-print max-time-viability_pa
  file-type "# Amonium " file-print Amonium
  file-type "# total-length-world " file-print total-length-world
  file-type "# file-id " file-print file-id
  file-type "# file-ir " file-print file-ir
  file-print "#-------------------Beginning of the population data---------------"
  file-type "#time(hrs) "   file-print "population"
  file-type "0   "  file-print Microrganism
  file-close

  setup-thermodynamics-pantoea
  setup-model                    ; setup general model parameters to obtain the starting points accord to ORNL experimental data
  setup-time                     ; scaled time setting
  setup-medium                   ; formulation of the initial culture medium
  setup-pantoea
  setup-balance
end

to write-outputfile
  file-open output_file
  file-type time_now file-type "   "  carefully [ file-print (count b-aceptors) ] [ file-print "0" ]
  file-close
end

to go
  tick
  maintenance-requirements
  start_pantoea
  ask patches [color-nutrient count-turtles]
  diffuse glucose-medium diffusion-coefficient
  diffuse amonium-medium diffusion-coefficient
  diffuse oxygen-medium diffusion-coefficient
  diffuse bicarbonate-medium diffusion-coefficient
  setup-monitors
  do-plotting
  general_mass_balance

  ;if time_now >= 12 [images file-close stop]                             ; Stop simulation if number of individulas is higher than 500000 and close output-file
  if (sum [glucose-medium] of patches / world) <= 0 [file-close stop]     ; Stop simulation if the glucose global concentration is lower to zero and close output-file
  if (sum [amonium-medium] of patches / world) <= 0 [file-close stop]     ; Stop simulation if the amonium global concentration is lower to zero and close output-file
  if (sum [oxygen-medium] of patches / world) <= 0 [file-close stop]      ; Stop simulation if the oxygen global concentration is lower to zero and close output-file
  if (sum [bicarbonate-medium] of patches / world) <= 0 [file-close stop] ; Stop simulation if the global bicarbonate concentration is lower to zero and close output-file
  write-outputfile

end

to start_pantoea
ask b-aceptors
  [
    uptake_pa
    cellular-maintenance_pa
    metabolism_end_pa
    reproduce_pa
    if time-viability_pa >= max-time-viability_pa [set time-viability_pa max-time-viability_pa]
    ifelse time-viability_pa > Min/steptime [ set inactive_pa 0 ] [ set inactive_pa 1 ]
    if time-viability_pa <= 0 [set time-viability_pa 0]
    if time-viability_pa <= Min/steptime [ die ]  ; Command used to cause cell die
  ]
end

to resize
  set size (pa-biomass * molecular-weight / density_pa) ^ 0.33333 * 2
end

to setup-thermodynamics-pantoea
 ;--------------------------------------- procedure to setup e-acceptor and microorganims half-reactions matrix -------------------------------
 set column_matrix 12 ; the number of chemical species involved
 set row_matrix 1     ; this number is the number of reduction-half-reactions programmed
 set rd matrix:make-constant row_matrix column_matrix 0
 set ra matrix:make-constant row_matrix column_matrix 0
 set rc matrix:make-constant row_matrix column_matrix 0
 set energy_reaction matrix:make-constant row_matrix column_matrix 0
 set synthesis_reaction matrix:make-constant row_matrix column_matrix 0
 set global_reaction matrix:make-constant row_matrix column_matrix 0

; 0 1 ("CO2")     position for CO2 in matrix
; 0 2 ("H+")      position for H+ in matrix reactions
; 0 3 ("e-")      position for e- in matrix reactions
; 0 4 ("C6H12O6") positionfor glucose in matrix reactions
; 0 5 ("H2O")     position for H2O in matrix reactions
; 0 6 ("O2")      position for O2 in matrix reactions
; 0 7 ("HCO3-")   position for HCO3- in matrix reactions
; 0 8 ("NH4+")    position for NH4+ in matrix reactions
; 0 9 ("pantoea")  position for pantoea biomas in matrix reactions
; 0 10 ("∆G")     position for ∆G in matrix reactions

  ;  (rd) Electron donor --> Glucose :
  ;+ 0.25 CO2 + 1 H+ + 1 e- --> + 0.0417 C6H12O6 + 0.25 H2O [ ∆G = 41.35 KJ/e-eq ]

   matrix:set rd 0 1 (1 / 4)    ; CO2
   matrix:set rd 0 2 (1)        ; H+
   matrix:set rd 0 3 (1)        ; e-
   matrix:set rd 0 4 (-1 / 24)  ; C6H12O6
   matrix:set rd 0 5 (-1 / 4)   ; H2O
   matrix:set rd 0 10 (41.35)   ; ∆G

  ;  (ra) Electron acceptor --> O2 -> H2O :
  ;+ 0.25 O2 + 1 H+ + 1 e- --> + 0.5 H2O [ ∆G = -78.72 KJ/e-eq ]
   matrix:set ra 0 2 (1)        ; H+
   matrix:set ra 0 3 (1)        ; e-
   matrix:set ra 0 5 (-1 / 2)   ; H2O
   matrix:set ra 0 6 (1 / 4)    ; O2
   matrix:set ra 0 10 (-78.72)  ; ∆G

 ;  (rc) Biomass half reaction : C4.17H8O1.75N - bacteria general(pantoea), glucose , N-Source : NH4+

 ;0.1744 CO2 + 0.055 HCO3- + 0.055 NH4+ + 1 H+ + 1 e- --> 0.055 C'4.17'H'8'O'1.75'N'1 + 0.4175 H2O [ dG = 18.3112 KJ/e-eq ]

  set carbon 4.17 set hydrogen 8 set oxygen 1.75 set b-nitrogen 1                      ; pantoea empirical formula C4.17H8O1.75N
  set molecular-weight 12.011 * carbon + 1.008 * hydrogen + 16 * oxygen + 14 * b-nitrogen ; molecular weight of bacteria from empiral formula
  set e-eeq/mol (4 * carbon + hydrogen - 2 * oxygen - 3 * b-nitrogen)
  let b-co2 ((carbon - b-nitrogen) / e-eeq/mol)
  let b-hco3 (b-nitrogen / e-eeq/mol)
  let b-nh4 (b-nitrogen / e-eeq/mol)
  let b-bio (1 / e-eeq/mol)
  let b-h2o ((2 * carbon - oxygen + b-nitrogen) / e-eeq/mol)
  set dGpc (3.324 * molecular-weight / e-eeq/mol)

   matrix:set rc 0 1 (b-co2)      ; CO2
   matrix:set rc 0 2 (1)          ; H+
   matrix:set rc 0 3 (1)          ; e-
   matrix:set rc 0 5 (-1 * b-h2o) ; H2O
   matrix:set rc 0 7 (b-hco3)     ; HCO3-
   matrix:set rc 0 8 (b-nh4)      ; NH4+
   matrix:set rc 0 9 (-1 * b-bio) ; pantoea
   matrix:set rc 0 10 (dGpc)      ; ∆Gpc

; ---- Thermodynamic calculations to predict metabolic pathways using TEEM1 ---------------------------------------------------------------------

 set gamma_cells e-eeq/mol / carbon
 set dGp 35.09 - (matrix:get rd 0 10)
 ifelse dGp > 0 [set n 1][set n -1]
 set dGr (matrix:get ra 0 10) - (matrix:get rd 0 10)
 set A -1 * (((dGp / (efficiency ^ n))+(dGpc / efficiency))/(efficiency * dGr))
 set fs (1 / (1 + A))
 set fe (A / (1 + A))

 set energy_reaction ra matrix:- rd ; matrix with the stoichimetric coefficients to simulated celullar maintenance for pantoea
 set synthesis_reaction rc matrix:- rd
 set global_reaction ((fe matrix:* ra) matrix:+ (fs matrix:* rc) matrix:- (rd)) ; matrix that contain the stoichimetric coeficients for the pathway to simulated aerobic growth

;------ Calculations using stochiometric coefficients to stablish the maximun uptake rates for nutrients ----------------------------

set uglucose umax_pa * abs (matrix:get global_reaction 0 4 / matrix:get global_reaction 0 9) / carbon
set uoxygen umax_pa * abs (matrix:get global_reaction 0 6 / matrix:get global_reaction 0 9) / carbon
set uamonium umax_pa * abs (matrix:get global_reaction 0 8 / matrix:get global_reaction 0 9) / carbon
set ubicarbonate umax_pa * abs (matrix:get global_reaction 0 7 / matrix:get global_reaction 0 9) / carbon

end

to setup-model
 set world (total-length-world * total-length-world * depth) * 1e-15         ; Total Volume of World (L)
 set DST 0.10                                                   ; Standard Deviation for general model varaibles
 set nglucose 1 set noxygen 1 set namonium 1 set nbicarbonate 1

  set density_pa 1.1e-9                                                                         ; density in mg/um^3
  set width_big_pa 0.5                                                                          ; size of bacteria in um (from: Thermodynamic and kinetic controls on co-transport of Pantoea  agglomerans  cells and Zn through clean and iron oxide coated sand columns
  set width_small_pa width_big_pa / 2                                                       ; size of bacteria in um (from: Thermodynamic and kinetic controls on co-transport of Pantoea  agglomerans  cells and Zn through clean and iron oxide coated sand columns
  set length_big_pa 1.1                                                                     ; size of bacteria in um (from: Thermodynamic and kinetic controls on co-transport of Pantoea  agglomerans  cells and Zn through clean and iron oxide coated sand columns
  set length_small_pa length_big_pa  / 2                                                       ; size of bacteria in um (from: Thermodynamic and kinetic controls on co-transport of Pantoea  agglomerans  cells and Zn through clean and iron oxide coated sand columns
 set mmol_small_pa (pi * ((width_small_pa / 2) ^ 2) * length_small_pa) * density_pa / molecular-weight ; mmol
 set mmol_big_pa (pi * ((width_big_pa / 2) ^ 2) * length_big_pa) * density_pa / molecular-weight  ; ; mmol
 set ini_biomass_pa 0.75 * mmol_big_pa                                                       ; [mmol] initial biomass for each bacterium

end

to setup-time
  set steptime (Min/steptime / 60)       ; User define the number of minute for each step time using a slider labeled "Min_per_steptime"
end

to setup-medium
 set Glu random-normal (Glucose) (DST * Glucose)           ; Glucose initial concentration (user define using a slider labeled Succinate) [mM]
 set Amo random-normal (Amonium) (DST * Amonium)               ; Amonium initial concentration (user define using a slider labeled Amonium) [mM]
 set Oxy random-normal  (100) (0.1 * 100) ;                   ; Oxygen initial concentration (simulated saturated conditions) [mM]
 set Bic random-normal (1000) (DST * 1000)               ; Bicarbonate initial concentration (user defined) [mM]
 ask patches
   [
     set glucose-medium ((abs (random-normal (Glu) (DST * Glu)) * world) / (world-width * world-height)) ; Assign the glucose concentration to each cell according to the model's standar deviation
     set oxygen-medium ((abs (random-normal (Oxy) (DST * Oxy)) * world) / (world-width * world-height))              ; Assign the oxygen concentration to each cell according to the model's standar deviation
     set amonium-medium  ((abs (random-normal (Amo) (DST * Amo)) * world) / (world-width * world-height))             ; Assign the amonium concentration to each cell according to the model's standar deviation
     set bicarbonate-medium ((abs (random-normal (Bic) (DST * Bic)) * world) / (world-width * world-height))              ; Assign the bicarbonate concentration to each cell according to the model's standar deviation
     set dioxide 0                                                ; dioxide initial concentration  [mM]
     set water 0                                                       ;water initial concentration  [mM]
   ]

end

to setup-pantoea
  create-b-aceptors Microrganism                                                 ; User define the amount of initial active bacteria using a slider labaled Microrganism
   [
   set shape "bacteria"
   set color green         set heading (random 360)                              ; This color is for visual purpose only
   setxy random-xcor random-ycor                                                 ; Random inicial position on culture medium
   set pa-biomass abs random-normal (ini_biomass_pa * 0.75) (DST * ini_biomass_pa * 0.75); Initial biomass for each microorganism [mmol] according to model's standar deviation
   resize
   set time-viability_pa abs random-normal (max-time-viability_pa) (DST * max-time-viability_pa)
   set pa-biomass-reproduction abs random-normal (ini_biomass_pa) (DST * ini_biomass_pa)
   set repr_time_pa 0                                                           ; Set number if assuming seed in exponential growth behavior or <0> if is lag behavior
  ]
end

to setup-monitors
  set time_now ticks * steptime                                                                           ; Simulated time [ h ]
  set gbacteria_pa ((sum[pa-biomass] of b-aceptors / world) / 1000) * molecular-weight                    ; Biomass concentration [ mg/ml ]
  set gglucose (sum [glucose-medium]of patches / world)                                                   ; Glucose concentration [ mM ]
  set gamonium (sum [amonium-medium] of patches / world)                                                  ; Amonium concentration [ mM ]
  set goxygen (sum [oxygen-medium] of patches / world)                                                    ; Oxygen concentration [ mM ]
  set gbicarbonate (sum [bicarbonate-medium] of patches / world)                                          ; Bicarbonate concentration [ mM ]
  set gdioxide (sum [dioxide] of patches / world)                                                         ; Dioxide concentration [ mM ]
  set gwater (sum [water] of patches / world)                                                             ; water concentration [ mM ]

  reactor_balance
end

to maintenance-requirements
  ; according to Energy reaction: 0.0417 C6H12O6 + 0.25 O2 --> + 0.25 CO2 + 0.25 H2O
  set glc_pa (energy_maintenance_pa * steptime  * carbon / 6) ; [mol glucose / molBiomass]
  set ox_pa (glc_pa * abs (matrix:get energy_reaction 0 6) / abs (matrix:get energy_reaction 0 4)) ; [mol oxygen / molBiomass ]
end

to uptake_pa
  if inactive_pa = 0
  [
  let uptglc (abs (random-normal (uglucose) (DST * uglucose)) * pa-biomass * carbon * steptime)         ; [mmol] glucose by biomass uptake accord to INDISIM
  let uptox (abs (random-normal (uoxygen) (DST * uoxygen)) * pa-biomass * carbon * steptime)            ; [mmol] oxygen by biomass uptake accord to INDISIM
  let uptnh4 (abs (random-normal (uamonium) (DST * uamonium)) * pa-biomass * carbon * steptime)         ; [mmol] amonium by biomass uptake accord to INDISIM
  let upthco3 (abs (random-normal (ubicarbonate) (DST * ubicarbonate)) * pa-biomass * carbon * steptime)      ; [mmol] bicarbonate by biomass uptake accord to INDISIM

  let glc-available (abs (random-normal (nglucose) (DST * nglucose)) * glucose-medium * steptime)              ; INDISIM lineal model to determine glucose availability in culture medium [mmol]
  let nh4-available (abs (random-normal (namonium) (DST * namonium)) * amonium-medium * steptime)              ; INDISIM lineal model to determine amonium availability in culture medium [mmol]
  let ox-available (abs (random-normal (noxygen) (DST * noxygen)) * oxygen-medium * steptime)                  ; INDISIM lineal model to determine oxygen availability in culture medium [mmol]
  let hco3-available (abs (random-normal (nbicarbonate) (DST * nbicarbonate)) * bicarbonate-medium * steptime) ; INDISIM lineal model to determine bicarbonate availability in culture medium [mmol]

  ifelse glc-available <= uptglc [set glucose-useful-pa glc-available] [set glucose-useful-pa uptglc]             ; For glucose comparation between uptake by biomass and availability and takes the lowest value.
  ifelse nh4-available <= uptnh4 [set amonium-useful-pa nh4-available] [set amonium-useful-pa uptnh4]             ; For amonium comparation between uptake by biomass and availability and takes the lowest value.
  ifelse ox-available <= uptox [set oxygen-useful-pa ox-available] [set oxygen-useful-pa uptox]                   ; For oxygen comparation between uptake by biomass and availability and takes the lowest value.
  ifelse hco3-available <= upthco3 [set bicarbonate-useful-pa hco3-available] [set bicarbonate-useful-pa upthco3] ; For bicarbonate comparation between uptake by biomass and availability and takes the lowest value.

  set glucose-medium (glucose-medium - glucose-useful-pa)                         ; Update local quantity of glucose
  set amonium-medium (amonium-medium - amonium-useful-pa)                         ; Update local quantity of amonium
  set oxygen-medium (oxygen-medium - oxygen-useful-pa)                            ; Update local quantity of oxygen
  set bicarbonate-medium (bicarbonate-medium - bicarbonate-useful-pa)             ; Update local quantity of bicarbonate
  ]
end

to cellular-maintenance_pa
 if inactive_pa = 0
  [
  let glc-mant (glc_pa * pa-biomass)                                        ; Glucose celullar maintenance requirements by biomass [mmol]
  let ox-mant (ox_pa * pa-biomass)                                          ; Oxygen cellular maintenance requirements by biomass [mmol]

  let y min (list glc-mant glucose-useful-pa)                            ; takes the minimum value between maintenance requirements and uptake quantity for Succinate [mmol]
  let z min (list ox-mant oxygen-useful-pa)                             ; takes the minimum value between maintenance requirements and uptake quantity for Oxygen [mmol]

  let y1 (y / abs (matrix:get energy_reaction 0 4))                         ; glucose quantity by stoichometric coefficient
  let z1 (z / abs (matrix:get energy_reaction 0 6))                         ; oxygen quantity by stoichiometric coefficient

  let w1 min (list y1 z1)                                                ; limiting nutrient for aerobic maintenance

  let ndioxide (abs (matrix:get energy_reaction 0 1) * w1)            ; Carbon dioxide generation by cellular maintenance [mmol]
  let nwater (abs (matrix:get energy_reaction 0 5) * w1)              ; water generation by cellular maintenance [mmol]
  set dioxide (dioxide + ndioxide)                                    ; Carbon dioxide local quantity update
  set water (water + nwater)                           ; water local quantity update

  set glucose-useful-pa (glucose-useful-pa - abs (matrix:get energy_reaction 0 4) * w1)     ; Update of glucose uptaked  [mmol]
  set oxygen-useful-pa (oxygen-useful-pa - abs (matrix:get energy_reaction 0 6) * w1)     ; Update of Oxygen uptaked [mmol]

  ifelse (glucose-useful-pa > 0) and (oxygen-useful-pa > 0) and (amonium-useful-pa > 0) and (bicarbonate-useful-pa > 0)
      [aerobic_pa]  ; if uptaken updated quantities are higher than zero the microorganism could execute  aerobic parthway
  [
      set time-viability_pa time-viability_pa - Min/steptime
    ]
  ]
end


To aerobic_pa

  let a1 (glucose-useful-pa / abs (matrix:get global_reaction 0 4))      ; uptaken glucose updated divided by its stoichiometric coefficient
  let b (amonium-useful-pa / abs (matrix:get global_reaction 0 8))       ; uptaken amonium updated divided by its stoichiometric coefficient
  let c (oxygen-useful-pa / abs (matrix:get global_reaction 0 6))        ; uptaken oxygen updated divided by its stoichiometric coefficient
  let d (bicarbonate-useful-pa / abs (matrix:get global_reaction 0 7))   ; uptaken bicarbonate updated divided by its stoichiometric coefficient

  let x1 min (list a1 b c d)                                             ; limiting nutrient for pathway 1 in aerobic phase

  let nbiomass abs (matrix:get global_reaction 0 9) * x1                 ; Biomass generation by aerobic pathway 1 [mmol]
  let ndioxide abs (matrix:get global_reaction 0 1) * x1                 ; Carbon dioxide generation by aerobic pathway 1 [mmol]
  let nwater abs (matrix:get global_reaction 0 5) * x1                   ; water generation by aerobic pathway 1 [mmol]
  set pa-biomass (pa-biomass + nbiomass)                                 ; Individual biomass update
  set dioxide (dioxide + ndioxide)                                       ; Carbon dioxide local quantity update
  set water (water + nwater)                                             ; water local quantity update

  set glucose-useful-pa (glucose-useful-pa - abs (matrix:get global_reaction 0 4) * x1)       ; update electron aceptor bacterium quantity
  set amonium-useful-pa (amonium-useful-pa - abs (matrix:get global_reaction 0 8) * x1)       ; update N-source bacterium quantity
  set oxygen-useful-pa (oxygen-useful-pa - abs (matrix:get global_reaction 0 6) * x1)         ; update electron aceptor bacterium quantity
  set bicarbonate-useful-pa (bicarbonate-useful-pa - abs (matrix:get global_reaction 0 7) * x1)         ; update bicarbonate bacterium quantity

end

to metabolism_end_pa
  ifelse glucose-useful-pa > 0 [set glucose-medium (glucose-medium + glucose-useful-pa) set glucose-useful-pa 0][set glucose-useful-pa 0]   ; Release to culture medium the quantity not used
  ifelse amonium-useful-pa > 0 [set amonium-medium (amonium-medium + amonium-useful-pa) set amonium-useful-pa 0][set amonium-useful-pa 0]       ; Release to culture medium the quantity not used
  ifelse oxygen-useful-pa > 0 [set oxygen-medium (oxygen-medium + oxygen-useful-pa) set oxygen-useful-pa 0][set oxygen-useful-pa 0]         ; Release to culture medium the quantity not used
  ifelse bicarbonate-useful-pa > 0 [set bicarbonate-medium (bicarbonate-medium + bicarbonate-useful-pa) set bicarbonate-useful-pa 0][set bicarbonate-useful-pa 0]   ; Release to culture medium the quantity not used
end

to reproduce_pa

  if inactive_pa = 0
  [
   if pa-biomass >= pa-biomass-reproduction                         ; cell division happens when a treshold value is reached
    [
    set repr_time_pa repr_time_pa + Min/steptime
    if repr_time_pa > abs (random-normal rep_pa (0.1 * rep_pa))
      [
   ifelse count turtles-here <= pmax - 1                            ; "pmax" is the max number of turtles allowed in a patch
        [
        let pa-biomass-hatch pa-biomass * abs (random-normal 0.5 (0.1 * 0.5) ) ;; the mass from the original microorganism is not divided equally between the original and de newborn, it depends on a stochastic factor
        set pa-biomass pa-biomass - pa-biomass-hatch
        set repr_time_pa 0
          set pa-biomass-reproduction abs random-normal (ini_biomass_pa) (DST * ini_biomass_pa)
          hatch 1 [  set color green set heading (random 360 ) set repr_time_pa 0  set pa-biomass pa-biomass-hatch ]
        ]

        [ if count turtles-here >= pmax

          [let new_position one-of neighbors with [count turtles-here <  pmax ]
          if new_position != nobody

          [let pa-biomass-hatch pa-biomass * abs (random-normal 0.5 (0.1 * 0.5) ) ;; the mass from the original microorganism is not divided equally between the original and de newborn, it depends on a stochastic factor
          set pa-biomass pa-biomass - pa-biomass-hatch
                    set pa-biomass-reproduction abs random-normal (ini_biomass_pa) (DST * ini_biomass_pa)
          set repr_time_pa 0
          hatch 1 [  set color green set heading (random 360 ) set repr_time_pa 0  set pa-biomass pa-biomass-hatch move-to new_position ]
          ]
        ]]

    ]]]

end

to color-nutrient
  set pcolor black
end

to count-turtles
  ifelse any? b-aceptors-here
  [ set num_microorg_pa ( count b-aceptors-here) ]
  [ set num_microorg_pa 0 ]
end

; ____________________________ General mass balance proceedings _____________________________________________________________________________

to setup-balance
  set ib-aceptors (sum [pa-biomass] of b-aceptors)                                         ; Initial Biomass [mmol]
  set iglucose (sum[glucose-medium]of patches)                                             ; Initial Glucose [mmol]
  set iamonium (sum [amonium-medium] of patches)                                           ; Initial Amonium [mmol]
  set ioxygen (sum [oxygen-medium] of patches)                                             ; Initial Oxygen [mmol]
  set iwater (sum [water] of patches)                                                      ; Initial water [mmol]
  set idioxide (sum [dioxide] of patches)                                                  ; Initial Carbon dioxide [mmol]
  set ibicarbonate (sum [bicarbonate-medium] of patches)                                   ; Initial Bicarbonate [mmol]
end

to reactor_balance

  set fb-aceptors (gbacteria_pa * 1000 * world / molecular-weight)                         ; Actual Biomass [mmol]
  set fglucose (gglucose * world)                                                          ; Actual Glucose [mmol]
  set famonium (gamonium * world)                                                          ; Actual Amonium [mmol]
  set foxygen (goxygen * world)                                                            ; Actual Oxygen [mmol]
  set fwater (gwater * world)                                                              ; Actual water [mmol]
  set fdioxide (gdioxide * world)                                                          ; Actual Carbon dioxide [mmol]
  set fbicarbonate (gbicarbonate * world)                                                  ; Actual Bicarbonate [mmol]

end

to general_mass_balance
  set icb-aceptors (ib-aceptors * carbon * 12)                                             ; Initial C-Mic [mg]
  set inb-aceptors (ib-aceptors * b-nitrogen * 14)                                         ; Initial N-Mic [mg]
  set fcb-aceptors (fb-aceptors * carbon * 12)                                             ; Final C-Mic [mg]
  set fnb-aceptors (fb-aceptors * b-nitrogen * 14)                                         ; Final N-Mic [mg]
  set c-b-aceptors (fcb-aceptors - icb-aceptors)                                           ; Microorganism Carbon balance
  set n-b-aceptors (fnb-aceptors - inb-aceptors)                                           ; Microorganims Nitrogen balance


  set c-glucose ((fglucose - iglucose)  * 6 * 12 )                                         ; Glucose-carbon balance
  set n-amonium ((famonium - iamonium)  * 1 * 14 )                                         ; Amonium-nitrogen balance
  set c-co2 ((fdioxide - idioxide)* 1 * 12 )                                               ; Carbon dixoide-carbon balance
  set c-hco3 ((fbicarbonate - ibicarbonate) * 1 * 12)                                      ; Bicarbonate-carbon balance

  set n-obtained abs (n-b-aceptors )                                                       ; Total nitrogen obtained
  set n-delivered abs (n-amonium )                                                         ; Total nitrogen delivered
  set c-obtained abs (c-co2 + c-b-aceptors )                                               ; Total carbon obtained
  set c-delivered abs (c-glucose + c-hco3 )                                                ; Total carbon delivered

  set %_e_c abs ((c-delivered - c-obtained) / c-delivered ) * 100                          ; Carbon general mass balance porcentual error
 set %_e_n abs ((n-delivered - n-obtained) / n-delivered ) * 100                           ; Nitrogen general mass balance porcentual error
  set yc/c abs ((c-b-aceptors ) / c-glucose)                                               ; Population Yield growth

end

to do-plotting
  setup-monitors
  set-current-plot "Biomass"
  set-current-plot-pen "pantoea"
  carefully [ plotxy time_now count b-aceptors ] [  plotxy time_now 0 ]
end
@#$#@#$#@
GRAPHICS-WINDOW
586
10
960
385
-1
-1
1.46
1
10
1
1
1
0
0
0
1
-125
125
-125
125
1
1
1
Step time
30.0

BUTTON
431
211
497
244
NIL
Setup
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

SLIDER
9
14
181
47
efficiency
efficiency
0.01
1
0.37
0.01
1
NIL
HORIZONTAL

MONITOR
12
440
245
485
uGlucose [mol Glu/ molCcell h]
uglucose
4
1
11

MONITOR
248
440
450
485
uOxygen [mol O2/molCcell h]
uoxygen
4
1
11

MONITOR
451
440
660
485
uNH4+ [molNH4+/molCcells h]
uamonium
4
1
11

MONITOR
662
440
885
485
uHCO3- [molHCO3-/molCcells h]
ubicarbonate
4
1
11

SLIDER
8
46
180
79
umax_pa
umax_pa
0.1
10
2.0
0.1
1
NIL
HORIZONTAL

SLIDER
8
78
180
111
Min/steptime
Min/steptime
0.1
2
0.5
.5
1
[min]
HORIZONTAL

SLIDER
380
94
570
127
Glucose
Glucose
0
100
2.8
.1
1
[mM]
HORIZONTAL

SLIDER
8
271
180
304
Amonium
Amonium
0
100
18.7
.1
1
[mM]
HORIZONTAL

MONITOR
12
485
107
530
[cells] mg/ml
gbacteria_pa + gbacteria_ec
3
1
11

MONITOR
107
485
203
530
[glucose] mM
gglucose
2
1
11

MONITOR
203
485
288
530
[NH4+] mM
gamonium
2
1
11

MONITOR
288
485
354
530
[O2] mM
goxygen
4
1
11

MONITOR
354
485
446
530
[HCO3-] mM
gbicarbonate
2
1
11

MONITOR
446
485
521
530
[CO2] mM
gdioxide
2
1
11

MONITOR
521
485
596
530
[H2O] mM
gwater
2
1
11

BUTTON
498
211
561
244
Go
go
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

SLIDER
7
143
245
176
energy_maintenance_pa
energy_maintenance_pa
0.0001
0.01
0.0015
0.0001
1
[gC/gCm h]
HORIZONTAL

MONITOR
453
159
546
204
Time [hours]
time_now
2
1
11

SLIDER
7
111
191
144
diffusion-coefficient
diffusion-coefficient
0
1
0.05
0.01
1
NIL
HORIZONTAL

MONITOR
12
395
124
440
NIL
count b-aceptors
1
1
11

MONITOR
649
395
706
440
%_e_c
%_e_c
17
1
11

MONITOR
533
395
590
440
%_e_n
%_e_n
17
1
11

MONITOR
592
395
649
440
yc/c
yc/c
17
1
11

MONITOR
227
395
328
440
c-obtained [mg]
c-obtained
8
1
11

MONITOR
430
396
533
441
c-delivered [mg]
c-delivered
8
1
11

MONITOR
124
395
226
440
n-obtained [mg]
n-obtained
8
1
11

MONITOR
327
395
431
440
n-delivered [mg]
n-delivered
8
1
11

SLIDER
7
175
179
208
rep_pa
rep_pa
1
100
26.0
1
1
NIL
HORIZONTAL

PLOT
220
242
420
392
Biomass
NIL
NIL
0.0
0.5
0.0
0.5
true
true
"" ""
PENS
"pantoea" 1.0 0 -14454117 true "" ""

SLIDER
7
207
211
240
max-time-viability_pa
max-time-viability_pa
0
400
83.0
1
1
[min]
HORIZONTAL

SLIDER
380
126
569
159
Microrganism
Microrganism
10
1000
100.0
1
1
NIL
HORIZONTAL

SLIDER
7
239
179
272
pmax
pmax
1
10
1.0
1
1
NIL
HORIZONTAL

SLIDER
379
27
571
60
total-length-world
total-length-world
0
1000
605.0
1
1
um
HORIZONTAL

SLIDER
379
61
570
94
depth
depth
10
50
50.0
1
1
um
HORIZONTAL

MONITOR
380
159
454
204
Volume [nl]
total-length-world * total-length-world * depth / 1e6
2
1
11

TEXTBOX
403
10
553
28
Real world conditions
11
0.0
1

INPUTBOX
9
304
77
364
file-id
0.0
1
0
Number

INPUTBOX
77
304
142
364
file-ir
0.0
1
0
Number

@#$#@#$#@
## WHAT IS IT?

Here we developed a metamodel in NetLogo to describe the growth of a microbial population consisting of Pantoea. We applied 13 parameters that defined the model and actively changed seven of the parameters to modulate the evolution of the population curve in response to these changes.

## HOW IT WORKS

This model uses many of the implementations in IBM-INDISIM (Ginovart, López et al. 2002), an ABM model which has been successfully used to simulate microbial cultures in situations such as fermentation, multi-species composting, and yeast dynamics in aerobic media and denitrification processes (Ginovart and Cañadas 2008, Prats, Ferrer et al. 2010, Portell, Gras et al. 2014, Banitz, Gras et al. 2015, Araujo Granda, Gras et al. 2016). As in IBM-INDISIM, we use a thermodynamic approach to describe microbial metabolism, e.g., cellular maintenance and mass production, in the ABM for Pantoea. This approach is called the Thermodynamic Equivalent Electron model, TEEM, and it relies on a set of thermo-chemical reactions that account for the Gibbs free energy involved in the overall metabolism, catabolism and anabolism (McCarty 2007). 

## HOW TO USE IT

Changing input parameters (sliders) is possible to test different response of bacteria 
in the model.

## BRIEF DESCRIPTION

Two types of entities exist in our ABM: the bacterium and the square patches where it moves and grows. Each bacterium has its own unique identification number, biomass, metabolism and reproduction parameters, as well as viability and location coordinates, and each one performs the following actions: nutrient uptake, cellular maintenance, biomass synthesis, product generation, and bipartition (reproduction). The time-dependent variables are calculated and updated in each time-step according to a time scale of 0.1 to 1.6 min per time-step. The square patches contain the R2A agar growth medium, and growth medium actions include changes in metabolite concentration and nutrient consumption on each patch. Carbon and nitrogen consumption, as well as their generation and accumulation, are controlled to ensure mass balance. Additionally, the model includes behavior actions that control the overlap of bacteria and the interaction after bipartition.
  The empirical formula C4.17H8O1.75N (Araujo Granda, Gras et al. 2016) was used to describe the elemental composition of Pantoea. The bacterium’s dimensions were assumed to have length and width at 1.1 ± 0.5 µm and 0.55 ± 0.25 µm (mean ± standard deviation), respectively (Kapetas, Ngwenya et al. 2012). The individual mass was deduced from the volume and an assumed density of 1.1 g/cm3 (Gras, Ginovart et al. 2011). A two-dimensional lattice grid was used to describe the environment where the bacteria grow. Each grid-cell represents a volume that can be tuned by changing the world dimensions, depth and length. A volume was calculated using a portion of the wells in the experimental setup with dimensions of 605 µm x 605 µm x depth, where depth ranges from 10 µm to 50 µm.  The initial concentration of glucose and ammonium were estimated from commercial R2A agar composition. 


## CREDITS AND REFERENCES

- Wilensky, U. (1999). NetLogo. http://ccl.northwestern.edu/netlogo/. Center for Connected Learning and Computer-Based Modeling, Northwestern University, Evanston, IL.

-Ginovart, M., et al. (2002). "INDISIM, An Individual-based Discrete Simulation Model to Study Bacterial Cultures." Journal of Theoretical Biology 214(2): 305-319.

-Ginovart, M. and J. C. Cañadas (2008). "INDISIM-YEAST: an individual-based simulator on a website for experimenting and investigating diverse dynamics of yeast populations in liquid media." Journal of Industrial Microbiology and Biotechnology 35(11): 1359-1359.

-Prats, C., et al. (2010). "Individual-based modelling and simulation of microbial processes: yeast fermentation and multi-species composting." Journal of Industrial Microbiology and Biotechnology 16(6): 489-510.

-Portell, X., et al. (2014). "INDISIM-Saccha, an individual-based model to tackle Saccharomyces cerevisiae fermentations." Ecological Modelling 279: 12-23.

-Banitz, T., et al. (2015). "Individual-based modeling of soil organic matter in NetLogo: Transparent, user-friendly, and open." Environmental Modelling & Software 71: 39-45.

-Araujo Granda, P., et al. (2016). "MbT-Tool: An open-access tool based on Thermodynamic Electron Equivalents Model to obtain microbial-metabolic reactions to be used in biotechnological process." Computational and structural biotechnology journal 14: 325-332.

-Araujo Granda, P., et al. (2016). "INDISIM-Paracoccus, an individual-based and thermodynamic model for a denitrifying bacterium." J Theor Biol 403: 45-58.

-Kapetas, L., et al. (2012). "Thermodynamic and Kinetic Controls on Cotransport of Pantoea agglomerans Cells and Zn through Clean and Iron Oxide Coated Sand Columns." Environmental Science & Technology 46(24): 13193-13201.

-Gras, A., et al. (2011). "Individual-based modelling of carbon and nitrogen dynamics in soils: Parameterization and sensitivity analysis of microbial components." Ecological Modelling 222(12): 1998-2010.
@#$#@#$#@
default
true
0
Polygon -7500403 true true 150 5 40 250 150 205 260 250

airplane
true
0
Polygon -7500403 true true 150 0 135 15 120 60 120 105 15 165 15 195 120 180 135 240 105 270 120 285 150 270 180 285 210 270 165 240 180 180 285 195 285 165 180 105 180 60 165 15

arrow
true
0
Polygon -7500403 true true 150 0 0 150 105 150 105 293 195 293 195 150 300 150

bacteria
true
0
Circle -7500403 true true 90 0 120
Circle -7500403 true true 89 179 122
Rectangle -7500403 true true 90 60 210 240

box
false
0
Polygon -7500403 true true 150 285 285 225 285 75 150 135
Polygon -7500403 true true 150 135 15 75 150 15 285 75
Polygon -7500403 true true 15 75 15 225 150 285 150 135
Line -16777216 false 150 285 150 135
Line -16777216 false 150 135 15 75
Line -16777216 false 150 135 285 75

bug
true
0
Circle -7500403 true true 96 182 108
Circle -7500403 true true 110 127 80
Circle -7500403 true true 110 75 80
Line -7500403 true 150 100 80 30
Line -7500403 true 150 100 220 30

butterfly
true
0
Polygon -7500403 true true 150 165 209 199 225 225 225 255 195 270 165 255 150 240
Polygon -7500403 true true 150 165 89 198 75 225 75 255 105 270 135 255 150 240
Polygon -7500403 true true 139 148 100 105 55 90 25 90 10 105 10 135 25 180 40 195 85 194 139 163
Polygon -7500403 true true 162 150 200 105 245 90 275 90 290 105 290 135 275 180 260 195 215 195 162 165
Polygon -16777216 true false 150 255 135 225 120 150 135 120 150 105 165 120 180 150 165 225
Circle -16777216 true false 135 90 30
Line -16777216 false 150 105 195 60
Line -16777216 false 150 105 105 60

car
false
0
Polygon -7500403 true true 300 180 279 164 261 144 240 135 226 132 213 106 203 84 185 63 159 50 135 50 75 60 0 150 0 165 0 225 300 225 300 180
Circle -16777216 true false 180 180 90
Circle -16777216 true false 30 180 90
Polygon -16777216 true false 162 80 132 78 134 135 209 135 194 105 189 96 180 89
Circle -7500403 true true 47 195 58
Circle -7500403 true true 195 195 58

circle
false
0
Circle -7500403 true true 0 0 300

circle 2
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240

cow
false
0
Polygon -7500403 true true 200 193 197 249 179 249 177 196 166 187 140 189 93 191 78 179 72 211 49 209 48 181 37 149 25 120 25 89 45 72 103 84 179 75 198 76 252 64 272 81 293 103 285 121 255 121 242 118 224 167
Polygon -7500403 true true 73 210 86 251 62 249 48 208
Polygon -7500403 true true 25 114 16 195 9 204 23 213 25 200 39 123

cylinder
false
0
Circle -7500403 true true 0 0 300

dot
false
0
Circle -7500403 true true 90 90 120

face happy
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 255 90 239 62 213 47 191 67 179 90 203 109 218 150 225 192 218 210 203 227 181 251 194 236 217 212 240

face neutral
false
0
Circle -7500403 true true 8 7 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Rectangle -16777216 true false 60 195 240 225

face sad
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 168 90 184 62 210 47 232 67 244 90 220 109 205 150 198 192 205 210 220 227 242 251 229 236 206 212 183

fish
false
0
Polygon -1 true false 44 131 21 87 15 86 0 120 15 150 0 180 13 214 20 212 45 166
Polygon -1 true false 135 195 119 235 95 218 76 210 46 204 60 165
Polygon -1 true false 75 45 83 77 71 103 86 114 166 78 135 60
Polygon -7500403 true true 30 136 151 77 226 81 280 119 292 146 292 160 287 170 270 195 195 210 151 212 30 166
Circle -16777216 true false 215 106 30

flag
false
0
Rectangle -7500403 true true 60 15 75 300
Polygon -7500403 true true 90 150 270 90 90 30
Line -7500403 true 75 135 90 135
Line -7500403 true 75 45 90 45

flower
false
0
Polygon -10899396 true false 135 120 165 165 180 210 180 240 150 300 165 300 195 240 195 195 165 135
Circle -7500403 true true 85 132 38
Circle -7500403 true true 130 147 38
Circle -7500403 true true 192 85 38
Circle -7500403 true true 85 40 38
Circle -7500403 true true 177 40 38
Circle -7500403 true true 177 132 38
Circle -7500403 true true 70 85 38
Circle -7500403 true true 130 25 38
Circle -7500403 true true 96 51 108
Circle -16777216 true false 113 68 74
Polygon -10899396 true false 189 233 219 188 249 173 279 188 234 218
Polygon -10899396 true false 180 255 150 210 105 210 75 240 135 240

house
false
0
Rectangle -7500403 true true 45 120 255 285
Rectangle -16777216 true false 120 210 180 285
Polygon -7500403 true true 15 120 150 15 285 120
Line -16777216 false 30 120 270 120

leaf
false
0
Polygon -7500403 true true 150 210 135 195 120 210 60 210 30 195 60 180 60 165 15 135 30 120 15 105 40 104 45 90 60 90 90 105 105 120 120 120 105 60 120 60 135 30 150 15 165 30 180 60 195 60 180 120 195 120 210 105 240 90 255 90 263 104 285 105 270 120 285 135 240 165 240 180 270 195 240 210 180 210 165 195
Polygon -7500403 true true 135 195 135 240 120 255 105 255 105 285 135 285 165 240 165 195

line
true
0
Line -7500403 true 150 0 150 300

line half
true
0
Line -7500403 true 150 0 150 150

pentagon
false
0
Polygon -7500403 true true 150 15 15 120 60 285 240 285 285 120

person
false
0
Circle -7500403 true true 110 5 80
Polygon -7500403 true true 105 90 120 195 90 285 105 300 135 300 150 225 165 300 195 300 210 285 180 195 195 90
Rectangle -7500403 true true 127 79 172 94
Polygon -7500403 true true 195 90 240 150 225 180 165 105
Polygon -7500403 true true 105 90 60 150 75 180 135 105

plant
false
0
Rectangle -7500403 true true 135 90 165 300
Polygon -7500403 true true 135 255 90 210 45 195 75 255 135 285
Polygon -7500403 true true 165 255 210 210 255 195 225 255 165 285
Polygon -7500403 true true 135 180 90 135 45 120 75 180 135 210
Polygon -7500403 true true 165 180 165 210 225 180 255 120 210 135
Polygon -7500403 true true 135 105 90 60 45 45 75 105 135 135
Polygon -7500403 true true 165 105 165 135 225 105 255 45 210 60
Polygon -7500403 true true 135 90 120 45 150 15 180 45 165 90

sheep
false
15
Circle -1 true true 203 65 88
Circle -1 true true 70 65 162
Circle -1 true true 150 105 120
Polygon -7500403 true false 218 120 240 165 255 165 278 120
Circle -7500403 true false 214 72 67
Rectangle -1 true true 164 223 179 298
Polygon -1 true true 45 285 30 285 30 240 15 195 45 210
Circle -1 true true 3 83 150
Rectangle -1 true true 65 221 80 296
Polygon -1 true true 195 285 210 285 210 240 240 210 195 210
Polygon -7500403 true false 276 85 285 105 302 99 294 83
Polygon -7500403 true false 219 85 210 105 193 99 201 83

square
false
0
Rectangle -7500403 true true 30 30 270 270

square 2
false
0
Rectangle -7500403 true true 30 30 270 270
Rectangle -16777216 true false 60 60 240 240

star
false
0
Polygon -7500403 true true 151 1 185 108 298 108 207 175 242 282 151 216 59 282 94 175 3 108 116 108

target
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240
Circle -7500403 true true 60 60 180
Circle -16777216 true false 90 90 120
Circle -7500403 true true 120 120 60

tree
false
0
Circle -7500403 true true 118 3 94
Rectangle -6459832 true false 120 195 180 300
Circle -7500403 true true 65 21 108
Circle -7500403 true true 116 41 127
Circle -7500403 true true 45 90 120
Circle -7500403 true true 104 74 152

triangle
false
0
Polygon -7500403 true true 150 30 15 255 285 255

triangle 2
false
0
Polygon -7500403 true true 150 30 15 255 285 255
Polygon -16777216 true false 151 99 225 223 75 224

truck
false
0
Rectangle -7500403 true true 4 45 195 187
Polygon -7500403 true true 296 193 296 150 259 134 244 104 208 104 207 194
Rectangle -1 true false 195 60 195 105
Polygon -16777216 true false 238 112 252 141 219 141 218 112
Circle -16777216 true false 234 174 42
Rectangle -7500403 true true 181 185 214 194
Circle -16777216 true false 144 174 42
Circle -16777216 true false 24 174 42
Circle -7500403 false true 24 174 42
Circle -7500403 false true 144 174 42
Circle -7500403 false true 234 174 42

turtle
true
0
Polygon -10899396 true false 215 204 240 233 246 254 228 266 215 252 193 210
Polygon -10899396 true false 195 90 225 75 245 75 260 89 269 108 261 124 240 105 225 105 210 105
Polygon -10899396 true false 105 90 75 75 55 75 40 89 31 108 39 124 60 105 75 105 90 105
Polygon -10899396 true false 132 85 134 64 107 51 108 17 150 2 192 18 192 52 169 65 172 87
Polygon -10899396 true false 85 204 60 233 54 254 72 266 85 252 107 210
Polygon -7500403 true true 119 75 179 75 209 101 224 135 220 225 175 261 128 261 81 224 74 135 88 99

wheel
false
0
Circle -7500403 true true 3 3 294
Circle -16777216 true false 30 30 240
Line -7500403 true 150 285 150 15
Line -7500403 true 15 150 285 150
Circle -7500403 true true 120 120 60
Line -7500403 true 216 40 79 269
Line -7500403 true 40 84 269 221
Line -7500403 true 40 216 269 79
Line -7500403 true 84 40 221 269

wolf
false
0
Polygon -16777216 true false 253 133 245 131 245 133
Polygon -7500403 true true 2 194 13 197 30 191 38 193 38 205 20 226 20 257 27 265 38 266 40 260 31 253 31 230 60 206 68 198 75 209 66 228 65 243 82 261 84 268 100 267 103 261 77 239 79 231 100 207 98 196 119 201 143 202 160 195 166 210 172 213 173 238 167 251 160 248 154 265 169 264 178 247 186 240 198 260 200 271 217 271 219 262 207 258 195 230 192 198 210 184 227 164 242 144 259 145 284 151 277 141 293 140 299 134 297 127 273 119 270 105
Polygon -7500403 true true -1 195 14 180 36 166 40 153 53 140 82 131 134 133 159 126 188 115 227 108 236 102 238 98 268 86 269 92 281 87 269 103 269 113

x
false
0
Polygon -7500403 true true 270 75 225 30 30 225 75 270
Polygon -7500403 true true 30 75 75 30 270 225 225 270
@#$#@#$#@
NetLogo 6.0.4
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
<experiments>
  <experiment name="pantoea" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <timeLimit steps="50"/>
    <metric>time_now</metric>
    <enumeratedValueSet variable="Microrganism">
      <value value="100"/>
      <value value="150"/>
    </enumeratedValueSet>
  </experiment>
</experiments>
@#$#@#$#@
@#$#@#$#@
default
0.0
-0.2 0 0.0 1.0
0.0 1 1.0 0.0
0.2 0 0.0 1.0
link direction
true
0
Line -7500403 true 150 150 90 180
Line -7500403 true 150 150 210 180
@#$#@#$#@
0
@#$#@#$#@
