# matom_shell

Macro-atom shell data and related scripts

### Folders

* scripts_useful: useful scripts not directly related
* m-shell: full dataset
* c-test: Carbon dataset 




## Data Formats

### LINELIST

Z, ION, LINES
LINE INDEX, LOWER LEVEL, HIGHER LEVEL, A VALUE

     6       4      18
       1       1       2   2.6300000e+08
       2       1       3   2.6500000e+08
       3       1       5   4.5500001e+09

### COLLLIST

Bound-bound collisions

electrons cooler than excitation energy? take low temperature limie


Z, ION, COLLLINES
LINE INDEX, LOWER LEVEL, HIGHER LEVEL, UPSILON
      6       4      14
       1       1       2       2.8340001
       2       1       3       5.6020002
       3       1       4      0.65560001



### PHIXS
Z, ION UPPER, LEVEL UPPER, ION LOWER, LEVEL LOWER, NUMBER OF RECORDS
EV ABOVE THRESHOLD, SIGMA IN MEGABARNS
6           5       1           4       1      50



### AUGER
Z ION NJUMPS
INDEX, LOWER ION LEVEL, UPPER ION LEVEL,??,A VALUE
	26	22	350
     1     514       1      4.75500  7.97000e+12


STATE 514 of Fe 22 autoionizes to ground state of Fe 23, with A_r = 7.97000e+12
Some have flag due in upper ion column, no info about target ions- -1 is a flag, where
we have a Total Auger ionization state, but no idea about where we go -> make a kpkt.

Autoionization can be way higher than radiative de-excitation so need the rate
to get line strength right.

When we autoionize we have a probability of making a kpkt, with probability equal to the 
energy which flows into each process.


###input_kedge.txt
assumption- all ion has K-shell occupied, so don't worry if upper levels occupied

NRECORDS
Z, ION, principal quantum number, l quantum number, E0 (eV), fit coefficients 

Next lines, what to do with jump!
HOW MANY K VACANCY STATES I COULD JUMP TO (2 TARGET LEVELS)
NEXT LINES ARE WHICH LEVELS, TO GO TO
-1 -1 -1 -1 --- KAASTRA AND MAWE -- FOUR NUMBERS SAY SINGLY< DOUBLEY< TRIPLY - EFFECTIVE FLUORESCENCE YIELDS.



```
97  
6 	4	1	0 	3.522e+02 	8.412e+01 	8.111e+01 	7.459e+01 	1.428e+00 	0.000e+00
2
2
3
-1 -1 -1 -1
```
