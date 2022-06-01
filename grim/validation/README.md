
# graph-imputation-match
  
Graph-based imputation and matching.


| Subdirectory    | Description |
| :-------------- | :---------- |
| g2gl            | conversion of genotypes to genotype-list (deprecated) |
| haplogic        | CT validation data used to evaluate haplogic |
| simulation      | simulated datasets |
| wmda            | WMDA validation datasets |



# Explanation of `priority` parameter.

### delta:
True prior SR (race1, race1) +=delta and (race2, race2)+=delta
           
### alpha: 
only taking the right combination- prior(race1, race2)+=alpha and prior(race2, race1)+=alpha

### eta: 
allowing for all Multi-Race combination- for each i,j - prior(i,j)+=eta. if eta=1- is like with no prior

### beta:
allows for all single race populations- for each i - prior(i,i)+=beta

### gamma:
the whole column and row of race1 and race2 get prior- gamma



examples (pops= AAFA, NAMER, CARB):

1. priority = {'alpha':0, 'eta':0, 'beta':0, 'gamma':0, 'delta':1}
   
   race1=AAFA, race2=NAMER
            AAFA NAMER CARB
    AAFA      1     0   0
    NAMER     0     1   0
    CARB      0     0   0
    
    race1=AAFA, race2=AAFA
            AAFA NAMER CARB
    AAFA      1     0   0
    NAMER     0     0   0
    CARB      0     0   0
    
2. priority = {'alpha':0, 'eta':1, 'beta':0, 'gamma':0, 'delta':0}

   race1=AAFA, race2=NAMER / race1=AAFA, race2=AAFA
            AAFA NAMER CARB
    AAFA      1     1   1
    NAMER     1     1   1
    CARB      1     1   1
    
3. priority = {'alpha':1, 'eta':0, 'beta':0, 'gamma':0, ' 'delta':0}

    race1=AAFA, race2=NAMER
            AAFA NAMER CARB
    AAFA      0     1   0
    NAMER     1     0   0
    CARB      0     0   0
    
    race1=AAFA, race2=AAFA
            AAFA NAMER CARB
    AAFA      1     0   0
    NAMER     0     0   0
    CARB      0     0   0
    
4.priority = {'alpha':0, 'eta':0, 'beta':1, 'gamma':0, 'delta':0}

   race1=AAFA, race2=NAMER / race1=AAFA, race2=AAFA
            AAFA NAMER CARB
    AAFA      1     0   0
    NAMER     0     1   0
    CARB      0     0   1
    
5. priority = {'alpha':0, 'eta':0, 'beta':0, 'gamma':1, 'delta':0}

    race1=AAFA, race2=NAMER
            AAFA NAMER CARB
    AAFA     1    1    1
    NAMER    1    1    1
    CARB     1    1    0
    
    race1=AAFA, race2=AAFA
            AAFA NAMER CARB
    AAFA      1     1   1
    NAMER     1     0   0
    CARB      1     0   0

6. priority = {'alpha': 0.4999999, 'eta': 0, 'beta': 1e-07, 'gamma': 1e-07, 'delta': 0.4999999}


    race1=AAFA, race2=NAMER
             AAFA        NAMER        CARB
    AAFA     0.5000001   0.5          1e-07
    NAMER    0.5         0.5000001    1e-07
    CARB     1e-07       1e-07        1e-07

    race1=AAFA, race2=AAFA
             AAFA     NAMER    CARB
    AAFA     1        1e-07    1e-07
    NAMER    1e-07    1e-07    0
    CARB     1e-07    0        1e-07

