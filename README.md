# Transfer-Learning-of-Stochastic-Kriging-for-Individualized-Prediction_T-PAMI-

To replicate the simulation results for the boxplot in **Fig. 5** of the manuscript, simply run the `box_pd.R` script. This one-click execution integrates all benchmark methods and the proposed approach.

## Required R Packages

Make sure the following packages are installed:
- `Matrix`
- `nloptr`
- `minna`
- `optimx`
- `nlme`
- `mvtnorm`
- `MASS`

## Interpret the Result

The simulation is based on **Signal Setting I**. The `box_pd.R` script uses the `source()` function to run individual files that implement the proposed method and each benchmark.  

- The RMSE results are saved in a `100 x 6` array named `RMSE`, where:  
  - Each **row** corresponds to a single repetition of prediction and RMSE calculation.  
  - Each **column** corresponds to one of the following methods:  
    1. Proposed  
    2. Unpenalized  
    3. Minimal Transfer  
    4. i.i.d. MGP  
    5. Individualized SK
    6. i.i.d Enforeced 

The default number of repetitions is **100** (`rep_time = 100` in line 2 of `box_pd.R`). The simulation takes approximately **15 hours** to complete across all methods.

## Reproduce Other RMSE Results

To replicate the RMSE results for **Signal Setting II** (shown in Fig. 8 of the manuscript):
1. Open `TrainData.R`
2. In **line 4**, change the sourced file from `"fun1.R"` to `"fun2.R"`

This will switch the signal generation from Setting I (Eq. 19) to Setting II (Eq. 21), which uses segmentally collected data.

## File Structure and Function Summary

The `demo` folder includes the following **9 files**:

### Primary Scripts

1. **box_pd.R**  
   - Integrates RMSE results from proposed and benchmark methods  
   - Generates the RMSE boxplot for comparison  

2. **TrainData.R**  
   - Sources either `fun1.R` or `fun2.R` to create heterogeneous replications  
   - Computes sample means and noise scaling factors \( R_{i,j} \)

### Method Implementations

3. **MGP_stoch.R** — Implements the **Proposed** method  
4. **MGP_stoch(unpenalized).R** — Implements the **Unpenalized** method  
5. **Minimal Transfer.R** — Implements the **Minimal Transfer** method  
6. **Multioutput with identity.R** — Implements the **i.i.d. MGP** method  
7. **Single SK.R** — Implements the **Individualized SK** method  
8. **FGP_19.R** — Implements the **iid Enforced** method  

### Signal Generation Files

9. **fun1.R**  
   - Constructs signals for **Setting I** (Eq. 19) with randomly collected data  

10. **fun2.R**  
   - Constructs signals for **Setting II** (Eq. 21) with segmentally collected data  
