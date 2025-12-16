Epidemiological Model: Intervention Timing vs Coverage
================
Saud
2025-12-16

- [Model Definition](#model-definition)
- [Simulation Execution](#simulation-execution)
- [Visualization](#visualization)

``` r
# odin is used to write and compile differential equations
library(odin)
# tidyverse includes ggplot2 (for plotting) and dplyr (for data manipulation)
library(tidyverse)
```

# Model Definition

``` r
model_code <- "
  ## --- 1. Transmission Rates (Force of Infection) ---
  # The 'Force of Infection' (lambda) represents the rate at which susceptible individuals become exposed.
  # It depends on the number of infectious people (Ia, Iu_c) and the contact rates (beta).
  
  # Lambda for Children: Risks from Adults (Ia) and other Children (T_c + Iu_c)
  lambda_c <- (beta_ac * (Ia / total_adults)) + (beta_cc * ((T_c + Iu_c) / total_children)) 
  
  # Lambda for Adults: Risks from Children (Iu_c + T_c) and other Adults (Ia)
  lambda_a <- (beta_ca * (Iu_c + T_c) / total_children) + (beta_aa * (Ia / total_adults)) 

  ## --- 2. Parameters (User-defined) ---
  # 'user()' tells odin that these values can be changed later when we run the simulation
  
  theta <- user(1 / 7)    # Latency rate: Children take ~7 days to become infectious (E -> I)
  theta_a <- user(1 / 10) # Latency rate: Adults take ~10 days to become infectious (E -> I)

  # Transmission coefficients (Probability of transmission per contact)
  beta_cc <- user(0.8)    # Child-to-Child contact
  beta_ac <- user(0.5)    # Adult-to-Child contact
  beta_ca <- beta_ac      # Child-to-Adult (assumed symmetric)
  beta_aa <- user(0.6)    # Adult-to-Adult contact

  ## --- 3. Treatment Strategy Logic ---
  # We are modeling a trade-off:
  # - Cheap Medicine = High Coverage (active >= 0.5) but Low Efficacy (alpha = 2)
  # - Expensive Medicine = Low Coverage (active < 0.5) but High Efficacy (alpha = 4)
  
  alpha_low_cov <- user(2) 
  alpha_high_cov <- user(4) 
  cov_active <- user(0.6) # The % of children reached by the treatment
  
  # IF coverage is low (< 50%), use the high efficacy rate (Expensive Medicine)
  # ELSE use the low efficacy rate (Cheap Medicine)
  alpha <- if (cov_active < 0.5) alpha_high_cov else alpha_low_cov

  # Pulse Treatment Logic:
  # Treatment is only active for ONE DAY, starting at 't_start'
  t_start <- user(14) 
  cov <- if (t >= t_start && t <= (t_start + 1)) cov_active else 0 

  ## --- 4. Initial Conditions ---
  # Sets the starting population numbers at Day 0
  initial(S_c) <- 100   # 100 Healthy children
  initial(E_c) <- 0     # 0 Exposed
  initial(T_c) <- 0     # 0 Treated
  initial(Iu_c) <- 1    # 1 Infectious child (Patient Zero)
  initial(R_c) <- 0     # 0 Recovered
  
  initial(S_a) <- 20    # 20 Healthy adults
  initial(E_a) <- 0     # 0 Exposed
  initial(Ia) <- 0      # 0 Infectious

  ## --- 5. Differential Equations: Children ---
  # deriv() defines the rate of change over time (dS/dt, dE/dt, etc.)
  
  deriv(S_c) <- -lambda_c * S_c  # Susceptibles decrease as they get infected
  deriv(E_c) <- lambda_c * S_c - theta * E_c # Exposed increase by infection, decrease by becoming infectious
  
  # Treated Children (T_c):
  # A fraction 'cov' of Exposed children get treatment. They recover at rate 'alpha'.
  deriv(T_c) <- cov * theta * E_c - alpha * T_c 
  
  # Untreated Infectious Children (Iu_c):
  # The remaining fraction (1 - cov) do not get treatment.
  deriv(Iu_c) <- (1 - cov) * theta * E_c 
  
  deriv(R_c) <- alpha * T_c # Recovered children accumulate here

  ## --- 6. Differential Equations: Adults ---
  deriv(S_a) <- -lambda_a * S_a
  deriv(E_a) <- lambda_a * S_a - theta_a * E_a
  deriv(Ia) <- theta_a * E_a 

  ## --- 7. Output Calculations ---
  # Helper variables to track total population sizes to keep math consistent
  total_children <- S_c + E_c + T_c + Iu_c + R_c
  total_adults <- S_a + E_a + Ia
  
  # We want to track the number of currently active infections (Treated + Untreated)
  output(total_infected_children) <- T_c + Iu_c 
"

# Compile the model logic into C code (this makes it run very fast)
model <- odin::odin(model_code)
```

    ## ── R CMD INSTALL ───────────────────────────────────────────────────────────────
    ##   ─  installing *source* package 'odin9a2c293f' ... (335ms)
    ##      ** using staged installation
    ##      ** libs
    ##      using C compiler: 'gcc.exe (GCC) 13.3.0'
    ##      gcc  -I"C:/PROGRA~1/R/R-44~1.2/include" -DNDEBUG     -I"C:/RBuildTools/4.4/x86_64-w64-mingw32.static.posix/include"     -O2 -Wall -gdwarf-2 -mfpmath=sse -msse2 -mstackrealign  -UNDEBUG -Wall -pedantic -g -O0 -c odin.c -o odin.o
    ##      odin.c: In function 'odin_metadata':
    ##      odin.c:170:18: warning: unused variable 'internal' [-Wunused-variable]
    ##      170 |   odin_internal *internal = odin_get_internal(internal_p, 1);
    ##          |                  ^~~~~~~~
    ##    odin.c: In function 'odin_output_dde':
    ##    odin.c:258:18: warning: unused variable 'internal' [-Wunused-variable]
    ##      258 |   odin_internal *internal = (odin_internal*) internal_p;
    ##          |                  ^~~~~~~~
    ##      odin.c: In function 'user_get_scalar_int':
    ##    odin.c:302:47: warning: format '%d' expects argument of type 'int', but argument 2 has type 'const char *' [-Wformat=]
    ##      302 |       Rf_error("Expected scalar integer for '%d'", name);
    ##          |                                              ~^    ~~~~
    ##          |                                               |    |
    ##          |                                               int  const char *
    ##          |                                              %s
    ##      gcc  -I"C:/PROGRA~1/R/R-44~1.2/include" -DNDEBUG     -I"C:/RBuildTools/4.4/x86_64-w64-mingw32.static.posix/include"     -O2 -Wall -gdwarf-2 -mfpmath=sse -msse2 -mstackrealign  -UNDEBUG -Wall -pedantic -g -O0 -c registration.c -o registration.o
    ##      gcc -shared -static-libgcc -o odin9a2c293f.dll tmp.def odin.o registration.o -LC:/RBuildTools/4.4/x86_64-w64-mingw32.static.posix/lib/x64 -LC:/RBuildTools/4.4/x86_64-w64-mingw32.static.posix/lib -LC:/PROGRA~1/R/R-44~1.2/bin/x64 -lR
    ##      installing to C:/Users/SAUDHA~1/AppData/Local/Temp/RtmpiSDdpn/devtools_install_585c6f86773e/00LOCK-file585c5964536d/00new/odin9a2c293f/libs/x64
    ##   ─  DONE (odin9a2c293f)
    ## 

# Simulation Execution

``` r
# 1. Define the ranges we want to test
t_start_range <- seq(1, 50, by = 1)      # Try every start day from 1 to 50
cov_active_range <- seq(0.1, 1, by = 0.1) # Try coverage from 10% to 100%

# 2. Create the Parameter Grid
# expand.grid creates a table with every possible combination of start_day and coverage
param_grid <- expand.grid(t_start = t_start_range, cov_active = cov_active_range)
times <- seq(0, 50, by = 1) # Run each simulation for 50 days

# 3. Run Loop (Using lapply for efficiency)
# We iterate through every row (i) in our parameter grid
results_list <- lapply(1:nrow(param_grid), function(i) {
  
  # Extract the specific values for this run
  current_t <- param_grid$t_start[i]
  current_cov <- param_grid$cov_active[i]
  
  # List all parameters to send to the model
  params <- list(
    t_start = current_t,
    cov_active = current_cov,
    theta = 1 / 7,
    theta_a = 1 / 10,
    beta_cc = 0.8,
    beta_ac = 0.5,
    beta_aa = 0.6,
    alpha_low_cov = 2, 
    alpha_high_cov = 4
  )
  
  # Initialize the model with these specific parameters
  mod_instance <- model$new(user = params)
  
  # Run the simulation over time
  sim <- mod_instance$run(times)
  
  # Extract the result we care about: The HIGHEST number of infected children (Peak)
  peak_infected <- max(sim[, "total_infected_children"])
  
  # Return a clean data frame row with the inputs and the result
  data.frame(
    t_start = current_t,
    cov_active = current_cov,
    peak_infected_children = peak_infected
  )
})

# Combine the list of results into one big data frame
results <- do.call(rbind, results_list)

# 4. Categorize Strategy for Plotting
# We add a label to distinguish the two strategy types based on coverage
results$medicine_type <- ifelse(results$cov_active < 0.5, 
                                "Expensive Medicine (Low Coverage)", 
                                "Cheap Medicine (High Coverage)")
```

# Visualization

``` r
ggplot(results, aes(x = t_start, y = peak_infected_children, color = medicine_type, group = medicine_type)) +
  geom_line(size = 1) + 
  labs(title = "Peak Infected Children vs. Start of Intervention",
       subtitle = "Comparison of Medicine Strategies",
       x = "Start Day of Intervention (t_start)",
       y = "Peak Infected Children (Count)",
       color = "Strategy") +
  # Custom colors: Blue for high coverage, Red for low coverage
  scale_color_manual(values = c("Cheap Medicine (High Coverage)" = "blue", 
                                "Expensive Medicine (Low Coverage)" = "red")) +
  theme_minimal() +
  theme(legend.position = "top")
```

![](IFDM_project_files/figure-gfm/unnamed-chunk-4-1.png)<!-- -->

``` r
ggplot(results, aes(x = t_start, y = cov_active, fill = peak_infected_children)) +
  geom_tile() + 
  scale_fill_gradientn(
    colors = c("green", "yellow", "orange"),
    name = "Total Infections"
  ) +
  labs(
    title = "Total Infections vs. Treatment Start and Coverage",
    x = "Treatment Start Day (t_start)",
    y = "Coverage Rate (cov_active)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5), # Center the title
    legend.position = "right"
  )
```

![](IFDM_project_files/figure-gfm/unnamed-chunk-4-2.png)<!-- -->
