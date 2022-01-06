#### parameters ####

# plant-only and single virus
params_def1 <- c(a_n_lo = 1.1e-6,
                 a_n_hi = 5.6e-5,
                 a_p_lo = 1.6e-7,
                 a_p_hi = 8.2e-6,
                 u_n = 2.9e-3,
                 u_p = 1.4e-3,
                 k_n = 4.9e-5,
                 k_p = 3.0e-5,
                 Qmin_n = 1.1e-3,
                 Qmin_p = 7.4e-5,
                 q_n = 1.1e-3,
                 q_p = 7.4e-5,
                 z_n = 1.7e-18,
                 z_p = 2.6e-19)

# two viruses
params_def2 <- c(a_n_lo = 1.1e-6,
                 a_n_hi = 5.6e-5,
                 a_p_lo = 1.6e-7,
                 a_p_hi = 8.2e-6,
                 u_n = 2.9e-3,
                 u_p = 1.4e-3,
                 k_n = 4.9e-5,
                 k_p = 3.0e-5,
                 Qmin_n = 1.1e-3,
                 Qmin_p = 7.4e-5,
                 q_nb = 1.1e-3,
                 q_pb = 7.4e-5,
                 q_nc = 1.1e-3,
                 q_pc = 7.4e-5,
                 z_nb = 1.6e-18,
                 z_pb = 2.6e-19,
                 z_nc = 1.7e-18,
                 z_pc = 2.6e-19,
                 m = 0.000860,
                 g = 0.136,
                 c_b = 0.00544,
                 c_c = 0.0000175,
                 r_b = 0.0961,
                 r_c = 0.221)


#### initial values ####

# initial values
H0 <- 1e-2
Q0_n <- as.numeric(params_def1["Qmin_n"]) * 10
Q0_p <- as.numeric(params_def1["Qmin_p"]) * 10
R0_n_lo <- as.numeric(params_def1["a_n_hi"]) * 7 + as.numeric(params_def1["a_n_lo"]) * 3
R0_n_hi <- as.numeric(params_def1["a_n_hi"]) * 10
R0_p_lo <- as.numeric(params_def1["a_p_hi"]) * 7 + as.numeric(params_def1["a_p_lo"]) * 3
R0_p_hi <- as.numeric(params_def1["a_p_hi"]) * 10
V0 <- 5e5

# plant-only and single virus
init_def1 <- c(R_n_low = R0_n_lo,
               R_p_low = R0_p_lo,
               Q_n_low = Q0_n,
               Q_p_low = Q0_p,
               H_low = H0,
               R_n_n = R0_n_hi,
               R_p_n = R0_p_lo,
               Q_n_n = Q0_n,
               Q_p_n = Q0_p,
               H_n = H0,
               R_n_p = R0_n_lo,
               R_p_p = R0_p_hi,
               Q_n_p = Q0_n,
               Q_p_p = Q0_p,
               H_p = H0,
               R_n_np = R0_n_hi,
               R_p_np = R0_p_hi,
               Q_n_np = Q0_n,
               Q_p_np = Q0_p,
               H_np = H0)

init_def2 <- c(R_n_low = R0_n_lo,
               R_p_low = R0_p_lo,
               Q_n_low = Q0_n,
               Q_p_low = Q0_p,
               H_low = H0,
               V_b_low = 0, 
               V_c_low = 0,
               R_n_n = R0_n_hi,
               R_p_n = R0_p_lo,
               Q_n_n = Q0_n,
               Q_p_n = Q0_p,
               H_n = H0,
               V_b_n = 0, 
               V_c_n = 0,
               R_n_p = R0_n_lo,
               R_p_p = R0_p_hi,
               Q_n_p = Q0_n,
               Q_p_p = Q0_p,
               H_p = H0,
               V_b_p = 0, 
               V_c_p = 0,
               R_n_np = R0_n_hi,
               R_p_np = R0_p_hi,
               Q_n_np = Q0_n,
               Q_p_np = Q0_p,
               H_np = H0,
               V_b_np = 0, 
               V_c_np = 0)


#### plant only model ####

plant_model = function (t, yy, parms) { 
  
  with(as.list(c(yy, parms)), {
    
    # model
    dR_n_low <- a_n_lo - (u_n * R_n_low * H_low) / (R_n_low + k_n)
    dR_p_low <- a_p_lo - (u_p * R_p_low * H_low) / (R_p_low + k_p)
    dQ_n_low <- (u_n * R_n_low) / (R_n_low + k_n) - min((1 - Qmin_n / Q_n_low), (1 - Qmin_p / Q_p_low)) * g * Q_n_low
    dQ_p_low <- (u_p * R_p_low) / (R_p_low + k_p) - min((1 - Qmin_n / Q_n_low), (1 - Qmin_p / Q_p_low)) * g * Q_p_low
    dH_low <- min((1 - Qmin_n / Q_n_low), (1 - Qmin_p / Q_p_low)) * g * H_low - m * H_low
    
    dR_n_n <- a_n_hi - (u_n * R_n_n * H_n) / (R_n_n + k_n)
    dR_p_n <- a_p_lo - (u_p * R_p_n * H_n) / (R_p_n + k_p)
    dQ_n_n <- (u_n * R_n_n) / (R_n_n + k_n) - min((1 - Qmin_n / Q_n_n), (1 - Qmin_p / Q_p_n)) * g * Q_n_n
    dQ_p_n <- (u_p * R_p_n) / (R_p_n + k_p) - min((1 - Qmin_n / Q_n_n), (1 - Qmin_p / Q_p_n)) * g * Q_p_n
    dH_n <- min((1 - Qmin_n / Q_n_n), (1 - Qmin_p / Q_p_n)) * g * H_n - m * H_n
    
    dR_n_p <- a_n_lo - (u_n * R_n_p * H_p) / (R_n_p + k_n)
    dR_p_p <- a_p_hi - (u_p * R_p_p * H_p) / (R_p_p + k_p)
    dQ_n_p <- (u_n * R_n_p) / (R_n_p + k_n) - min((1 - Qmin_n / Q_n_p), (1 - Qmin_p / Q_p_p)) * g * Q_n_p
    dQ_p_p <- (u_p * R_p_p) / (R_p_p + k_p) - min((1 - Qmin_n / Q_n_p), (1 - Qmin_p / Q_p_p)) * g * Q_p_p
    dH_p <- min((1 - Qmin_n / Q_n_p), (1 - Qmin_p / Q_p_p)) * g * H_p - m * H_p
    
    dR_n_np <- a_n_hi - (u_n * R_n_np * H_np) / (R_n_np + k_n)
    dR_p_np <- a_p_hi - (u_p * R_p_np * H_np) / (R_p_np + k_p)
    dQ_n_np <- (u_n * R_n_np) / (R_n_np + k_n) - min((1 - Qmin_n / Q_n_np), (1 - Qmin_p / Q_p_np)) * g * Q_n_np
    dQ_p_np <- (u_p * R_p_np) / (R_p_np + k_p) - min((1 - Qmin_n / Q_n_np), (1 - Qmin_p / Q_p_np)) * g * Q_p_np
    dH_np <- min((1 - Qmin_n / Q_n_np), (1 - Qmin_p / Q_p_np)) * g * H_np - m * H_np
    
    return(list(c(dR_n_low, dR_p_low, dQ_n_low, dQ_p_low, dH_low,
                  dR_n_n, dR_p_n, dQ_n_n, dQ_p_n, dH_n,
                  dR_n_p, dR_p_p, dQ_n_p, dQ_p_p, dH_p,
                  dR_n_np, dR_p_np, dQ_n_np, dQ_p_np, dH_np)))
  })
}

plant_model_format <- function(mod_in){
  
  mod_out <- mod_in %>%
    as_tibble() %>%
    mutate(across(everything(), as.double)) %>%
    pivot_longer(cols = -time,
                 names_to = "variable",
                 values_to = "value") %>%
    mutate(nutrient = case_when(str_ends(variable, "low") == T ~ "low",
                                str_ends(variable, "np") == T ~ "+N+P",
                                str_ends(variable, "n") == T ~ "+N",
                                str_ends(variable, "p") == T ~ "+P") %>%
             fct_relevel("low", "+N", "+P"),
           variable2 = case_when(str_starts(variable, "R_n") == T ~ "R_n",
                                 str_starts(variable, "R_p") == T ~ "R_p",
                                 str_starts(variable, "Q_n") == T ~ "Q_n",
                                 str_starts(variable, "Q_p") == T ~ "Q_p",
                                 str_starts(variable, "H") == T ~ "H"))
  
  return(mod_out)
  
}


#### single virus model ####

virus1_model <- function (t, yy, parms) { 
  
  with(as.list(c(yy, parms)), {
  
  # model
  dR_n_low <- a_n_lo - (u_n * R_n_low * H_low) / (R_n_low + k_n);
  dR_p_low <- a_p_lo - (u_p * R_p_low * H_low) / (R_p_low + k_p);
  dQ_n_low <- (u_n * R_n_low) / (R_n_low + k_n) - min((1 - Qmin_n / Q_n_low), (1 - Qmin_p / Q_p_low)) * g * Q_n_low - min((1 - q_n / Q_n_low), (1 - q_p / Q_p_low)) * z_n * r * V_low
  dQ_p_low <- (u_p * R_p_low) / (R_p_low + k_p) - min((1 - Qmin_n / Q_n_low), (1 - Qmin_p / Q_p_low)) * g * Q_p_low - min((1 - q_n / Q_n_low), (1 - q_p / Q_p_low)) * z_p * r * V_low
  dH_low <- min((1 - Qmin_n / Q_n_low), (1 - Qmin_p / Q_p_low)) * g * H_low - m * H_low
  dV_low <- min((1 - q_n / Q_n_low), (1 - q_p / Q_p_low)) * r * V_low - c * V_low
  
  dR_n_n <- a_n_hi - (u_n * R_n_n * H_n) / (R_n_n + k_n);
  dR_p_n <- a_p_lo - (u_p * R_p_n * H_n) / (R_p_n + k_p);
  dQ_n_n <- (u_n * R_n_n) / (R_n_n + k_n) - min((1 - Qmin_n / Q_n_n), (1 - Qmin_p / Q_p_n)) * g * Q_n_n - min((1 - q_n / Q_n_n), (1 - q_p / Q_p_n)) * z_n * r * V_n
  dQ_p_n <- (u_p * R_p_n) / (R_p_n + k_p) - min((1 - Qmin_n / Q_n_n), (1 - Qmin_p / Q_p_n)) * g * Q_p_n - min((1 - q_n / Q_n_n), (1 - q_p / Q_p_n)) * z_p * r * V_n
  dH_n <- min((1 - Qmin_n / Q_n_n), (1 - Qmin_p / Q_p_n)) * g * H_n - m * H_n
  dV_n <- min((1 - q_n / Q_n_n), (1 - q_p / Q_p_n)) * r * V_n - c * V_n
  
  dR_n_p <- a_n_lo - (u_n * R_n_p * H_p) / (R_n_p + k_n);
  dR_p_p <- a_p_hi - (u_p * R_p_p * H_p) / (R_p_p + k_p);
  dQ_n_p <- (u_n * R_n_p) / (R_n_p + k_n) - min((1 - Qmin_n / Q_n_p), (1 - Qmin_p / Q_p_p)) * g * Q_n_p - min((1 - q_n / Q_n_p), (1 - q_p / Q_p_p)) * z_n * r * V_p
  dQ_p_p <- (u_p * R_p_p) / (R_p_p + k_p) - min((1 - Qmin_n / Q_n_p), (1 - Qmin_p / Q_p_p)) * g * Q_p_p - min((1 - q_n / Q_n_p), (1 - q_p / Q_p_p)) * z_p * r * V_p
  dH_p <- min((1 - Qmin_n / Q_n_p), (1 - Qmin_p / Q_p_p)) * g * H_p - m * H_p
  dV_p <- min((1 - q_n / Q_n_p), (1 - q_p / Q_p_p)) * r * V_p - c * V_p
  
  dR_n_np <- a_n_hi - (u_n * R_n_np * H_np) / (R_n_np + k_n);
  dR_p_np <- a_p_hi - (u_p * R_p_np * H_np) / (R_p_np + k_p);
  dQ_n_np <- (u_n * R_n_np) / (R_n_np + k_n) - min((1 - Qmin_n / Q_n_np), (1 - Qmin_p / Q_p_np)) * g * Q_n_np - min((1 - q_n / Q_n_np), (1 - q_p / Q_p_np)) * z_n * r * V_np
  dQ_p_np <- (u_p * R_p_np) / (R_p_np + k_p) - min((1 - Qmin_n / Q_n_np), (1 - Qmin_p / Q_p_np)) * g * Q_p_np - min((1 - q_n / Q_n_np), (1 - q_p / Q_p_np)) * z_p * r * V_np
  dH_np <- min((1 - Qmin_n / Q_n_np), (1 - Qmin_p / Q_p_np)) * g * H_np - m * H_np
  dV_np <- min((1 - q_n / Q_n_np), (1 - q_p / Q_p_np)) * r * V_np - c * V_np
  
  return(list(c(dR_n_low, dR_p_low, dQ_n_low, dQ_p_low, dH_low, dV_low,
                dR_n_n, dR_p_n, dQ_n_n, dQ_p_n, dH_n, dV_n,
                dR_n_p, dR_p_p, dQ_n_p, dQ_p_p, dH_p, dV_p,
                dR_n_np, dR_p_np, dQ_n_np, dQ_p_np, dH_np, dV_np)))
  })
}

virus1_model_init <- function(){
  
  # simulate plant growth
  plant_init <- ode(init_def1, times = seq(0, plant_days, length.out = 100), 
                        plant_model, params_def1) %>%
    as_tibble() %>%
    mutate(across(everything(), as.double)) %>%
    filter(time == max(time))
  
  # extract values
  y0 <- c(R_n_low = pull(plant_init, R_n_low), 
          R_p_low = pull(plant_init, R_p_low), 
          Q_n_low = pull(plant_init, Q_n_low), 
          Q_p_low = pull(plant_init, Q_p_low), 
          H_low = pull(plant_init, H_low), 
          V_low = V0,
          R_n_n = pull(plant_init, R_n_n), 
          R_p_n = pull(plant_init, R_p_n), 
          Q_n_n = pull(plant_init, Q_n_n), 
          Q_p_n = pull(plant_init, Q_p_n), 
          H_n = pull(plant_init, H_n), 
          V_n = V0,
          R_n_p = pull(plant_init, R_n_p), 
          R_p_p = pull(plant_init, R_p_p), 
          Q_n_p = pull(plant_init, Q_n_p), 
          Q_p_p = pull(plant_init, Q_p_p), 
          H_p = pull(plant_init, H_p), 
          V_p = V0,
          R_n_np = pull(plant_init, R_n_np), 
          R_p_np = pull(plant_init, R_p_np), 
          Q_n_np = pull(plant_init, Q_n_np), 
          Q_p_np = pull(plant_init, Q_p_np), 
          H_np = pull(plant_init, H_np), 
          V_np = V0)
  
  return(y0)
}

virus1_model_format <- function(mod_in){
  
  mod_out <- mod_in  %>%
    as_tibble() %>%
    mutate(across(everything(), as.double)) %>%
    pivot_longer(cols = -time,
                 names_to = "variable",
                 values_to = "value") %>%
    mutate(nutrient = case_when(str_ends(variable, "low") == T ~ "low",
                                str_ends(variable, "np") == T ~ "+N+P",
                                str_ends(variable, "n") == T ~ "+N",
                                str_ends(variable, "p") == T ~ "+P") %>%
             fct_relevel("low", "+N", "+P"),
           variable2 = case_when(str_starts(variable, "R_n") == T ~ "R_n",
                                 str_starts(variable, "R_p") == T ~ "R_p",
                                 str_starts(variable, "Q_n") == T ~ "Q_n",
                                 str_starts(variable, "Q_p") == T ~ "Q_p",
                                 str_starts(variable, "H") == T ~ "H",
                                 str_starts(variable, "V") == T ~ "V"))
  
  return(mod_out)
  
}


#### two virus model ####

virus2_model = function (t, yy, parms) {
  with(as.list(c(yy, parms)), {
    
    dR_n_low <- a_n_lo - (u_n * R_n_low * H_low) / (R_n_low + k_n);
    dR_p_low <- a_p_lo - (u_p * R_p_low * H_low) / (R_p_low + k_p);
    dQ_n_low <- (u_n * R_n_low) / (R_n_low + k_n) - min((1 - Qmin_n / Q_n_low), (1 - Qmin_p / Q_p_low)) * g * Q_n_low - min((1 - q_nb / Q_n_low), (1 - q_pb / Q_p_low)) * z_nb * r_b * V_b_low - min((1 - q_nc / Q_n_low), (1 - q_pc / Q_p_low)) * z_nc * r_c * V_c_low
    dQ_p_low <- (u_p * R_p_low) / (R_p_low + k_p) - min((1 - Qmin_n / Q_n_low), (1 - Qmin_p / Q_p_low)) * g * Q_p_low - min((1 - q_nb / Q_n_low), (1 - q_pb / Q_p_low)) * z_pb * r_b * V_b_low - min((1 - q_nc / Q_n_low), (1 - q_pc / Q_p_low)) * z_pc * r_c * V_c_low
    dH_low <- min((1 - Qmin_n / Q_n_low), (1 - Qmin_p / Q_p_low)) * g * H_low - m * H_low
    dV_b_low <- min((1 - q_nb / Q_n_low), (1 - q_pb / Q_p_low)) * r_b * V_b_low - c_b * V_b_low
    dV_c_low <- min((1 - q_nc / Q_n_low), (1 - q_pc / Q_p_low)) * r_c * V_c_low - c_c * V_c_low
    
    dR_n_n <- a_n_hi - (u_n * R_n_n * H_n) / (R_n_n + k_n);
    dR_p_n <- a_p_lo - (u_p * R_p_n * H_n) / (R_p_n + k_p);
    dQ_n_n <- (u_n * R_n_n) / (R_n_n + k_n) - min((1 - Qmin_n / Q_n_n), (1 - Qmin_p / Q_p_n)) * g * Q_n_n - min((1 - q_nb / Q_n_n), (1 - q_pb / Q_p_n)) * z_nb * r_b * V_b_n - min((1 - q_nc / Q_n_n), (1 - q_pc / Q_p_n)) * z_nc * r_c * V_c_n
    dQ_p_n <- (u_p * R_p_n) / (R_p_n + k_p) - min((1 - Qmin_n / Q_n_n), (1 - Qmin_p / Q_p_n)) * g * Q_p_n - min((1 - q_nb / Q_n_n), (1 - q_pb / Q_p_n)) * z_pb * r_b * V_b_n - min((1 - q_nc / Q_n_n), (1 - q_pc / Q_p_n)) * z_pc * r_c * V_c_n
    dH_n <- min((1 - Qmin_n / Q_n_n), (1 - Qmin_p / Q_p_n)) * g * H_n - m * H_n
    dV_b_n <- min((1 - q_nb / Q_n_n), (1 - q_pb / Q_p_n)) * r_b * V_b_n - c_b * V_b_n
    dV_c_n <- min((1 - q_nc / Q_n_n), (1 - q_pc / Q_p_n)) * r_c * V_c_n - c_c * V_c_n
    
    dR_n_p <- a_n_lo - (u_n * R_n_p * H_p) / (R_n_p + k_n);
    dR_p_p <- a_p_hi - (u_p * R_p_p * H_p) / (R_p_p + k_p);
    dQ_n_p <- (u_n * R_n_p) / (R_n_p + k_n) - min((1 - Qmin_n / Q_n_p), (1 - Qmin_p / Q_p_p)) * g * Q_n_p - min((1 - q_nb / Q_n_p), (1 - q_pb / Q_p_p)) * z_nb * r_b * V_b_p - min((1 - q_nc / Q_n_p), (1 - q_pc / Q_p_p)) * z_nc * r_c * V_c_p
    dQ_p_p <- (u_p * R_p_p) / (R_p_p + k_p) - min((1 - Qmin_n / Q_n_p), (1 - Qmin_p / Q_p_p)) * g * Q_p_p - min((1 - q_nb / Q_n_p), (1 - q_pb / Q_p_p)) * z_pb * r_b * V_b_p - min((1 - q_nc / Q_n_p), (1 - q_pc / Q_p_p)) * z_pc * r_c * V_c_p
    dH_p <- min((1 - Qmin_n / Q_n_p), (1 - Qmin_p / Q_p_p)) * g * H_p - m * H_p
    dV_b_p <- min((1 - q_nb / Q_n_p), (1 - q_pb / Q_p_p)) * r_b * V_b_p - c_b * V_b_p
    dV_c_p <- min((1 - q_nc / Q_n_p), (1 - q_pc / Q_p_p)) * r_c * V_c_p - c_c * V_c_p
    
    dR_n_np <- a_n_hi - (u_n * R_n_np * H_np) / (R_n_np + k_n);
    dR_p_np <- a_p_hi - (u_p * R_p_np * H_np) / (R_p_np + k_p);
    dQ_n_np <- (u_n * R_n_np) / (R_n_np + k_n) - min((1 - Qmin_n / Q_n_np), (1 - Qmin_p / Q_p_np)) * g * Q_n_np - min((1 - q_nb / Q_n_np), (1 - q_pb / Q_p_np)) * z_nb * r_b * V_b_np - min((1 - q_nc / Q_n_np), (1 - q_pc / Q_p_np)) * z_nc * r_c * V_c_np
    dQ_p_np <- (u_p * R_p_np) / (R_p_np + k_p) - min((1 - Qmin_n / Q_n_np), (1 - Qmin_p / Q_p_np)) * g * Q_p_np - min((1 - q_nb / Q_n_np), (1 - q_pb / Q_p_np)) * z_pb * r_b * V_b_np - min((1 - q_nc / Q_n_np), (1 - q_pc / Q_p_np)) * z_pc * r_c * V_c_np
    dH_np <- min((1 - Qmin_n / Q_n_np), (1 - Qmin_p / Q_p_np)) * g * H_np - m * H_np
    dV_b_np <- min((1 - q_nb / Q_n_np), (1 - q_pb / Q_p_np)) * r_b * V_b_np - c_b * V_b_np
    dV_c_np <- min((1 - q_nc / Q_n_np), (1 - q_pc / Q_p_np)) * r_c * V_c_np - c_c * V_c_np
    
    return(list(c(dR_n_low, dR_p_low, dQ_n_low, dQ_p_low, dH_low, dV_b_low, dV_c_low,
                  dR_n_n, dR_p_n, dQ_n_n, dQ_p_n, dH_n, dV_b_n, dV_c_n,
                  dR_n_p, dR_p_p, dQ_n_p, dQ_p_p, dH_p, dV_b_p, dV_c_p,
                  dR_n_np, dR_p_np, dQ_n_np, dQ_p_np, dH_np, dV_b_np, dV_c_np)))
    
  })
}

virus2_model_sim <- function(params, first_virus, V0_b, V0_c, plant_time, res_time, inv_time, inits = init_def2){
  
  # simulate plant growth
  plant_mod <- ode(inits, times = seq(0, plant_time, length.out = 100),
                   virus2_model, params) %>%
    as_tibble() %>%
    mutate(across(everything(), as.double)) 
  
  plant_init <- plant_mod %>%
    filter(time == plant_time)
  
  y0_first_virus <- c(R_n_low = pull(plant_init, R_n_low), 
                      R_p_low = pull(plant_init, R_p_low), 
                      Q_n_low = pull(plant_init, Q_n_low), 
                      Q_p_low = pull(plant_init, Q_p_low), 
                      H_low = pull(plant_init, H_low), 
                      V_b_low = if_else(first_virus == "PAV", V0_b, 0),
                      V_c_low = if_else(first_virus == "RPV", V0_c, 0), 
                      R_n_n = pull(plant_init, R_n_n), 
                      R_p_n = pull(plant_init, R_p_n), 
                      Q_n_n = pull(plant_init, Q_n_n), 
                      Q_p_n = pull(plant_init, Q_p_n), 
                      H_n = pull(plant_init, H_n), 
                      V_b_n = if_else(first_virus == "PAV", V0_b, 0),
                      V_c_n = if_else(first_virus == "RPV", V0_c, 0), 
                      R_n_p = pull(plant_init, R_n_p), 
                      R_p_p = pull(plant_init, R_p_p), 
                      Q_n_p = pull(plant_init, Q_n_p), 
                      Q_p_p = pull(plant_init, Q_p_p), 
                      H_p = pull(plant_init, H_p), 
                      V_b_p = if_else(first_virus == "PAV", V0_b, 0),
                      V_c_p = if_else(first_virus == "RPV", V0_c, 0), 
                      R_n_np = pull(plant_init, R_n_np), 
                      R_p_np = pull(plant_init, R_p_np), 
                      Q_n_np = pull(plant_init, Q_n_np), 
                      Q_p_np = pull(plant_init, Q_p_np), 
                      H_np = pull(plant_init, H_np), 
                      V_b_np = if_else(first_virus == "PAV", V0_b, 0),
                      V_c_np = if_else(first_virus == "RPV", V0_c, 0))
  
  # first virus model
  first_virus_mod <- ode(y0_first_virus, seq(0, res_time, length.out = 100),
                         virus2_model, params) %>%
    as_tibble() %>%
    mutate(across(everything(), as.double)) 
  
  first_virus_init <- first_virus_mod %>%
    filter(time == res_time)
  
  # virus initial conditions
  y0_second_virus <- c(R_n_low = pull(first_virus_init, R_n_low), 
                       R_p_low = pull(first_virus_init, R_p_low), 
                       Q_n_low = pull(first_virus_init, Q_n_low), 
                       Q_p_low = pull(first_virus_init, Q_p_low), 
                       H_low = pull(first_virus_init, H_low), 
                       V_b_low = if_else(first_virus == "PAV", 
                                         pull(first_virus_init, V_b_low), V0_b),
                       V_c_low = if_else(first_virus == "RPV", 
                                         pull(first_virus_init, V_c_low), V0_c), 
                       R_n_n = pull(first_virus_init, R_n_n), 
                       R_p_n = pull(first_virus_init, R_p_n), 
                       Q_n_n = pull(first_virus_init, Q_n_n), 
                       Q_p_n = pull(first_virus_init, Q_p_n), 
                       H_n = pull(first_virus_init, H_n), 
                       V_b_n = if_else(first_virus == "PAV", 
                                       pull(first_virus_init, V_b_n), V0_b),
                       V_c_n = if_else(first_virus == "RPV", 
                                       pull(first_virus_init, V_c_n), V0_c),
                       R_n_p = pull(first_virus_init, R_n_p), 
                       R_p_p = pull(first_virus_init, R_p_p), 
                       Q_n_p = pull(first_virus_init, Q_n_p), 
                       Q_p_p = pull(first_virus_init, Q_p_p), 
                       H_p = pull(first_virus_init, H_p), 
                       V_b_p = if_else(first_virus == "PAV", 
                                       pull(first_virus_init, V_b_p), V0_b),
                       V_c_p = if_else(first_virus == "RPV", 
                                       pull(first_virus_init, V_c_p), V0_c),
                       R_n_np = pull(first_virus_init, R_n_np), 
                       R_p_np = pull(first_virus_init, R_p_np), 
                       Q_n_np = pull(first_virus_init, Q_n_np), 
                       Q_p_np = pull(first_virus_init, Q_p_np), 
                       H_np = pull(first_virus_init, H_np), 
                       V_b_np = if_else(first_virus == "PAV", 
                                        pull(first_virus_init, V_b_np), V0_b),
                       V_c_np = if_else(first_virus == "RPV", 
                                        pull(first_virus_init, V_c_np), V0_c))
  
  # first virus model
  second_virus_mod <- ode(y0_second_virus, seq(0, inv_time, length.out = 100),
                          virus2_model, params) %>%
    as_tibble() %>%
    mutate(across(everything(), as.double))
  
  # combine simulations
  comb_mod <- plant_mod %>%
    filter(time != max(time)) %>%
    full_join(first_virus_mod %>%
                filter(time != max(time)) %>%
                mutate(time = time + plant_time)) %>%
    full_join(second_virus_mod %>%
                mutate(time = time + plant_time + res_time))
  
  return(comb_mod)
  
}

virus2_model_format <- function(mod_in, params){
  
  # min nutrient conc.
  Qmin_n <- as.numeric(params["Qmin_n"])
  Qmin_p <- as.numeric(params["Qmin_p"])
  
  # format
  mod_out <- mod_in  %>%
    as_tibble() %>%
    mutate(across(everything(), as.double)) %>%
    mutate(VbH_low = V_b_low * H_low,
           VbH_n = V_b_n * H_n,
           VbH_p = V_b_p * H_p,
           VbH_np = V_b_np * H_np,
           VcH_low = V_c_low * H_low,
           VcH_n = V_c_n * H_n,
           VcH_p = V_c_p * H_p,
           VcH_np = V_c_np * H_np,
           Qlim_n_low = Qmin_n / Q_n_low,
           Qlim_n_n = Qmin_n / Q_n_n,
           Qlim_n_p = Qmin_n / Q_n_p,
           Qlim_n_np = Qmin_n / Q_n_np,
           Qlim_p_low = Qmin_p / Q_p_low,
           Qlim_p_n = Qmin_p / Q_p_n,
           Qlim_p_p = Qmin_p / Q_p_p,
           Qlim_p_np = Qmin_p / Q_p_np) %>%
    rowwise() %>%
    mutate(LimQ_low = max(Qlim_n_low, Qlim_p_low),
           LimQ_n = max(Qlim_n_n, Qlim_p_n),
           LimQ_p = max(Qlim_n_p, Qlim_p_p),
           LimQ_np = max(Qlim_n_np, Qlim_p_np)) %>%
    ungroup() %>%
    mutate(lim_N_H_low = case_when(LimQ_low == Qlim_n_low ~ 1,
                                   TRUE ~ 0),
           lim_N_H_n = case_when(LimQ_n == Qlim_n_n ~ 1,
                                 TRUE ~ 0),
           lim_N_H_p = case_when(LimQ_p == Qlim_n_p ~ 1,
                                 TRUE ~ 0),
           lim_N_H_np = case_when(LimQ_np == Qlim_n_np ~ 1,
                                  TRUE ~ 0)) %>%
    pivot_longer(cols = -time,
                 names_to = "variable",
                 values_to = "value") %>%
    mutate(nutrient = case_when(str_ends(variable, "low") == T ~ "low",
                                str_ends(variable, "np") == T ~ "+N+P",
                                str_ends(variable, "n") == T ~ "+N",
                                str_ends(variable, "p") == T ~ "+P") %>%
             fct_relevel("low", "+N", "+P"),
           variable2 = case_when(str_starts(variable, "R_n") == T ~ "R_n",
                                 str_starts(variable, "R_p") == T ~ "R_p",
                                 str_starts(variable, "Q_n") == T ~ "Q_n",
                                 str_starts(variable, "Q_p") == T ~ "Q_p",
                                 str_starts(variable, "H") == T ~ "H",
                                 str_starts(variable, "V_b") == T ~ "PAV_conc",
                                 str_starts(variable, "V_c") == T ~ "RPV_conc",
                                 str_starts(variable, "VbH") == T ~ "PAV_pop",
                                 str_starts(variable, "VcH") == T ~ "RPV_pop",
                                 str_starts(variable, "Qlim_n") == T ~ "Qlim_n",
                                 str_starts(variable, "Qlim_p") == T ~ "Qlim_p",
                                 str_starts(variable, "LimQ") == T ~ "Qlim",
                                 str_starts(variable, "lim_N") == T ~ "lim_N"))
  
  # add limiting nutrient to all rows
  mod_out2 <- mod_out %>%
    filter(!(variable2 %in% c("Qlim_n", "Qlim_p", "LimQ", "lim_N"))) %>%
    full_join(mod_out %>%
                filter(variable2 == "lim_N") %>%
                mutate(value = fct_recode(as.character(value), 
                                          "Q[N]" = "1", "Q[P]" = "0") %>%
                         fct_relevel("Q[N]")) %>%
                rename("lim_nut_H" = "value") %>%
                select(nutrient, time, lim_nut_H))
  
  return(mod_out2)
  
}

plant_fig_fun <- function(mod_dat, params, q_adj_n, q_adj_p){
  
  # min nutrient conc.
  Qmin_n <- as.numeric(params["Qmin_n"])
  Qmin_p <- as.numeric(params["Qmin_p"])
  
  # figure labels
  Qn_lab <- tibble(time = 0, 
                   value = Qmin_n, 
                   label = "Q['min,N']")
  Qp_lab <- tibble(time = 0, 
                   value = Qmin_p, 
                   label = "Q['min,P']")
  
  # panels
  plant_H_fig <- filter(mod_dat, variable2 == "H") %>%
    ggplot(aes(time, value, linetype = nutrient, color = nutrient)) +
    geom_line(aes(size = lim_nut_H)) +
    scale_color_viridis_d(end = 0.7, name = "Nutrient\nsupply") +
    scale_linetype(name = "Nutrient\nsupply") +
    scale_size_manual(values = c(1.2, 0.6), guide= "none") +
    labs(x = "Time (days)", y = "Plant biomass (g)", title = "(A)") +
    fig_theme +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank()) +
    guides(color = guide_legend(override.aes = list(size = 1.2)))
  
  plant_lim_fig <- filter(mod_dat, variable2 == "Qlim") %>%
    ggplot(aes(time, value, linetype = nutrient, color = nutrient)) +
    geom_line(aes(size = lim_nut_H)) +
    scale_color_viridis_d(end = 0.7, guide = "none") +
    scale_linetype(guide = "none") +
    scale_size_manual(values = c(1.2, 0.6), name = "Limiting\nnutrient",
                      labels = c("N", "P")) +
    labs(x = "Time (days)", title = "(B)",
         y = expression(paste("Limiting nutrient ratio (", Q[min], "/Q)", sep = ""))) +
    fig_theme +
    theme(axis.title.x = element_blank(),
          axis.text.x = element_blank())
  
  plant_Rn_fig <- filter(mod_dat, variable2 == "R_n") %>%
    ggplot(aes(time, value, linetype = nutrient, color = nutrient)) +
    geom_line(aes(size = lim_nut_H)) +
    scale_color_viridis_d(end = 0.7) +
    scale_size_manual(values = c(1.2, 0.6)) +
    labs(x = "Time (days)", y = "Environment N (g)", title = "(C)") +
    fig_theme +
    theme(legend.position = "none",
          axis.title.x = element_blank(),
          axis.text.x = element_blank())
  
  plant_Rp_fig <- filter(mod_dat, variable2 == "R_p") %>%
    ggplot(aes(time, value, linetype = nutrient, color = nutrient)) +
    geom_line(aes(size = lim_nut_H)) +
    scale_color_viridis_d(end = 0.7) +
    scale_size_manual(values = c(1.2, 0.6)) +
    labs(x = "", y = "Environment P (g)", title = "(D)") +
    fig_theme +
    theme(legend.position = "none")
  
  plant_Qn_fig <- filter(mod_dat, variable2 == "Q_n") %>%
    ggplot(aes(time, value)) +
    geom_hline(yintercept = Qmin_n, color = "black") +
    geom_text(data = Qn_lab, aes(label = label), parse = T, color = "black", fontface = "italic",
              size = 3, hjust = 0, vjust = 0, nudge_y = q_adj_n) +
    geom_line(aes(linetype = nutrient, color = nutrient, size = lim_nut_H)) +
    scale_color_viridis_d(end = 0.7) +
    scale_size_manual(values = c(1.2, 0.6)) +
    labs(x = "Time (days)", y = "Plant N concentration", title = "(E)") +
    fig_theme +
    theme(legend.position = "none")
  
  plant_Qp_fig <- filter(mod_dat, variable2 == "Q_p") %>%
    ggplot(aes(time, value)) +
    geom_hline(yintercept = Qmin_p, color = "black") +
    geom_text(data = Qp_lab, aes(label = label), parse = T, color = "black", fontface = "italic",
              size = 3, hjust = 0, vjust = 0, nudge_y = q_adj_p) +
    geom_line(aes(linetype = nutrient, color = nutrient, size = lim_nut_H)) +
    scale_color_viridis_d(end = 0.7) +
    scale_size_manual(values = c(1.2, 0.6)) +
    labs(x = "", y = "Plant P concentration", title = "(F)") +
    fig_theme +
    theme(legend.position = "none")
  
  # extract legends
  nut_leg <- get_legend(plant_H_fig)
  lim_leg <- get_legend(plant_lim_fig)
  
  # combine figures
  plant_top_fig <- plot_grid(plant_H_fig + theme(legend.position = "none"),
                             plant_lim_fig + theme(legend.position = "none"),
                             plant_Rn_fig,
                             nut_leg,
                             nrow = 1,
                             rel_widths = c(1, 1, 1, 0.35))
  plant_bot_fig <- plot_grid(plant_Rp_fig, plant_Qn_fig, plant_Qp_fig,
                             lim_leg,
                             nrow = 1,
                             rel_widths = c(1, 1, 1, 0.35))
  
  # output
  return(plot_grid(plant_top_fig, plant_bot_fig,
                   nrow = 2, rel_heights = c(0.85, 1)))
}


#### figure settings ####

fig_theme <- theme_bw() +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.spacing.x = unit(0,"line"),
        axis.text.y = element_text(size = 8, color = "black"),
        axis.text.x = element_text(size = 8, color = "black"),
        axis.title = element_text(size = 8),
        axis.line = element_line(color = "black"),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 8),
        legend.background = element_blank(),
        legend.key.width = unit(7, "mm"),
        legend.key.heigh = unit(5, "mm"),
        strip.background = element_blank(),
        strip.text = element_text(size = 8, hjust = 0.5),
        strip.placement = "outside",
        plot.title = element_text(size = 8, vjust = 0, face = "bold"))
