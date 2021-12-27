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
                 m = 0.006,
                 g = 0.144,
                 c_b = 0.002,
                 c_c = 0.005,
                 r_b = 0.09,
                 r_c = 0.27)


#### initial values ####

# initial values
H0 <- 1e-2
Q0_n <- as.numeric(params_def1["Qmin_n"]) * 10
Q0_p <- as.numeric(params_def1["Qmin_p"]) * 10
R0_n_lo <- as.numeric(params_def1["a_n_hi"]) * 7 + as.numeric(params_def1["a_n_lo"]) * 3
R0_n_hi <- as.numeric(params_def1["a_n_hi"]) * 10
R0_p_lo <- as.numeric(params_def1["a_p_hi"]) * 7 + as.numeric(params_def1["a_p_lo"]) * 3
R0_p_hi <- as.numeric(params_def1["a_p_hi"]) * 10
V0_init <- 5e5

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
                                 str_starts(variable, "H") == T ~ "H"))
  
  return(mod_out)
  
}


#### start here: plant and single virus model ####


#### model ####

plant_virus_model = function (t, yy, parms) { 
  
  # supply rates
  a_n = parms[1];
  a_p = parms[2];
  
  # set initial values
  R_n = yy[1];
  R_p = yy[2];
  Q_n = yy[3];
  Q_p = yy[4];
  H = yy[5];
  V_b = yy[6];
  V_c = yy[7];
  
  # model
  dR_n = a_n - (u_n * R_n * H) / (R_n + k_n);
  dR_p = a_p - (u_p * R_p * H) / (R_p + k_p);
  dQ_n = (u_n * R_n) / (R_n + k_n) - min((1 - Qmin_n / Q_n), (1 - Qmin_p / Q_p)) * g * Q_n - min((1 - q_nb / Q_n), (1 - q_pb / Q_p)) * z_nb * r_b * V_b - min((1 - q_nc / Q_n), (1 - q_pc / Q_p)) * z_nc * r_c * V_c
  dQ_p = (u_p * R_p) / (R_p + k_p) - min((1 - Qmin_n / Q_n), (1 - Qmin_p / Q_p)) * g * Q_p - min((1 - q_nb / Q_n), (1 - q_pb / Q_p)) * z_pb * r_b * V_b - min((1 - q_nc / Q_n), (1 - q_pc / Q_p)) * z_pc * r_c * V_c
  dH = min((1 - Qmin_n / Q_n), (1 - Qmin_p / Q_p)) * g * H - m * H
  dV_b = min((1 - q_nb / Q_n), (1 - q_pb / Q_p)) * r_b * V_b - c_b * V_b
  dV_c = min((1 - q_nc / Q_n), (1 - q_pc / Q_p)) * r_c * V_c - c_c * V_c
  
  return(list(c(dR_n, dR_p, dQ_n, dQ_p, dH, dV_b, dV_c)))
}


#### simulation wrapper ####

sim_fun <- function(N_level, P_level, nut_trt, first_virus, plant_time, res_time, inv_time){
  
  # initial resource values
  if(N_level == "low"){
    R0_n <- R0_n_lo
    a_n <- a_n_lo
  }else{
    R0_n <- R0_n_hi
    a_n <- a_n_hi
  }
  
  if(P_level == "low"){
    R0_p <- R0_p_lo
    a_p <- a_p_lo
  }else{
    R0_p <- R0_p_hi
    a_p <- a_p_hi
  }
  
  # plant model
  plant_mod <- ode(c(R_n = R0_n, R_p = R0_p, Q_n = Q0_n, Q_p = Q0_p, H = H0, V_b = 0, V_c = 0),
                   seq(0, plant_time, length.out = 100), plant_virus_model, c(a_n = a_n, a_p = a_p)) %>%
    as_tibble() %>%
    mutate(across(everything(), as.double))
  
  # final time
  plant_init <- plant_mod %>%
    filter(time == plant_time)
  
  # virus initial conditions
  y0_first_virus <- c(R_n = pull(plant_init, R_n), 
                      R_p = pull(plant_init, R_p), 
                      Q_n = pull(plant_init, Q_n), 
                      Q_p = pull(plant_init, Q_p), 
                      H = pull(plant_init, H), 
                      V_b = if_else(first_virus == "PAV", V0_b, 0),
                      V_c = if_else(first_virus == "RPV", V0_c, 0))
  
  # first virus model
  first_virus_mod <- ode(y0_first_virus, seq(0, res_time, length.out = 100),
                         plant_virus_model, c(a_n = a_n, a_p = a_p)) %>%
    as_tibble() %>%
    mutate(across(everything(), as.double))
  
  # final time
  first_virus_init <- first_virus_mod %>%
    filter(time == res_time)
  
  # virus initial conditions
  y0_second_virus <- c(R_n = pull(first_virus_init, R_n),
                       R_p = pull(first_virus_init, R_p),
                       Q_n = pull(first_virus_init, Q_n),
                       Q_p = pull(first_virus_init, Q_p),
                       H = pull(first_virus_init, H),
                       V_b = if_else(first_virus == "PAV", pull(first_virus_init, V_b), V0_b),
                       V_c = if_else(first_virus == "RPV", pull(first_virus_init, V_c), V0_c))
  
  # first virus model
  second_virus_mod <- ode(y0_second_virus, seq(0, inv_time, length.out = 100),
                          plant_virus_model, c(a_n = a_n, a_p = a_p)) %>%
    as_tibble() %>%
    mutate(across(everything(), as.double))
  
  # combine models
  mod_out <- plant_mod %>%
    filter(time != max(time)) %>%
    full_join(first_virus_mod %>%
                filter(time != max(time)) %>%
                mutate(time = time + plant_time)) %>%
    full_join(second_virus_mod %>%
                mutate(time = time + plant_time + res_time)) %>%
    mutate(nutrient = nut_trt,
           resident = first_virus,
           invader = if_else(resident == "PAV", "RPV", "PAV"))
  
}


#### model with more flexible parameters ####

plant_virus_param_model = function (t, yy, parms) {
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
