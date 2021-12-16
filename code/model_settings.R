#### parameters ####

# supply rates
a_n_lo <- 1.1e-6
a_n_hi <- 5.6e-5
a_p_lo <- 1.6e-7
a_p_hi <- 8.2e-6

# other parameters
u_n <- 2.9e-3
u_p <- 1.4e-3
k_n <- 4.9e-5
k_p <- 3.0e-5
Qmin_n <- 1.1e-3
Qmin_p <- 7.4e-5
q_n <- q_nb <- q_nc <- 1.1e-3
q_p <- q_pb <- q_pc <- 7.4e-5
z_nb <- 1.6e-18
z_pb <- 2.6e-19
z_nc <- 1.7e-18
z_pc <- 2.6e-19

m <- 0.007
g <- 0.144
c_b <- 0.002
c_c <- 0.005
r_b <- 0.09
r_c <- 0.27

# initial values
H0 <- 1e-2
Q0_const <- 10
Q0_n <- Qmin_n * Q0_const
Q0_p <- Qmin_p * Q0_const
R0_const <- 10
R0_n_lo <- a_n_hi * R0_const
R0_n_hi <- a_n_hi * R0_const
R0_p_lo <- a_p_hi * R0_const
R0_p_hi <- a_p_hi * R0_const
V0_init <- 5e5


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
    full_join(first_virus_mod %>%
                mutate(time = time + plant_time)) %>%
    full_join(second_virus_mod %>%
                mutate(time = time + plant_time + res_time)) %>%
    mutate(nutrient = nut_trt,
           resident = first_virus,
           invader = if_else(resident == "PAV", "RPV", "PAV"))
  
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
