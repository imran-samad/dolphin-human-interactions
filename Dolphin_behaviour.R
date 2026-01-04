############################################
## Libraries
############################################
library(tidyverse)
library(mclogit)
library(caret)
library(pscl)
library(performance)

############################################
## Load data
############################################
tracks_filtered = read.csv("Data/dolphin_tracks_processed.csv")

############################################
## DATA PREPARATION
############################################

## ---- For modelling STATE TRANSITIONS (trans == 1)
tracks_trans = tracks_filtered |> 
  filter(trans == 1) |> 
  mutate(cpue = scale(cpue))

## ---- For modelling STATE PROBABILITIES (trans == 0)
tracks_state = tracks_filtered |> 
  filter(trans == 0) |> 
  mutate(cpue = scale(cpue))

## ---- Segregate treatments (state probabilities) # do the same for TRANSITION probabilities
tracks_state_1 = tracks_state |> filter(treatment == "no_boats")
tracks_state_2 = tracks_state |> filter(treatment == "TB")
tracks_state_3 = tracks_state |> filter(treatment == "FB")
tracks_state_4 = tracks_state |> filter(treatment == "both")

############################################
## MULTINOMIAL LOGISTIC MODELS (STATE PROBABILITY)
############################################

## ---- No boats
mm1 = mblogit(
  behaviour_f ~ (pod_size + season) + neighbour,
  data = tracks_state_1,
  estimator = "ML"
)
summary(mm1)

## ---- Tourist boats
mm2 = mblogit(
  behaviour_f ~ (pod_size + closest_boat_dist + n_boats_200m + season) + neighbour,
  data = tracks_state_2,
  estimator = "ML"
)
summary(mm2)

## ---- Fishing boats
mm3 = mblogit(
  behaviour_f ~ (pod_size + closest_net_dist + season) + neighbour,
  data = tracks_state_3,
  estimator = "ML"
)
summary(mm3)

## ---- Both tourist + fishing boats
mm4 = mblogit(
  behaviour_f ~ (pod_size + closest_boat_dist + n_boats_200m +
                   closest_net_dist + season) + neighbour,
  data = tracks_state_4,
  estimator = "ML"
)
summary(mm4)

############################################
## MODEL DIAGNOSTICS (example: mm1)
############################################

## ---- Predicted vs observed
pred = predict(mm1, type = "response")
pred_class = colnames(pred)[apply(pred, 1, which.max)]

used_rows = as.numeric(rownames(pred))
obs_used  = tracks_state_1$behaviour_f[used_rows]

pred_f = factor(pred_class, levels = levels(obs_used))
obs_f  = factor(obs_used,  levels = levels(obs_used))

cm = confusionMatrix(pred_f, obs_f)
cm$overall["Kappa"]
mean(pred_f == obs_f)
cm$byClass

## ---- Pseudo-R2 and residual checks
pR2(mm1)
check_model(mm1)

############################################
## NULL MODELS FOR MODEL COMPARISON
############################################

mm10 = mblogit(behaviour_f ~ 1, data = tracks_state_1, estimator = "ML")

vars_used = c("behaviour_f","pod_size","closest_boat_dist",
              "n_boats_200m","season","neighbour")
mm20 = mblogit(
  behaviour_f ~ 1,
  data = tracks_state_2[complete.cases(tracks_state_2[, vars_used]),],
  estimator = "ML"
)

vars_used = c("behaviour_f","pod_size","closest_net_dist","season","neighbour")
mm30 = mblogit(
  behaviour_f ~ 1,
  data = tracks_state_3[complete.cases(tracks_state_3[, vars_used]),],
  estimator = "ML"
)

vars_used = c("behaviour_f","pod_size","closest_net_dist",
              "closest_boat_dist","n_boats_200m","season","neighbour")
mm40 = mblogit(
  behaviour_f ~ 1,
  data = tracks_state_4[complete.cases(tracks_state_4[, vars_used]),],
  estimator = "ML"
)

anova(mm10, mm1, test = "Chisq")
anova(mm20, mm2, test = "Chisq")
anova(mm30, mm3, test = "Chisq")
anova(mm40, mm4, test = "Chisq")

## ---- McFadden's R2 (explicit)
1 - (logLik(mm1) / logLik(mm10))

############################################
## EXTRACTING COEFFICIENTS
############################################

extract_parms = function(model, model_id) {
  summ = summary(model)
  df = data.frame(summ$coefficients)
  df$response = sapply(strsplit(rownames(df), "~"), `[`, 1)
  df$cov      = sapply(strsplit(rownames(df), "~"), `[`, 2)
  rownames(df) = NULL
  df$model = model_id
  df
}

parms = bind_rows(
  extract_parms(mm1, 1),
  extract_parms(mm2, 2),
  extract_parms(mm3, 3),
  extract_parms(mm4, 4)
) |> 
  arrange(model, response)

parms[,1:4] = sapply(parms[,1:4], round, 3)

parms = parms |>
  mutate(
    cov = case_when(
      cov == "pod_size" ~ "Pod size",
      cov == "seasonpre" ~ "Season: low fish catch",
      cov == "neighbour" ~ "Neighbour presence",
      cov == "closest_boat_dist" ~ "Distance to closest boat",
      cov == "n_boats_200m" ~ "Number of surrounding boats",
      cov == "closest_net_dist" ~ "Distance to closest fishing net",
      TRUE ~ cov
    ),
    cov = factor(
      cov,
      levels = c(
        "Pod size",
        "Neighbour presence",
        "Distance to closest boat",
        "Number of surrounding boats",
        "Distance to closest fishing net",
        "Season: low fish catch",
        "(Intercept)"
      )
    )
  )

############################################
## STATE PROBABILITIES (SEQUENCE-STANDARDISED)
############################################

calc_state_prob = function(dat) {
  dat |>
    select(seq_grp, behaviour) |>
    na.omit() |>
    group_by(seq_grp) |>
    count(behaviour) |>
    complete(
      behaviour = c("forage","socialise","travel","escape"),
      fill = list(n = 0)
    ) |>
    mutate(b_prob = n / sum(n)) |>
    ungroup() |>
    group_by(behaviour) |>
    reframe(
      mean = mean(b_prob),
      se   = sd(b_prob) / sqrt(n())
    ) |>
    mutate(tr = 4)
}

b1 = calc_state_prob(tracks_state_1)
b2 = calc_state_prob(tracks_state_2)
b3 = calc_state_prob(tracks_state_3)
b4 = calc_state_prob(tracks_state_4)

base_prob = bind_rows(b1, b2, b3, b4)

############################################
## TRANSITION PROBABILITIES
############################################

calc_trans_prob = function(dat) {
  dat |>
    select(seq_grp, behaviour_lag, behaviour) |>
    na.omit() |>
    group_by(seq_grp) |>
    count(behaviour_lag, behaviour) |>
    ungroup() |>
    group_by(seq_grp) |>
    complete(
      behaviour_lag = c("forage","socialise","travel","escape"),
      behaviour     = c("forage","socialise","travel","escape"),
      fill = list(n = 0)
    ) |>
    ungroup() |>
    group_by(behaviour_lag, behaviour) |>
    summarise(n = sum(n), .groups = "drop") |>
    group_by(behaviour_lag) |>
    mutate(b_prob_mean = n / sum(n), tr = 1)
}

t1 = calc_trans_prob(tracks_state_1)
t2 = calc_trans_prob(tracks_state_2)
t3 = calc_trans_prob(tracks_state_3)
t4 = calc_trans_prob(tracks_state_4)

trans_prob = bind_rows(t1, t2, t3, t4)

############################################
## SPLITTING & MERGING ANALYSIS
############################################

tracks = read.csv("Data/dolphin_tracks.csv")
tracks_m_s = data.frame()

for (i in unique(tracks$seq_grp)) {
  
  seq_i = tracks$seq[tracks$seq_grp == i][1]
  
  ## ---- Splitting
  if (!is.na(tracks$split[tracks$seq_grp == i][1])) {
    tm = tracks$tm_split[tracks$seq_grp == i][1]
    grp = tracks$cluster_prev[tracks$seq_grp == i][1]
    
    row_s = which(tracks$tm_seq == tm & tracks$seq == seq_i & tracks$group == grp)
    row_b = which(tracks$tm_seq %in% (tm-10):tm &
                    tracks$seq == seq_i & tracks$group == grp)[1]
    
    tracks_m_s = rbind(
      tracks_m_s,
      cbind(tracks[c(row_s, row_b),],
            splitting = c(1, 0),
            merging   = c(NA, NA))
    )
  }
  
  ## ---- Merging
  if (!is.na(tracks$merge[tracks$seq_grp == i][1])) {
    tm = tracks$tm_merge[tracks$seq_grp == i][1]
    grp = tracks$group[tracks$seq_grp == i][1]
    
    row_m = which(tracks$tm_seq == tm & tracks$seq == seq_i & tracks$group == grp)
    row_b = which(tracks$tm_seq %in% (tm-10):tm &
                    tracks$seq == seq_i & tracks$group == grp)[1]
    
    tracks_m_s = rbind(
      tracks_m_s,
      cbind(tracks[c(row_m, row_b),],
            splitting = c(NA, NA),
            merging   = c(1, 0))
    )
  }
}

############################################
## Group SPLITTING / MERGING dynamics
############################################

tracks_m_s_t1 = tracks_m_s |>
  filter(treatment == "both") |> # change for different pseudo-treatments (no_boats, TB, FB, both)
  mutate(across(
    c(pod_cohesion, pod_size, closest_boat_dist, n_boats_200m,
      closest_dol_grp_dist, closest_dol_grp_size, closest_net_dist),
    ~ as.numeric(scale(.))
  ))

summary(glm(splitting ~ closest_net_dist,
            data = tracks_m_s_t1,
            family = "binomial"))

summary(glm(merging ~ n_boats_200m,
            data = tracks_m_s_t1,
            family = "binomial"))
#Assessment: the odds of neither merging, nor splitting are impacted by any of the covariates


############################################
## ADDITIONAL QUESTIONS
############################################

## ---- Between-group distance
summary(lm(
  closest_dol_grp_dist ~ n_boats_200m,
  data = tracks |> filter(treatment == "both") # change for different pseudo-treatments
))
#Assessment: groups moved closer when n_boats increased (as well as when boats approached closer)

## ---- Within-group cohesion
summary(lm(
  pod_cohesion ~ n_boats_200m,
  data = tracks |> filter(treatment == "TB") # change for different pseudo-treatments
))
#Assessment: dolphins within a group moved closer when n_boats increased (as well as when boats approached closer)