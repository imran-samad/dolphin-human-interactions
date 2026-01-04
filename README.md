# Data and code for replicating the results from Samad, I., Patil, H., Cantor M., Farine D. R., Sutaria, D. & Shanker, K. (in revision). Human activities drive novel behaviours and transitions in dolphins. Marine Ecology Progress Series. https://doi.org/10.1101/2025.07.10.663819

## Project overview

This repository contains data processing pipelines and statistical models used to quantify how Indo-Pacific humpback dolphins respond behaviourally to tourism and fisheries activities. Using high-resolution drone-derived movement tracks, we model:

* **State probabilities** (likelihood of dolphins being in behavioural states such as forage, travel, socialise, escape)
* **State transition probabilities** (how dolphins switch between behaviours)
* **Group dynamics** (probability of group splitting and merging)

The analyses rely on multinomial logistic regression and binomial GLMs, stratified by treatment (no boats, tourist boats, fishing boats, both), while accounting for sequence structure and neighbourhood effects.
Note: Dolphin group data was extracted from drone videos using the package FastGeoRef https://github.com/imran-samad/fastgeoref

---

## Input data files

### 1. `Data/dolphin_tracks.csv`

Raw, second‑by‑second dolphin tracking data derived from drone footage. Each row represents a dolphin group at a given time step.

#### Key identifiers

* **`seq_grp`**: Unique identifier for a continuous tracking sequence of a dolphin group. This is the primary unit for standardising probabilities.
* **`seq`**: Survey or flight identifier.
* **`tm_seq`**: Time index (seconds) within a sequence.
* **`file`**: Drone video file name.
* **`group`**: Dolphin group ID within a sequence.

#### Group composition and structure

* **`pod_size`**: Number of dolphins in the group.
* **`pod_calves`**: Presence/number of calves (if observed).
* **`pod_cohesion`**: Measure of within‑group spacing (lower values = tighter group).

#### Neighbouring dolphin groups

* **`closest_dol_grp_dist`**: Distance to nearest other dolphin group (m).
* **`closest_dol_grp_no`**: ID of the closest neighbouring group.
* **`closest_dol_grp_size`**: Size of the closest neighbouring group.

#### Human activity covariates

* **`closest_boat_dist`**: Distance to nearest tourist boat (m).
* **`closest_boat_type`**: Type of boat (e.g. tourist boat).
* **`n_boats_200m`**: Number of boats within 200 m radius.
* **`closest_net_dist`**: Distance to nearest fishing net (m).

#### Movement metrics

* **`dist`**: Distance moved between time steps (m).
* **`orient`**: Orientation/turning angle (degrees).
* **`speed`**: Dolphin movement speed (m/s).
* **`dist_drone_move`, `orient_drone_move`, `speed_drone_move`**: Drone‑corrected movement metrics.
* **`drone_moving`**: Indicator of drone movement.

#### Behavioural states

* **`behaviour`**: Observed behavioural state (forage, travel, socialise, escape).
* **`behaviour_lag`**: Behaviour at previous time step.
* **`diving` / `diving_behaviour`**: Indicators of diving activity.

#### Group dynamics

* **`split`**: Indicator that a group split occurred.
* **`merge`**: Indicator that a group merge occurred.
* **`tm_split`, `tm_merge`**: Time indices of splitting/merging events.
* **`cluster_prev`, `cluster_next`**: Group identity before/after split or merge.

#### Treatment variables

* **`treatment`**: Presence/type of boats (`no_boats`, `TB`, `FB`, `both`).
* **`treatment_dol`**: Presence of neighbouring dolphin groups.

---

### 2. `Data/dolphin_tracks_processed.csv`

Processed and model-ready dataset derived from `dolphin_tracks.csv`. This file is used directly in all statistical models.

Each row corresponds to a valid modelling observation after filtering, scaling, and factor construction.

#### Key modelling variables

* **`behaviour_f`**: Observed behavioural state at a time step.
* **`behaviour_lag_f`**: Behaviour at previous time step.

#### Covariates (scaled)

* **`pod_size`**: Standardised group size.
* **`closest_boat_dist`**: Standardised distance to nearest boat.
* **`n_boats_200m`**: Standardised boat density.
* **`closest_net_dist`**: Standardised distance to fishing nets.

#### Contextual variables

* **`neighbour`**: Binary indicator of neighbouring dolphin group presence.
* **`season`**: Season category (e.g. pre / post fishing season).
* **`Date`**: Survey date.
* **`cpue`**: Catch per unit effort proxy for fish availability.

#### Sequence identifiers

* **`seq_f`**: Sequence identifiers - primary focal group.
* **`seq_grp_f`**: Sequence group identifier - all the groups tracked across primary group follows.

#### Transition indicator

* **`trans`**:

  * `0` → used for **state probability** models
  * `1` → used for **state transition probability** models

---

## Analytical workflow

1. **Data filtering and scaling**

   * Split processed data into `trans == 0` (state probabilities) and `trans == 1` (transitions)
   * Standardise continuous covariates

2. **Multinomial logistic regression (`mclogit::mblogit`)**

   * Separate models for each treatment (no boats, tourist boats, fishing boats, both)
   * Estimate covariate effects on behavioural state probabilities

3. **Model diagnostics**

   * Confusion matrices and classification accuracy
   * McFadden’s pseudo‑R²
   * Residual diagnostics

4. **Null model comparisons**

   * Likelihood ratio tests against intercept‑only models

5. **Probability estimation**

   * Sequence‑standardised state probabilities
   * Empirical transition probability matrices

6. **Group dynamics models**

   * Binomial GLMs for splitting and merging events
   * Univariate models due to small sample sizes

---

## Outputs

* Estimated coefficients for each multinomial model
* State and transition probability tables
* Inference on how tourism and fisheries influence dolphin behaviour and group structure

---
