## ===============================
## Package installation & loading
## This step is to load all the necessary packages for future code
## ===============================

# CRAN packages
cran_packages <- c(
  "arrow", "tidyr", "purrr", "readr", "igraph", "readxl",
  "httr", "jsonlite", "webchem", "pheatmap", "tibble",
  "stringr", "parallel", "pbmcapply", "furrr", "progressr",
  "ggplot2", "dplyr", "rlang",
  "writexl", "readr","outliers"
)

# Get installed packages once
installed_pkgs <- rownames(installed.packages())

# Install missing CRAN packages
cran_to_install <- setdiff(cran_packages, installed_pkgs)

if (length(cran_to_install) > 0) {
  install.packages(cran_to_install, dependencies = TRUE)
}

# Load all packages
invisible(
  lapply(cran_packages, library, character.only = TRUE)
)

## ===============================
## Functions
## ===============================

# Mapping: original column name -> standardized column name
safe_rename_map <- function() {
  c(
    "Drug/Code" = "Drug",
    "Concentration/xCmax" = "Concentration",
    "ALBUMIN_CV Analysis" = "ALB_CV",
    "ALT_CV Analysis" = "ALT_CV",
    "Percent Viability" = "Viability",
    "Fold Viability (normalized to VEH)" = "Fold_Viability",
    "ATP (% normalized to VEH)" = "ATP_normalized",
    "ALT (ng/day/million cells)" = "ALT",
    "Albumin (µg/day/million cells)" = "ALB",
    "ATP Concentration (µM)" = "ATP"
  )
}

load_and_clean_raw <- function(filepath) {
  raw_data <- readxl::read_xlsx(filepath)
  
  # filter out unwanted Drug/Code categories (only if column exists)
  if ("Drug/Code" %in% names(raw_data)) {
    raw_data <- raw_data %>%
      dplyr::filter(!`Drug/Code` %in% c("OF","UF","Bubble","leak", "Pre-treatment"))
  }
  
  # Safe rename: rename only columns that exist
  mp <- safe_rename_map()
  exist_old <- intersect(names(mp), names(raw_data))
  if (length(exist_old) > 0) {
    raw_data <- raw_data %>%
      dplyr::rename(!!!stats::setNames(exist_old, mp[exist_old]))
  }
  
  # Default drop columns (safe)
  data <- raw_data %>%
    dplyr::select(-dplyr::any_of(c(
      "Type", "Condition", "Donor", "Cell count per chip",
      "ALB_OD1", "ALB_OD2", "ALT_OD1", "ALT_OD2"
    )))
  
  # ------------------------------
  # Build Treatment column
  # ------------------------------
  # If Vehicle control doesn't exist, create it as NA so code won't crash
  if (!("Vehicle control" %in% names(data))) {
    data[["Vehicle control"]] <- NA_character_
  }
  
  data <- data %>%
    dplyr::mutate(
      Drug          = stringr::str_trim(as.character(Drug)),
      Concentration = stringr::str_trim(as.character(Concentration)),
      
      # Remove trailing "XCmax" if present
      conc_clean = dplyr::if_else(
        is.na(Concentration), NA_character_,
        stringr::str_trim(stringr::str_replace(Concentration, "(?i)\\s*XCmax\\s*$", ""))
      ),
      
      # Remove all whitespaces in drug names
      drug_clean = stringr::str_replace_all(Drug, "\\s+", ""),
      
      # Construct treatment label
      treatment = dplyr::case_when(
        is.na(conc_clean) | conc_clean == "" ~
          drug_clean,
        
        `Vehicle control` == "Negative control" ~
          stringr::str_c(
            drug_clean,
            paste0(as.numeric(conc_clean) * 100, "%"),
            sep = "_"
          ),
        
        TRUE ~
          stringr::str_c(drug_clean, conc_clean, sep = "_")
      )
    ) %>%
    dplyr::select(-conc_clean, -drug_clean)
  
    data <- data %>%
    dplyr::mutate(
      ALB_CV = as.numeric(stringr::str_remove(ALB_CV, "%")) ,
      ALT_CV = as.numeric(stringr::str_remove(ALT_CV, "%"))
    )
  
  # ------------------------------
  # Convert assay columns to numeric (safe if missing)
  # ------------------------------
  num_cols <- intersect(names(data), c("ALB", "ALT", "ATP", "Viability", "Fold_Viability", "ATP_normalized"))
  if (length(num_cols) > 0) {
    data <- data %>%
      dplyr::mutate(dplyr::across(dplyr::all_of(num_cols), ~ suppressWarnings(as.numeric(.x))))
  }
  
  data
}


# Assay definitions (standardized names AFTER rename)
assay_definitions <- function() {
  list(
    ALB = c("ALB", "ALB_CV"),
    ALT = c("ALT", "ALT_CV"),
    ATP = c("ATP", "ATP_normalized"),
    Viability = c("Viability", "Fold_Viability")
  )
}

# Detect which assays are available based on existing columns
# Rule: if ANY related column exists, assay is considered available
detect_available_assays <- function(df) {
  defs <- assay_definitions()
  present <- names(df)
  names(defs)[vapply(defs, function(cols) any(cols %in% present), logical(1))]
}

# Suggest columns to keep based on selected assays
suggest_cols_by_assays <- function(df, selected_assays) {
  defs <- assay_definitions()
  present <- names(df)
  wanted <- unique(unlist(defs[selected_assays], use.names = FALSE))
  intersect(present, wanted)
}

# Pretty print columns with numbers
print_numbered <- function(x) {
  for (i in seq_along(x)) cat(sprintf("  [%d] %s\n", i, x[i]))
}


# Outlier detection 
dixon_iter_rm <- function(x, alpha = 0.05, min_n = 3, mad_thresh = 3.5) {
  x <- suppressWarnings(as.numeric(x))
  n <- length(x)
  removed <- rep(FALSE, n)
  
  if (sum(!is.na(x)) < min_n) return(removed)
  
  if (sum(!is.na(x)) >= 25) {
    xx  <- x
    med <- median(xx, na.rm = TRUE)
    md  <- mad(xx, na.rm = TRUE)
    
    if (is.na(md) || md == 0) {
      q <- quantile(xx, probs = c(0.25, 0.75), na.rm = TRUE)
      iqr <- q[2] - q[1]
      if (is.na(iqr) || iqr == 0) return(removed)
      removed <- (xx < (q[1] - 1.5 * iqr)) | (xx > (q[2] + 1.5 * iqr))
    } else {
      z <- abs(xx - med) / md
      removed <- z > mad_thresh
      removed[is.na(removed)] <- FALSE
    }
    return(removed)
  }
  
  repeat {
    idx <- which(!removed & !is.na(x))
    if (length(idx) < min_n) break
    
    xs <- x[idx]
    
    tst <- tryCatch(outliers::dixon.test(xs, two.sided = TRUE),
                    error = function(e) NULL)
    if (is.null(tst)) break
    
    if (!is.na(tst$p.value) && tst$p.value < alpha) {
      alt <- as.character(tst$alternative)
      if (grepl("lowest|min", alt, ignore.case = TRUE)) {
        remove_idx <- idx[which.min(xs)]
      } else if (grepl("highest|max", alt, ignore.case = TRUE)) {
        remove_idx <- idx[which.max(xs)]
      } else {
        break
      }
      removed[remove_idx] <- TRUE
    } else {
      break
    }
  }
  removed
}

# A general pipeline: optional CV filter + day filter + outlier flagging by treatment
run_outlier_pipeline <- function(
    data,
    value_col,        # e.g. "ALB"
    day_value,        # e.g. 14
    cv_col = NULL,    # e.g. "ALB_CV"
    cv_thresh = NULL, # e.g. 0.2
    group_col = "treatment",
    alpha = 0.05,
    min_n = 3,
    mad_thresh = 3.5
) {
  if (!value_col %in% names(data)) stop("value_col not found in data: ", value_col)
  if (!"Day" %in% names(data)) stop("Column 'Day' not found in data.")
  if (!group_col %in% names(data)) stop("group_col not found in data: ", group_col)
  
  df <- data %>%
    dplyr::filter(!is.na(.data[[value_col]]), Day == day_value)
  
  # Optional CV restriction
  if (!is.null(cv_col) && !is.null(cv_thresh)) {
    if (!cv_col %in% names(df)) {
      stop("cv_col not found in data: ", cv_col)
    }
    df <- df %>%
      dplyr::filter(!is.na(.data[[cv_col]]), .data[[cv_col]] <= cv_thresh)
  }
  
  # Flag outliers within each treatment (group_col)
  flagged <- df %>%
    dplyr::group_by(.data[[group_col]]) %>%
    dplyr::group_modify(~ dplyr::mutate(.x, dixon_outlier = dixon_iter_rm(.x[[value_col]], alpha, min_n, mad_thresh))) %>%
    dplyr::ungroup()
  
  outliers_df <- flagged %>% dplyr::filter(dixon_outlier)
  clean_df    <- flagged %>% dplyr::filter(!dixon_outlier) %>% dplyr::select(-dixon_outlier)
  
  list(
    flagged = flagged,
    outliers = outliers_df,
    clean = clean_df
  )
}

make_level_wide <- function(df, value_col, prefix) {
  if (!("treatment" %in% names(df))) stop("Column 'treatment' not found in df.")
  if (!(value_col %in% names(df))) stop("value_col not found in df: ", value_col)
  
  val_sym <- rlang::sym(value_col)
  
  df %>%
    dplyr::mutate(
      Compound = sub("_[^_]*$", "", treatment),
      conc_str = sub(".*_", "", treatment),
      conc_num = suppressWarnings(as.numeric(conc_str)),
      conc_label = ifelse(is.na(conc_num), "Base", paste0(conc_str, "X"))
    ) %>%
    dplyr::arrange(Compound, conc_num) %>%
    dplyr::select(Compound, conc_label, value = !!val_sym) %>%
    tidyr::pivot_wider(
      names_from   = conc_label,
      values_from  = value,
      names_prefix = paste0(prefix, "_"),
      names_sort   = FALSE
    )
}

ratio_overall_by_treatment_ttest <- function(data, response_col) {
  if (is.character(response_col) && length(response_col) == 1) {
    resp_name <- response_col
  } else {
    resp <- rlang::ensym(response_col)
    resp_name <- rlang::as_string(resp)
  }
  
  if (!("treatment" %in% names(data))) stop("Column 'treatment' not found in data.")
  if (!(resp_name %in% names(data))) stop("response_col not found in data: ", resp_name)
  if (!("Drug" %in% names(data))) data[["Drug"]] <- NA_character_
  if (!("Vehicle control" %in% names(data))) data[["Vehicle control"]] <- NA_character_
  
  data <- data %>%
    dplyr::mutate(
      `Vehicle control` = stringr::str_trim(as.character(`Vehicle control`)),
      Drug              = stringr::str_trim(as.character(Drug))
    )
  
  out <- data %>%
    dplyr::group_by(treatment) %>%
    dplyr::group_modify(~{
      x <- .x[[resp_name]]
      x <- x[!is.na(x)]
      n_treat    <- length(x)
      mean_treat <- if (n_treat > 0) mean(x) else NA_real_
      
      vc_raw <- unique(stringr::str_trim(as.character(.x[["Vehicle control"]])))
      vc_raw <- vc_raw[!is.na(vc_raw)]
      
      is_control_treat <- (length(vc_raw) == 1 && vc_raw == "Negative control")
      vc_name <- if (!is_control_treat) vc_raw[vc_raw != "Negative control"][1] else NA_character_
      
      ctrl_df <- if (!is_control_treat && !is.na(vc_name)) {
        data %>%
          dplyr::filter(
            `Vehicle control` == "Negative control",
            treatment == vc_name
          )
      } else {
        data[0, , drop = FALSE]
      }
      
      ctrl_vals <- ctrl_df[[resp_name]]
      ctrl_vals <- ctrl_vals[!is.na(ctrl_vals)]
      n_ctrl    <- length(ctrl_vals)
      mean_ctrl <- if (n_ctrl > 0) mean(ctrl_vals) else NA_real_
      
      if (is_control_treat || n_treat < 2 || n_ctrl < 2) {
        tibble::tibble(
          mean_treat = mean_treat,
          n_treat    = n_treat,
          control_name = if (!is_control_treat) {
            paste0("Vehicle[", vc_name, "]")
          } else {
            "Negative control"
          },
          !!paste0("mean_", resp_name, "_ctrl") := mean_ctrl,
          n_ctrl     = n_ctrl,
          t_stat     = NA_real_,
          p_value    = NA_real_
        )
      } else {
        tt <- try(stats::t.test(x, ctrl_vals, var.equal = FALSE), silent = TRUE)
        if (inherits(tt, "try-error")) {
          tibble::tibble(
            mean_treat = mean_treat,
            n_treat    = n_treat,
            control_name = paste0("Vehicle[", vc_name, "]"),
            !!paste0("mean_", resp_name, "_ctrl") := mean_ctrl,
            n_ctrl     = n_ctrl,
            t_stat     = NA_real_,
            p_value    = NA_real_
          )
        } else {
          tibble::tibble(
            mean_treat = mean_treat,
            n_treat    = n_treat,
            control_name = paste0("Vehicle[", vc_name, "]"),
            !!paste0("mean_", resp_name, "_ctrl") := mean_ctrl,
            n_ctrl     = n_ctrl,
            t_stat     = unname(tt$statistic),
            p_value    = tt$p.value
          )
        }
      }
    }) %>%
    dplyr::ungroup() %>%
    dplyr::rename(!!paste0("mean_", resp_name, "_treat") := mean_treat) %>%
    dplyr::arrange(treatment) %>%
    dplyr::filter(control_name != "Negative control")
  
  out
}

reorder_wide_columns <- function(df) {
  fixed_cols <- "Compound"
  data_cols <- setdiff(names(df), fixed_cols)
  
  sort_values <- sapply(data_cols, function(x) {
    val_str <- sub(".*_([0-9.]+)X$", "\\1", x)
    val_num <- suppressWarnings(as.numeric(val_str))
    ifelse(is.na(val_num), 0, val_num)
  })
  
  sorted_data_cols <- data_cols[order(sort_values)]
  
  df %>%
    dplyr::select(dplyr::all_of(fixed_cols), dplyr::all_of(sorted_data_cols))
}


toxify_t_stat <- function(df, direction = c("UP", "DOWN"),
                          t_col = "t_stat",
                          keep_raw = TRUE,
                          raw_col = "t_stat_raw") {
  
  direction <- match.arg(direction)
  
  if (!(t_col %in% names(df))) stop("toxify_t_stat: column not found: ", t_col)
  
  x <- df[[t_col]]
  
  if (keep_raw && !(raw_col %in% names(df))) {
    df[[raw_col]] <- x
  }
  
  if (direction == "UP") {
    # keep positives; negatives -> 0; NA stays NA
    df[[t_col]] <- ifelse(is.na(x), NA_real_, pmax(x, 0))
  } else {
    # positives -> 0; negatives -> abs; NA stays NA
    df[[t_col]] <- ifelse(is.na(x), NA_real_, pmax(-x, 0))
  }
  
  df
}

