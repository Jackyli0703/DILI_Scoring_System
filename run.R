if (!isatty(stdin())) {
  stop(
    "This script requires an interactive terminal.\n",
    "Please run it in a real Terminal with:\n",
    "  Rscript run.R"
  )
}

# Entry point — NO library() here
source(file.path("R", "functions.R"))

cat("\n=== Scoring System Program ===\n")

in_dir  <- "data_input"
out_dir <- "output"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# ---- list xlsx but ignore Excel temp files (~$...) ----
files <- list.files(in_dir, pattern = "^[^~].*\\.xlsx$", full.names = FALSE)
if (length(files) == 0) stop("No Excel file found in data_input/")

cat("\nAvailable input files:\n")
print_numbered(files)

ans <- trimws(readline("Select file number (default 1): "))
if (!nzchar(ans)) {
  idx <- 1
} else {
  idx <- suppressWarnings(as.integer(ans))
}
if (is.na(idx) || idx < 1 || idx > length(files)) stop("Invalid selection")

input_file <- files[idx]
input_path <- file.path(in_dir, input_file)

# ---- load + clean (includes treatment build & numeric conversion) ----
data <- load_and_clean_raw(input_path)

# ==============================
# Part 1: Output processed_data (always)
# ==============================
available_assays <- detect_available_assays(data)
picked_assays <- character(0)

if (length(available_assays) > 0) {
  cat("\nDetected assays:\n")
  print_numbered(available_assays)
  
  ans <- trimws(readline("Select assays (e.g. 1,3) or press Enter to skip: "))
  if (nzchar(ans)) {
    idxs <- suppressWarnings(as.integer(trimws(unlist(strsplit(ans, ",")))))
    idxs <- idxs[idxs %in% seq_along(available_assays)]
    picked_assays <- available_assays[idxs]
  }
}

all_cols <- names(data)
base_cols <- intersect(
  all_cols,
  c("Plate", "Device", "Chip", "Day", "Drug", "Concentration", "Vehicle control", "treatment")
)

assay_cols <- if (length(picked_assays) > 0) {
  suggest_cols_by_assays(data, picked_assays)
} else character(0)

suggested <- unique(c(base_cols, assay_cols))

cat("\n=== Column options preview ===\n")

cat("\n[A] All columns (", length(all_cols), "):\n", sep = "")
print_numbered(all_cols)

cat("\n[S] Suggested columns (", length(suggested), "):\n", sep = "")
if (length(suggested) == 0) {
  cat("  (none)\n")
} else {
  print_numbered(suggested)
}

cat("\nChoose column mode:\n")
cat("  (A) All columns\n")
cat("  (S) Suggested columns\n")
cat("  (M) Manual (pick by number from ALL columns list above)\n")

mode <- toupper(trimws(readline("Your choice (default S): ")))
if (!nzchar(mode)) mode <- "S"
if (!mode %in% c("A","S","M")) stop("Invalid choice. Use A/S/M.")

if (mode == "S") {
  final_cols <- suggested
} else if (mode == "M") {
  keep_raw <- trimws(readline("Select columns by number (e.g., 1,2,5): "))
  keep <- suppressWarnings(as.integer(trimws(unlist(strsplit(keep_raw, ",")))))
  keep <- keep[keep %in% seq_along(all_cols)]
  if (length(keep) == 0) stop("No valid columns selected.")
  final_cols <- all_cols[keep]
} else {
  final_cols <- all_cols
}

final_data <- data[, final_cols, drop = FALSE]

out_base <- sub("\\.xlsx$", "_processed_data", input_file)
out_xlsx <- file.path(out_dir, paste0(out_base, ".xlsx"))
out_csv  <- file.path(out_dir, paste0(out_base, ".csv"))

writexl::write_xlsx(final_data, out_xlsx)
readr::write_csv(final_data, out_csv)

cat("\n✅ Processed data saved:\n")
cat("  - ", out_xlsx, "\n", sep = "")
cat("  - ", out_csv,  "\n", sep = "")

# ==============================
# Part 2: Batch outlier detection (NO data carryover)
# - Multiple days per assay supported
# - All runs use ORIGINAL cleaned data with treatment (data)
# - Prevent duplicate assay+day runs
# - For each (assay, day): output outliers + clean + summary (ratio_overall + wide t_stat)
# - Combine all wide t_stat tables into one combined file
# ==============================
cat("\n=== Batch outlier detection ===\n")
cat("All outlier runs are based on the ORIGINAL cleaned data with treatment.\n")
cat("No run will modify the dataset for subsequent runs.\n")

do_batch <- toupper(trimws(readline("Run batch outlier detection? (Y/N, default N): ")))
if (!nzchar(do_batch)) do_batch <- "N"

if (do_batch == "Y") {
  
  # Ask ONCE: reorder wide columns?
  reorder_all <- toupper(trimws(readline("Reorder ALL wide t_stat tables by concentration? (Y/N, default N): ")))
  if (!nzchar(reorder_all)) reorder_all <- "N"
  
  # Candidate variables present in data
  candidate_vars <- intersect(names(data), c("ALB", "ALT", "ATP", "Viability"))
  if (length(candidate_vars) == 0) stop("No candidate assay columns found (ALB/ALT/ATP/Viability).")
  
  cat("\nAvailable assays for outlier detection (present in this file):\n")
  print_numbered(candidate_vars)
  
  var_ans <- trimws(readline("Select assays by number (e.g., 1,3,4). Press Enter to cancel: "))
  if (!nzchar(var_ans)) {
    cat("No assays selected. Skip outlier detection.\n")
  } else {
    
    var_idxs <- suppressWarnings(as.integer(trimws(unlist(strsplit(var_ans, ",")))))
    var_idxs <- var_idxs[var_idxs %in% seq_along(candidate_vars)]
    if (length(var_idxs) == 0) stop("No valid assay selection.")
    
    selected_vars <- candidate_vars[var_idxs]
    tox_dir <- list()
    
    cat("\nToxicity direction per assay:\n")
    cat("  UP   = higher value is more toxic (keep +t, set -t to 0)\n")
    cat("  DOWN = lower value is more toxic (set +t to 0, abs(-t))\n")
    
    for (a in selected_vars) {
      ans_dir <- toupper(trimws(readline(paste0("For ", a, ", which direction indicates toxicity? (UP/DOWN, default UP): "))))
      if (!nzchar(ans_dir)) ans_dir <- "UP"
      if (!ans_dir %in% c("UP", "DOWN")) ans_dir <- "UP"
      tox_dir[[a]] <- ans_dir
    }
    # Day choices available in the dataset
    if (!("Day" %in% names(data))) stop("Column 'Day' not found; cannot run outlier detection.")
    day_choices <- sort(unique(data$Day))
    cat("\nAvailable Day values in data:\n")
    print(day_choices)
    
    # Prevent duplicates across the whole batch
    processed_keys <- character(0)
    
    # Collect wide t_stat tables
    wide_tables <- list()
    wide_keys <- character(0)
    
    # CV mapping (only ALB/ALT typically)
    cv_map <- list(ALB = "ALB_CV", ALT = "ALT_CV")
    
    for (value_col in selected_vars) {
      
      # ---- multiple days per assay ----
      day_ans <- trimws(readline(paste0("Days for ", value_col, " (comma-separated, e.g., 11,14): ")))
      if (!nzchar(day_ans)) {
        cat("No Day provided for ", value_col, ". Skipping.\n", sep = "")
        next
      }
      
      days_vec <- suppressWarnings(as.integer(trimws(unlist(strsplit(day_ans, ",")))))
      days_vec <- days_vec[!is.na(days_vec)]
      if (length(days_vec) == 0) {
        cat("Invalid Day input for ", value_col, ". Skipping.\n", sep = "")
        next
      }
      
      # Keep only days that exist in data
      days_vec <- unique(days_vec)
      valid_days <- days_vec[days_vec %in% day_choices]
      invalid_days <- setdiff(days_vec, valid_days)
      
      if (length(invalid_days) > 0) {
        cat("⚠️ Days not found in data for ", value_col, ": ", paste(invalid_days, collapse = ", "), "\n", sep = "")
      }
      if (length(valid_days) == 0) {
        cat("No valid Day values for ", value_col, ". Skipping.\n", sep = "")
        next
      }
      
      # ---- CV filter per assay (applies to all selected days for that assay) ----
      cv_col <- cv_map[[value_col]]
      use_cv <- "N"
      cv_thresh <- NULL
      
      if (!is.null(cv_col) && cv_col %in% names(data)) {
        use_cv <- toupper(trimws(readline(paste0("Apply CV filter for ", value_col, " using ", cv_col, "? (Y/N, default N): "))))
        if (!nzchar(use_cv)) use_cv <- "N"
        
        if (use_cv == "Y") {
          thr_ans <- trimws(readline("CV threshold (e.g., 0.2): "))
          cv_thresh <- suppressWarnings(as.numeric(thr_ans))
          if (is.na(cv_thresh)) {
            cat("Invalid CV threshold. Skipping CV filter for ", value_col, ".\n", sep = "")
            use_cv <- "N"
            cv_thresh <- NULL
          }
        }
      } else {
        cat("(No CV column for ", value_col, "; skipping CV filter.)\n", sep = "")
      }
      
      # ---- run outlier for each valid day (no carryover) ----
      for (day_value in valid_days) {
        
        key <- paste0(value_col, "_Day", day_value)
        if (key %in% processed_keys) {
          cat("⚠️ Already processed ", key, ". Skipping duplicate.\n", sep = "")
          next
        }
        processed_keys <- c(processed_keys, key)
        
        res_out <- run_outlier_pipeline(
          data = data,  # ORIGINAL cleaned data with treatment
          value_col = value_col,
          day_value = day_value,
          cv_col = if (use_cv == "Y") cv_col else NULL,
          cv_thresh = if (use_cv == "Y") cv_thresh else NULL
        )
        
        # Output files: outliers + clean
        out_tag <- paste0(
          tools::file_path_sans_ext(basename(input_file)),
          "_outlier_", value_col, "_Day", day_value,
          if (use_cv == "Y") paste0("_cv", cv_thresh) else ""
        )
        
        out_xlsx2 <- file.path(out_dir, paste0(out_tag, ".xlsx"))
        writexl::write_xlsx(
          list(
            outliers = res_out$outliers,
            clean_data = res_out$clean
          ),
          out_xlsx2
        )
        
        # ---- t-test summary + wide t_stat table ----
        tt_df <- ratio_overall_by_treatment_ttest(res_out$clean, value_col)
        
        # ---- NEW: Align t_stat to toxicity direction for this assay ----
        dir_now <- tox_dir[[value_col]]
        if (is.null(dir_now)) dir_now <- "UP"  # fallback
        tt_df <- toxify_t_stat(tt_df, direction = dir_now, t_col = "t_stat", keep_raw = TRUE)
        prefix <- paste0(toupper(value_col), "_DAY", day_value)
        wide_t <- make_level_wide(tt_df, "t_stat", prefix)
        
        if (reorder_all == "Y") {
          wide_t <- reorder_wide_columns(wide_t)
        }
        
        out_summary_xlsx <- file.path(out_dir, paste0(out_tag, "_summary.xlsx"))
        writexl::write_xlsx(
          list(
            ratio_overall = tt_df,
            wide_tstat    = wide_t
          ),
          out_summary_xlsx
        )
        
        # Collect for combined output
        wide_key <- paste0(prefix, if (use_cv == "Y") paste0("_cv", cv_thresh) else "")
        if (!(wide_key %in% wide_keys)) {
          wide_tables[[wide_key]] <- wide_t
          wide_keys <- c(wide_keys, wide_key)
        }
        
        cat("\n✅ Completed ", key, "\n", sep = "")
        cat("  - ", out_xlsx2, "\n", sep = "")
        cat("  - ", out_summary_xlsx, "\n", sep = "")
        cat("  Outliers rows: ", nrow(res_out$outliers), " | Clean rows: ", nrow(res_out$clean), "\n", sep = "")
      } # end for day_value
    }   # end for value_col
    
    cat("\nBatch outlier detection finished.\n")
    
    # ---- Combine all wide tables into one ----
    if (length(wide_tables) > 0) {
      combined_wide_T <- Reduce(
        function(x, y) dplyr::full_join(x, y, by = "Compound"),
        wide_tables
      )
      
      combined_tag <- paste0(
        tools::file_path_sans_ext(basename(input_file)),
        "_combined_wide_tstat"
      )
      combined_xlsx <- file.path(out_dir, paste0(combined_tag, ".xlsx"))
      writexl::write_xlsx(combined_wide_T, combined_xlsx)
      
      cat("\n✅ Combined wide t_stat table saved:\n")
      cat("  - ", combined_xlsx, "\n", sep = "")
    } else {
      cat("\n(No wide tables generated; combined output skipped.)\n")
    }
    
  } # end else (assays selected)
}   # end if do_batch == "Y"




    
