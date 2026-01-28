if (!isatty(stdin())) {
  stop(
    "This script requires an interactive terminal.\n",
    "Please run it in a real Terminal with:\n",
    "  Rscript heatmap.R"
  )
}

source(file.path("R", "functions.R"))

cat("\n=== Heatmap Generator ===\n")

in_dir  <- "data_input"
out_dir <- "output"

if (!dir.exists(in_dir)) stop("data_input/ not found.")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)


# 1) List candidate xlsx files data_input/
xlsx_files <- list.files(in_dir, pattern = "^[^~].*\\.xlsx$", full.names = FALSE)
if (length(xlsx_files) == 0) stop("No .xlsx found in data_input/")

cat("\nAvailable Excel files in data_input/:\n")
print_numbered(xlsx_files)

ans <- trimws(readline("Select file number (default 1): "))
idx <- if (!nzchar(ans)) 1 else suppressWarnings(as.integer(ans))
if (is.na(idx) || idx < 1 || idx > length(xlsx_files)) stop("Invalid selection")

xlsx_path <- file.path(in_dir, xlsx_files[idx])

# 2) Choose sheet
sheets <- readxl::excel_sheets(xlsx_path)
cat("\nSheets:\n")
print_numbered(sheets)

ans <- trimws(readline("Select sheet number (default 1): "))
sidx <- if (!nzchar(ans)) 1 else suppressWarnings(as.integer(ans))
if (is.na(sidx) || sidx < 1 || sidx > length(sheets)) stop("Invalid sheet selection")

sheet_name <- sheets[sidx]

df <- readxl::read_excel(xlsx_path, sheet = sheet_name)

# Basic cleanup: drop Total_Sum if present
if ("Total_Sum" %in% names(df)) {
  df <- dplyr::select(df, -Total_Sum)
}

# 3) Column selection (excluding Compound)
if (!("Compound" %in% names(df))) stop("Selected sheet must include 'Compound' column.")

all_cols <- setdiff(names(df), "Compound")

cat("\nColumns available for heatmap:\n")
print_numbered(all_cols)

cat("\nChoose column mode:\n")
cat("  (A) All columns\n")
cat("  (M) Manual pick by number\n")
mode <- toupper(trimws(readline("Your choice (default A): ")))
if (!nzchar(mode)) mode <- "A"
if (!mode %in% c("A","M")) stop("Invalid choice. Use A/M.")

if (mode == "A") {
  keep_cols <- all_cols
} else {
  raw <- trimws(readline("Select columns by number (e.g., 1,3,5): "))
  idxs <- suppressWarnings(as.integer(trimws(unlist(strsplit(raw, ",")))))
  idxs <- idxs[idxs %in% seq_along(all_cols)]
  if (length(idxs) == 0) stop("No valid columns selected.")
  keep_cols <- all_cols[idxs]
}

heatmap_df <- df %>% dplyr::select(Compound, dplyr::all_of(keep_cols))

# 4) limits + breaks
lim_in <- trimws(readline("Fill limits as min,max (default 0,30): "))
if (!nzchar(lim_in)) lims <- c(0, 30) else {
  lims <- suppressWarnings(as.numeric(trimws(unlist(strsplit(lim_in, ",")))))
  if (length(lims) != 2 || any(is.na(lims))) stop("Invalid limits. Use e.g. 0,30")
}

br_in <- trimws(readline("Rescale breakpoints (e.g., 0,5,30; default 0,5,30): "))
if (!nzchar(br_in)) brks <- c(0, 5, 30) else {
  brks <- suppressWarnings(as.numeric(trimws(unlist(strsplit(br_in, ",")))))
  brks <- brks[!is.na(brks)]
  if (length(brks) < 2) stop("Invalid breakpoints. Need at least 2 numbers.")
}

# 5) Titles + output filename
title <- trimws(readline("Plot title (default 'Compound Toxicity Summary'): "))
if (!nzchar(title)) title <- "Compound Toxicity Summary"

subtitle <- trimws(readline("Plot subtitle (default 'T-score heatmap'): "))
if (!nzchar(subtitle)) subtitle <- "T-score heatmap"

img_name <- trimws(readline("Output image filename (default 'Compound_Toxicity_Heatmap.png'): "))
if (!nzchar(img_name)) img_name <- "Compound_Toxicity_Heatmap.png"
img_path <- file.path(out_dir, img_name)

# 6) Plot & save
p <- plot_tox_heatmap_all(
  heatmap_df,
  title = title,
  subtitle = subtitle,
  limits = lims,
  breaks = brks
)

ggplot2::ggsave(
  filename = img_path,
  plot     = p,
  width    = 10,
  height   = 8,
  dpi      = 300
)

cat("\n✅ Heatmap saved:\n")
cat("  - ", img_path, "\n", sep = "")
