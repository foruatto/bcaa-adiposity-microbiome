argv <- commandArgs(trailingOnly = TRUE)

usage_text <- paste(
  "Usage:",
  "  Rscript UKB_analysis.R --input-dir <directory> --output-dir <directory>",
  "",
  "Required files:",
  "  baseline.csv or baseline.tsv",
  "  diagnoses.csv or diagnoses.tsv",
  "  imaging.csv or imaging.tsv",
  "",
  "Optional files:",
  "  diet.csv or diet.tsv",
  "  medication.csv or medication.tsv",
  "",
  "The diagnosis file must contain eid, icd10, and diagnosis_date.",
  "Dates must use YYYY-MM-DD. Imaging masses must be in grams, ASAT in litres, and liver fat in percent.",
  sep = "\n"
)

read_options <- function(x) {
  if (length(x) == 1L && x == "--help") {
    cat(usage_text, "\n")
    quit(save = "no", status = 0L)
  }
  keys <- c("--input-dir", "--output-dir")
  if (length(x) != 4L || !all(keys %in% x[c(1L, 3L)])) stop(usage_text, call. = FALSE)
  values <- setNames(x[c(2L, 4L)], x[c(1L, 3L)])
  list(input = values[["--input-dir"]], output = values[["--output-dir"]])
}

options_run <- read_options(argv)

required_packages <- c("survival", "cmprsk", "ggplot2")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1L), quietly = TRUE)]
if (length(missing_packages)) stop("Missing packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)

if (!dir.exists(options_run$input)) stop("Input directory does not exist.", call. = FALSE)
if (dir.exists(options_run$output) && length(list.files(options_run$output, all.files = TRUE, no.. = TRUE))) stop("Output directory must be empty.", call. = FALSE)
dir.create(options_run$output, recursive = TRUE, showWarnings = FALSE)

locate_table <- function(folder, stem, optional = FALSE) {
  candidates <- file.path(folder, paste0(stem, c(".csv", ".tsv")))
  found <- candidates[file.exists(candidates)]
  if (!length(found) && optional) return(NULL)
  if (length(found) != 1L) stop("Expected exactly one input table for ", stem, ".", call. = FALSE)
  found
}

load_table <- function(path) {
  separator <- if (grepl("[.]tsv$", path, ignore.case = TRUE)) "\t" else ","
  utils::read.table(path, header = TRUE, sep = separator, quote = "\"", comment.char = "", check.names = FALSE, stringsAsFactors = FALSE, na.strings = c("", "NA"))
}

assert_fields <- function(x, fields, label) {
  absent <- setdiff(fields, names(x))
  if (length(absent)) stop(label, " is missing: ", paste(absent, collapse = ", "), call. = FALSE)
}

assert_ids <- function(x, label) {
  assert_fields(x, "eid", label)
  x$eid <- as.character(x$eid)
  if (anyNA(x$eid) || any(!nzchar(x$eid)) || anyDuplicated(x$eid)) stop(label, " must contain unique nonmissing eid values.", call. = FALSE)
  x
}

merge_one <- function(left, right, label) {
  right <- assert_ids(right, label)
  overlap <- intersect(setdiff(names(left), "eid"), setdiff(names(right), "eid"))
  if (length(overlap)) stop("Duplicated fields in ", label, ": ", paste(overlap, collapse = ", "), call. = FALSE)
  merge(left, right, by = "eid", all.x = TRUE, sort = FALSE)
}

save_csv <- function(x, name) utils::write.csv(x, file.path(options_run$output, name), row.names = FALSE, na = "")

parse_date <- function(x, field) {
  result <- as.Date(x, format = "%Y-%m-%d")
  invalid <- !is.na(x) & nzchar(as.character(x)) & is.na(result)
  if (any(invalid)) stop("Invalid YYYY-MM-DD values in ", field, ".", call. = FALSE)
  result
}

standardize_code <- function(x) toupper(gsub("[^A-Z0-9]", "", as.character(x)))

first_diagnosis <- function(records, ids, expression) {
  code <- standardize_code(records$icd10)
  date <- parse_date(records$diagnosis_date, "diagnosis_date")
  selected <- grepl(expression, code) & !is.na(records$eid)
  flag_ids <- unique(as.character(records$eid[selected]))
  dated <- selected & !is.na(date)
  earliest <- tapply(as.numeric(date[dated]), as.character(records$eid[dated]), min)
  position <- match(ids, names(earliest))
  first_date <- as.Date(as.numeric(earliest[position]), origin = "1970-01-01")
  data.frame(flag = as.integer(ids %in% flag_ids), date = first_date)
}

collapse_diagnoses <- function(records, ids) {
  assert_fields(records, c("eid", "icd10", "diagnosis_date"), "diagnoses")
  records$eid <- as.character(records$eid)
  ihd <- first_diagnosis(records, ids, "^I2[0-5]")
  diabetes <- first_diagnosis(records, ids, "^E11")
  liver <- first_diagnosis(records, ids, "^(K760|K758)")
  data.frame(eid = ids, ihd = ihd$flag, ihd_date = ihd$date, diabetes = diabetes$flag, diabetes_date = diabetes$date, liver = liver$flag, liver_date = liver$date)
}

weekly_frequency <- function(x) {
  y <- trimws(as.character(x))
  dictionary <- c("Never" = 0, "Less than once a week" = 0.5, "Once a week" = 1, "2-4 times a week" = 3, "5-6 times a week" = 5.5, "Once or more daily" = 7)
  numeric_value <- suppressWarnings(as.numeric(y))
  mapped <- unname(dictionary[y])
  ifelse(!is.na(numeric_value), numeric_value, mapped)
}

numeric_quantity <- function(x) {
  y <- trimws(as.character(x))
  y[y == "Less than one"] <- "0.5"
  suppressWarnings(as.numeric(y))
}

build_diet_score <- function(x) {
  fields <- c("eid", "fresh_fruit_day", "dried_fruit_day", "cooked_vegetables_day", "raw_vegetables_day", "bread_week", "bread_type", "cereal_week", "cereal_type", "oily_fish", "non_oily_fish", "cheese", "milk_type", "spread_type", "spread_detail", "processed_meat", "poultry", "beef", "lamb", "pork", "sugar_free_selection")
  assert_fields(x, fields, "diet")
  x <- assert_ids(x, "diet")
  fruit <- numeric_quantity(x$fresh_fruit_day) + 0.5 * numeric_quantity(x$dried_fruit_day)
  vegetables <- numeric_quantity(x$cooked_vegetables_day) + numeric_quantity(x$raw_vegetables_day)
  bread_day <- numeric_quantity(x$bread_week) / 7
  cereal_day <- numeric_quantity(x$cereal_week) / 7
  whole_bread <- ifelse(x$bread_type %in% c("Wholemeal", "Wholemeal or wholegrain"), bread_day, 0)
  whole_cereal <- ifelse(x$cereal_type %in% c("Bran cereal", "Oat cereal", "Muesli"), cereal_day, 0)
  fish <- weekly_frequency(x$oily_fish) + weekly_frequency(x$non_oily_fish)
  dairy <- weekly_frequency(x$cheese) + as.integer(x$milk_type %in% c("Semi-skimmed", "Skimmed")) * 7
  healthy_spread <- x$spread_type %in% c("Flora Pro-Active/Benecol") | x$spread_detail %in% c("Flora Pro-Active or Benecol", "Olive oil based spread", "Polyunsaturated/sunflower oil based spread", "Other low or reduced fat spread", "Soft margarine")
  refined_bread <- ifelse(x$bread_type %in% c("White", "Brown", "Other"), bread_day, 0)
  refined_cereal <- ifelse(x$cereal_type %in% c("Biscuit cereal", "Other"), cereal_day, 0)
  unprocessed <- weekly_frequency(x$poultry) + weekly_frequency(x$beef) + weekly_frequency(x$lamb) + weekly_frequency(x$pork)
  components <- cbind(fruit >= 3, vegetables >= 3, whole_bread + whole_cereal >= 3, fish >= 2, dairy >= 3, bread_day >= 2 & healthy_spread, refined_bread + refined_cereal <= 2, weekly_frequency(x$processed_meat) <= 1, unprocessed <= 2, as.character(x$sugar_free_selection) %in% c("1", "yes", "Yes", "TRUE"))
  components[is.na(components)] <- FALSE
  data.frame(eid = x$eid, diet_score = rowSums(components))
}

medicine_dictionary <- list(
  antihypertensive = c("lisinopril", "captopril", "enalapril", "ramipril", "perindopril", "losartan", "valsartan", "irbesartan", "candesartan", "telmisartan", "olmesartan", "amlodipine", "felodipine", "nifedipine", "verapamil", "diltiazem", "atenolol", "bisoprolol", "metoprolol", "carvedilol", "propranolol", "indapamide", "hydrochlorothiazide", "chlortalidone", "furosemide", "bumetanide", "spironolactone"),
  glucose_lowering = c("metformin", "insulin", "gliclazide", "glipizide", "glibenclamide", "glimepiride", "tolbutamide", "pioglitazone", "rosiglitazone"),
  lipid_lowering = c("atorvastatin", "rosuvastatin", "simvastatin", "fluvastatin", "pravastatin", "bezafibrate", "fenofibrate", "clofibrate", "gemfibrozil", "cholestyramine", "colestipol", "acipimox")
)

detect_medicine <- function(x) {
  x <- assert_ids(x, "medication")
  columns <- setdiff(names(x), "eid")
  if (!length(columns)) stop("Medication input has no medicine columns.", call. = FALSE)
  joined <- apply(x[columns], 1L, function(z) paste(tolower(z[!is.na(z)]), collapse = " "))
  detect <- function(terms) as.integer(vapply(joined, function(z) any(vapply(terms, grepl, logical(1L), x = z, fixed = TRUE)), logical(1L)))
  blood_pressure <- detect(medicine_dictionary$antihypertensive)
  diabetes <- detect(medicine_dictionary$glucose_lowering)
  lipid <- detect(medicine_dictionary$lipid_lowering)
  data.frame(eid = x$eid, antihypertensive_use = blood_pressure, diabetes_medication_use = diabetes, lipid_medication_use = lipid, medication_use = as.integer(blood_pressure + diabetes + lipid > 0))
}

estimate_egfr <- function(creatinine, sex, age) {
  scr <- suppressWarnings(as.numeric(creatinine)) / 88.4
  female <- tolower(as.character(sex)) %in% c("female", "f", "0")
  male <- tolower(as.character(sex)) %in% c("male", "m", "1")
  valid <- is.finite(scr) & scr > 0 & is.finite(age) & age > 0 & (female | male)
  result <- rep(NA_real_, length(scr))
  kappa <- ifelse(female, 0.7, 0.9)
  alpha <- ifelse(female, -0.329, -0.411)
  coefficient <- ifelse(female, 144, 141)
  result[valid] <- coefficient[valid] * pmin(scr[valid] / kappa[valid], 1)^alpha[valid] * pmax(scr[valid] / kappa[valid], 1)^-1.209 * 0.993^age[valid]
  result
}

convert_covariates <- function(x) {
  factor_fields <- intersect(c("sex", "smoking", "drinking", "physical_activity", "education", "ethnicity", "medication_use"), names(x))
  for (field in factor_fields) x[[field]] <- factor(x[[field]])
  x$age_group <- factor(ifelse(x$age < 60, "Under 60", "60 or older"), levels = c("Under 60", "60 or older"))
  x$healthy_diet <- factor(ifelse(x$diet_score >= 5, "Yes", "No"), levels = c("No", "Yes"))
  x
}

make_exposures <- function(x, exposures) {
  details <- vector("list", length(exposures))
  for (i in seq_along(exposures)) {
    field <- exposures[i]
    deviation <- stats::sd(x[[field]], na.rm = TRUE)
    if (!is.finite(deviation) || deviation <= 0) stop("Invalid exposure standard deviation for ", field, ".", call. = FALSE)
    boundaries <- stats::quantile(x[[field]], seq(0, 1, 0.25), na.rm = TRUE, type = 7)
    if (anyDuplicated(boundaries)) stop("Nonunique quartile boundaries for ", field, ".", call. = FALSE)
    x[[paste0(field, "_per_sd")]] <- x[[field]] / deviation
    x[[paste0(field, "_quartile")]] <- cut(x[[field]], boundaries, include.lowest = TRUE, labels = paste0("Q", 1:4))
    details[[i]] <- data.frame(exposure = field, standard_deviation = deviation, minimum = boundaries[1], q1 = boundaries[2], median = boundaries[3], q3 = boundaries[4], maximum = boundaries[5])
  }
  attr(x, "exposure_details") <- do.call(rbind, details)
  x
}

baseline <- assert_ids(load_table(locate_table(options_run$input, "baseline")), "baseline")
diagnosis_records <- load_table(locate_table(options_run$input, "diagnoses"))
imaging <- assert_ids(load_table(locate_table(options_run$input, "imaging")), "imaging")

diet_path <- locate_table(options_run$input, "diet", optional = TRUE)
medication_path <- locate_table(options_run$input, "medication", optional = TRUE)
if (!is.null(diet_path)) baseline <- merge_one(baseline, build_diet_score(load_table(diet_path)), "diet")
if (!is.null(medication_path)) baseline <- merge_one(baseline, detect_medicine(load_table(medication_path)), "medication")
if (!"diet_score" %in% names(baseline)) stop("diet_score must exist in baseline when diet input is absent.", call. = FALSE)
if (!"medication_use" %in% names(baseline)) stop("medication_use must exist in baseline when medication input is absent.", call. = FALSE)

baseline <- merge_one(baseline, collapse_diagnoses(diagnosis_records, baseline$eid), "diagnosis summary")

exposures <- c("bcaa", "valine", "leucine", "isoleucine")
social_covariates <- c("sex", "age", "smoking", "drinking", "physical_activity", "education", "deprivation", "diet_score", "ethnicity")
clinical_covariates <- c("egfr", "ldl_cholesterol", "hba1c", "medication_use")
baseline_required <- unique(c("eid", exposures, social_covariates, "bmi", "creatinine", "fasting_glucose", "baseline_date", "death_date", "ihd", "ihd_date", "diabetes", "diabetes_date", "liver", "liver_date", "medication_use"))
assert_fields(baseline, baseline_required, "combined baseline data")

baseline$age <- suppressWarnings(as.numeric(baseline$age))
baseline$egfr <- estimate_egfr(baseline$creatinine, baseline$sex, baseline$age)
date_fields <- c("baseline_date", "death_date", "ihd_date", "diabetes_date", "liver_date")
for (field in date_fields) baseline[[field]] <- parse_date(baseline[[field]], field)

flow <- data.frame(step = "Input", participants = nrow(baseline))
complete_fields <- unique(c(exposures, social_covariates, clinical_covariates, "bmi", "fasting_glucose", "baseline_date"))
eligible <- baseline[stats::complete.cases(baseline[complete_fields]), , drop = FALSE]
flow <- rbind(flow, data.frame(step = "Complete analysis variables", participants = nrow(eligible)))

administrative_end <- as.Date("2023-03-31")
eligible$followup_end <- pmin(eligible$death_date, administrative_end, na.rm = TRUE)
eligible <- eligible[eligible$followup_end > eligible$baseline_date, , drop = FALSE]
flow <- rbind(flow, data.frame(step = "Positive follow-up", participants = nrow(eligible)))

prevalent <- (!is.na(eligible$ihd_date) & eligible$ihd_date <= eligible$baseline_date) | (!is.na(eligible$diabetes_date) & eligible$diabetes_date <= eligible$baseline_date) | (!is.na(eligible$liver_date) & eligible$liver_date <= eligible$baseline_date)
eligible <- eligible[!prevalent & eligible$hba1c < 48 & eligible$fasting_glucose < 7, , drop = FALSE]
flow <- rbind(flow, data.frame(step = "Incident-disease cohort", participants = nrow(eligible)))
if (nrow(eligible) < 20L) stop("Fewer than 20 eligible participants.", call. = FALSE)

endpoints <- c("ihd", "diabetes", "liver")
for (endpoint in endpoints) {
  diagnosis_date <- eligible[[paste0(endpoint, "_date")]]
  end_date <- pmin(diagnosis_date, eligible$death_date, administrative_end, na.rm = TRUE)
  eligible[[paste0(endpoint, "_event")]] <- as.integer(!is.na(diagnosis_date) & diagnosis_date <= end_date)
  eligible[[paste0(endpoint, "_years")]] <- as.numeric(end_date - eligible$baseline_date) / 365
}

eligible <- make_exposures(convert_covariates(eligible), exposures)
save_csv(flow, "cohort_flow.csv")
save_csv(attr(eligible, "exposure_details"), "exposure_scaling.csv")

extract_fit <- function(model, exponentiate) {
  coefficients <- stats::coef(model)
  errors <- sqrt(diag(stats::vcov(model)))
  critical <- if (inherits(model, "lm")) stats::qt(0.975, stats::df.residual(model)) else stats::qnorm(0.975)
  probability <- if (inherits(model, "lm")) 2 * stats::pt(-abs(coefficients / errors), stats::df.residual(model)) else 2 * stats::pnorm(-abs(coefficients / errors))
  lower <- coefficients - critical * errors
  upper <- coefficients + critical * errors
  data.frame(term = names(coefficients), coefficient = unname(coefficients), standard_error = unname(errors), estimate = if (exponentiate) exp(coefficients) else coefficients, lower_95 = if (exponentiate) exp(lower) else lower, upper_95 = if (exponentiate) exp(upper) else upper, p_value = unname(probability), sample_size = stats::nobs(model), row.names = NULL)
}

fit_cox <- function(data, endpoint, exposure, covariates) {
  formula <- stats::reformulate(c(exposure, covariates), response = paste0("survival::Surv(", endpoint, "_years, ", endpoint, "_event)"))
  survival::coxph(formula, data = data, ties = "efron", x = TRUE, model = TRUE)
}

model_sets <- list(Model_1 = character(), Model_2 = social_covariates, Model_3 = c(social_covariates, clinical_covariates))
cox_rows <- list()
for (endpoint in endpoints) for (parent in exposures) for (exposure in paste0(parent, c("_per_sd", "_quartile"))) for (model_name in names(model_sets)) {
  fit <- fit_cox(eligible, endpoint, exposure, model_sets[[model_name]])
  result <- extract_fit(fit, TRUE)
  result <- result[startsWith(result$term, exposure), , drop = FALSE]
  result$endpoint <- endpoint
  result$exposure <- exposure
  result$model <- model_name
  cox_rows[[length(cox_rows) + 1L]] <- result
}
cox_results <- do.call(rbind, cox_rows)
save_csv(cox_results, "cox_models.csv")

subgroup_fields <- c("sex", "smoking", "drinking", "physical_activity", "age_group", "healthy_diet")
subgroup_rows <- list()
for (endpoint in endpoints) for (parent in exposures) for (stratum in subgroup_fields) {
  exposure <- paste0(parent, "_per_sd")
  adjusted <- setdiff(c(social_covariates, clinical_covariates), c(stratum, if (stratum == "age_group") "age" else character(), if (stratum == "healthy_diet") "diet_score" else character()))
  variables <- unique(c(paste0(endpoint, c("_years", "_event")), exposure, stratum, adjusted))
  analysis_data <- droplevels(eligible[stats::complete.cases(eligible[variables]), variables, drop = FALSE])
  interaction_formula <- stats::as.formula(paste0("survival::Surv(", endpoint, "_years, ", endpoint, "_event) ~ ", exposure, " * ", stratum, " + ", paste(adjusted, collapse = " + ")))
  additive_formula <- stats::as.formula(paste0("survival::Surv(", endpoint, "_years, ", endpoint, "_event) ~ ", exposure, " + ", stratum, " + ", paste(adjusted, collapse = " + ")))
  p_interaction <- tryCatch({
    comparison <- stats::anova(survival::coxph(additive_formula, analysis_data), survival::coxph(interaction_formula, analysis_data), test = "Chisq")
    probability_column <- grep("^P|Pr", names(comparison), value = TRUE)
    if (!length(probability_column)) NA_real_ else as.numeric(comparison[2, probability_column[1]])
  }, error = function(e) NA_real_)
  for (level in levels(factor(analysis_data[[stratum]]))) {
    section <- droplevels(analysis_data[as.character(analysis_data[[stratum]]) == level, , drop = FALSE])
    usable <- adjusted[vapply(section[adjusted], function(z) length(unique(z)) > 1L, logical(1L))]
    fitted <- tryCatch(fit_cox(section, endpoint, exposure, usable), error = function(e) NULL)
    if (is.null(fitted)) next
    result <- extract_fit(fitted, TRUE)
    result <- result[result$term == exposure, , drop = FALSE]
    if (!nrow(result)) next
    result$endpoint <- endpoint
    result$exposure <- exposure
    result$stratum <- stratum
    result$level <- level
    result$p_interaction <- p_interaction
    subgroup_rows[[length(subgroup_rows) + 1L]] <- result
  }
}
subgroup_results <- if (length(subgroup_rows)) do.call(rbind, subgroup_rows) else data.frame()
save_csv(subgroup_results, "subgroup_models.csv")

if (nrow(subgroup_results)) {
  plot_data <- subgroup_results[is.finite(subgroup_results$estimate) & is.finite(subgroup_results$lower_95) & is.finite(subgroup_results$upper_95) & subgroup_results$estimate > 0 & subgroup_results$lower_95 > 0 & subgroup_results$upper_95 > 0, , drop = FALSE]
  plot_data$label <- paste(plot_data$endpoint, plot_data$exposure, plot_data$stratum, plot_data$level, sep = " | ")
  plot_data$label <- factor(plot_data$label, levels = rev(plot_data$label))
  forest <- ggplot2::ggplot(plot_data, ggplot2::aes(estimate, label)) + ggplot2::geom_vline(xintercept = 1, linetype = 2) + ggplot2::geom_errorbar(ggplot2::aes(xmin = lower_95, xmax = upper_95), width = 0, orientation = "y") + ggplot2::geom_point() + ggplot2::scale_x_log10() + ggplot2::labs(x = "Hazard ratio", y = NULL) + ggplot2::theme_bw()
  ggplot2::ggsave(file.path(options_run$output, "subgroup_forest.pdf"), forest, width = 10, height = max(7, nrow(plot_data) * 0.16), limitsize = FALSE)
}

spline_rows <- list()
for (endpoint in endpoints) for (parent in exposures) {
  knots <- stats::quantile(eligible[[parent]], c(0.1, 0.5, 0.9), na.rm = TRUE, type = 7)
  if (anyDuplicated(knots)) stop("Nonunique spline knots for ", parent, ".", call. = FALSE)
  basis <- splines::ns(eligible[[parent]], knots = knots[2], Boundary.knots = knots[c(1, 3)])
  spline_names <- paste0("spline_", seq_len(ncol(basis)))
  analysis_data <- eligible
  analysis_data[spline_names] <- basis
  full <- fit_cox(analysis_data, endpoint, spline_names, c(social_covariates, clinical_covariates))
  linear <- fit_cox(analysis_data, endpoint, parent, c(social_covariates, clinical_covariates))
  null <- fit_cox(analysis_data, endpoint, character(), c(social_covariates, clinical_covariates))
  overall <- stats::anova(null, full, test = "Chisq")
  nonlinear <- stats::anova(linear, full, test = "Chisq")
  overall_column <- grep("^P|Pr", names(overall), value = TRUE)
  nonlinear_column <- grep("^P|Pr", names(nonlinear), value = TRUE)
  overall_probability <- if (length(overall_column)) as.numeric(overall[2, overall_column[1]]) else NA_real_
  nonlinear_probability <- if (length(nonlinear_column)) as.numeric(nonlinear[2, nonlinear_column[1]]) else NA_real_
  spline_rows[[length(spline_rows) + 1L]] <- data.frame(endpoint = endpoint, exposure = parent, knot_10 = knots[1], knot_50 = knots[2], knot_90 = knots[3], p_overall = overall_probability, p_nonlinear = nonlinear_probability)
}
save_csv(do.call(rbind, spline_rows), "spline_tests.csv")

sensitivity_rows <- list()
for (endpoint in endpoints) {
  early_removed <- eligible[eligible[[paste0(endpoint, "_event")]] == 0L | eligible[[paste0(endpoint, "_years")]] > 2, , drop = FALSE]
  medication_removed <- droplevels(eligible[as.character(eligible$medication_use) == "0", , drop = FALSE])
  for (label in c("Exclude_early_events", "Exclude_baseline_medication")) {
    section <- if (label == "Exclude_early_events") early_removed else medication_removed
    covariates <- if (label == "Exclude_early_events") c(social_covariates, clinical_covariates) else c(social_covariates, setdiff(clinical_covariates, "medication_use"))
    for (parent in exposures) {
      exposure <- paste0(parent, "_per_sd")
      fitted <- fit_cox(section, endpoint, exposure, covariates)
      result <- extract_fit(fitted, TRUE)
      result <- result[result$term == exposure, , drop = FALSE]
      result$endpoint <- endpoint
      result$exposure <- exposure
      result$analysis <- label
      sensitivity_rows[[length(sensitivity_rows) + 1L]] <- result
    }
  }
}
save_csv(do.call(rbind, sensitivity_rows), "sensitivity_models.csv")

competing_data <- function(data, endpoint) {
  event_date <- data[[paste0(endpoint, "_date")]]
  event_first <- !is.na(event_date) & event_date < data$death_date & event_date <= administrative_end
  event_first[is.na(event_first)] <- !is.na(event_date[is.na(event_first)]) & event_date[is.na(event_first)] <= administrative_end
  death_first <- !is.na(data$death_date) & data$death_date <= administrative_end & (is.na(event_date) | data$death_date <= event_date)
  status <- ifelse(event_first, 1L, ifelse(death_first, 2L, 0L))
  finish <- ifelse(status == 1L, as.numeric(event_date), ifelse(status == 2L, as.numeric(data$death_date), as.numeric(administrative_end)))
  data$competing_status <- status
  data$competing_years <- (finish - as.numeric(data$baseline_date)) / 365
  data[data$competing_years > 0 & is.finite(data$competing_years), , drop = FALSE]
}

fine_gray_rows <- list()
for (endpoint in endpoints) {
  analysis_data <- competing_data(eligible, endpoint)
  for (parent in exposures) for (exposure in paste0(parent, c("_per_sd", "_quartile"))) {
    variables <- unique(c("competing_years", "competing_status", exposure, social_covariates, clinical_covariates))
    section <- droplevels(analysis_data[stats::complete.cases(analysis_data[variables]), variables, drop = FALSE])
    design <- stats::model.matrix(stats::reformulate(c(exposure, social_covariates, clinical_covariates)), section)[, -1, drop = FALSE]
    fitted <- cmprsk::crr(section$competing_years, section$competing_status, cov1 = design, failcode = 1, cencode = 0)
    coefficients <- fitted$coef
    errors <- sqrt(diag(fitted$var))
    result <- data.frame(term = names(coefficients), estimate = exp(coefficients), lower_95 = exp(coefficients - 1.96 * errors), upper_95 = exp(coefficients + 1.96 * errors), p_value = 2 * stats::pnorm(-abs(coefficients / errors)), sample_size = nrow(section))
    result <- result[startsWith(result$term, exposure), , drop = FALSE]
    result$endpoint <- endpoint
    result$exposure <- exposure
    fine_gray_rows[[length(fine_gray_rows) + 1L]] <- result
  }
  quartile <- "bcaa_quartile"
  cumulative <- cmprsk::cuminc(analysis_data$competing_years, analysis_data$competing_status, group = analysis_data[[quartile]], cencode = 0)
  grDevices::pdf(file.path(options_run$output, paste0("cumulative_incidence_", endpoint, ".pdf")), width = 7, height = 6)
  plot(cumulative, xlab = "Years", ylab = "Cumulative incidence", main = endpoint)
  grDevices::dev.off()
  times <- seq(0, floor(max(analysis_data$competing_years)), by = 3)
  risk <- expand.grid(year = times, group = levels(analysis_data[[quartile]]), stringsAsFactors = FALSE)
  risk$at_risk <- mapply(function(time, group) sum(analysis_data$competing_years >= time & analysis_data[[quartile]] == group, na.rm = TRUE), risk$year, risk$group)
  save_csv(risk, paste0("risk_table_", endpoint, ".csv"))
}
save_csv(do.call(rbind, fine_gray_rows), "fine_gray_models.csv")

imaging_required <- c("eid", "imaging_date", "liver_fat", "asat", "gynoid_fat", "android_fat", "trunk_fat", "visceral_fat", "leg_fat", "arm_fat")
assert_fields(imaging, imaging_required, "imaging")
fat_data <- merge_one(baseline, imaging, "imaging")
fat_data$imaging_date <- parse_date(fat_data$imaging_date, "imaging_date")
fat_data$baseline_date <- parse_date(fat_data$baseline_date, "baseline_date")
fat_data$imaging_years <- as.numeric(fat_data$imaging_date - fat_data$baseline_date) / 365
fat_fields <- c("liver_fat", "asat", "gynoid_fat", "android_fat", "trunk_fat", "visceral_fat", "leg_fat", "arm_fat")
fat_required <- unique(c(exposures, social_covariates, clinical_covariates, "bmi", "imaging_years", fat_fields))
fat_data$egfr <- estimate_egfr(fat_data$creatinine, fat_data$sex, suppressWarnings(as.numeric(fat_data$age)))
fat_data <- fat_data[stats::complete.cases(fat_data[fat_required]), , drop = FALSE]
positive <- Reduce(`&`, lapply(fat_data[fat_fields], function(z) is.finite(z) & z > 0))
fat_data <- fat_data[positive & fat_data$imaging_years > 0, , drop = FALSE]
for (field in setdiff(fat_fields, c("liver_fat", "asat"))) fat_data[[field]] <- fat_data[[field]] / 1000
for (field in setdiff(fat_fields, "liver_fat")) fat_data[[field]] <- log(fat_data[[field]])
fat_data <- make_exposures(convert_covariates(fat_data), exposures)

linear_rows <- list()
fat_adjustment <- c(social_covariates, "bmi", "imaging_years", clinical_covariates)
for (outcome in fat_fields) for (parent in exposures) for (exposure in paste0(parent, c("_per_sd", "_quartile"))) for (model_name in names(model_sets)) {
  covariates <- if (model_name == "Model_1") character() else if (model_name == "Model_2") c(social_covariates, "bmi", "imaging_years") else fat_adjustment
  fitted <- stats::lm(stats::reformulate(c(exposure, covariates), outcome), data = fat_data)
  result <- extract_fit(fitted, FALSE)
  result <- result[startsWith(result$term, exposure), , drop = FALSE]
  result$outcome <- outcome
  result$exposure <- exposure
  result$model <- model_name
  linear_rows[[length(linear_rows) + 1L]] <- result
}
save_csv(do.call(rbind, linear_rows), "adiposity_models.csv")

fat_subgroup_rows <- list()
for (outcome in fat_fields) for (parent in exposures) for (stratum in subgroup_fields) {
  exposure <- paste0(parent, "_per_sd")
  adjusted <- setdiff(fat_adjustment, c(stratum, if (stratum == "age_group") "age" else character(), if (stratum == "healthy_diet") "diet_score" else character()))
  variables <- unique(c(outcome, exposure, stratum, adjusted))
  analysis_data <- droplevels(fat_data[stats::complete.cases(fat_data[variables]), variables, drop = FALSE])
  additive_formula <- stats::reformulate(c(exposure, stratum, adjusted), outcome)
  interaction_formula <- stats::as.formula(paste0(outcome, " ~ ", exposure, " * ", stratum, " + ", paste(adjusted, collapse = " + ")))
  p_interaction <- tryCatch({
    comparison <- stats::anova(stats::lm(additive_formula, analysis_data), stats::lm(interaction_formula, analysis_data))
    as.numeric(comparison[2, "Pr(>F)"])
  }, error = function(e) NA_real_)
  for (level in levels(factor(analysis_data[[stratum]]))) {
    section <- droplevels(analysis_data[as.character(analysis_data[[stratum]]) == level, , drop = FALSE])
    usable <- adjusted[vapply(section[adjusted], function(z) length(unique(z)) > 1L, logical(1L))]
    fitted <- tryCatch(stats::lm(stats::reformulate(c(exposure, usable), outcome), section), error = function(e) NULL)
    if (is.null(fitted)) next
    result <- extract_fit(fitted, FALSE)
    result <- result[result$term == exposure, , drop = FALSE]
    if (!nrow(result)) next
    result$outcome <- outcome
    result$exposure <- exposure
    result$stratum <- stratum
    result$level <- level
    result$p_interaction <- p_interaction
    fat_subgroup_rows[[length(fat_subgroup_rows) + 1L]] <- result
  }
}
fat_subgroup_results <- if (length(fat_subgroup_rows)) do.call(rbind, fat_subgroup_rows) else data.frame()
save_csv(fat_subgroup_results, "adiposity_subgroup_models.csv")

if (nrow(fat_subgroup_results)) {
  fat_plot_data <- fat_subgroup_results[is.finite(fat_subgroup_results$estimate) & is.finite(fat_subgroup_results$lower_95) & is.finite(fat_subgroup_results$upper_95), , drop = FALSE]
  fat_plot_data$label <- factor(paste(fat_plot_data$outcome, fat_plot_data$exposure, fat_plot_data$stratum, fat_plot_data$level, sep = " | "), levels = rev(paste(fat_plot_data$outcome, fat_plot_data$exposure, fat_plot_data$stratum, fat_plot_data$level, sep = " | ")))
  fat_forest <- ggplot2::ggplot(fat_plot_data, ggplot2::aes(estimate, label)) + ggplot2::geom_vline(xintercept = 0, linetype = 2) + ggplot2::geom_errorbar(ggplot2::aes(xmin = lower_95, xmax = upper_95), width = 0, orientation = "y") + ggplot2::geom_point() + ggplot2::labs(x = "Regression coefficient", y = NULL) + ggplot2::theme_bw()
  ggplot2::ggsave(file.path(options_run$output, "adiposity_subgroup_forest.pdf"), fat_forest, width = 10, height = max(7, nrow(fat_plot_data) * 0.16), limitsize = FALSE)
}

fat_spline_rows <- list()
for (outcome in fat_fields) for (parent in exposures) {
  knots <- stats::quantile(fat_data[[parent]], c(0.1, 0.5, 0.9), na.rm = TRUE, type = 7)
  if (anyDuplicated(knots)) stop("Nonunique spline knots for ", parent, ".", call. = FALSE)
  basis <- splines::ns(fat_data[[parent]], knots = knots[2], Boundary.knots = knots[c(1, 3)])
  spline_names <- paste0("fat_spline_", seq_len(ncol(basis)))
  analysis_data <- fat_data
  analysis_data[spline_names] <- basis
  full <- stats::lm(stats::reformulate(c(spline_names, fat_adjustment), outcome), analysis_data)
  linear <- stats::lm(stats::reformulate(c(parent, fat_adjustment), outcome), analysis_data)
  null <- stats::lm(stats::reformulate(fat_adjustment, outcome), analysis_data)
  overall <- stats::anova(null, full)
  nonlinear <- stats::anova(linear, full)
  fat_spline_rows[[length(fat_spline_rows) + 1L]] <- data.frame(outcome = outcome, exposure = parent, knot_10 = knots[1], knot_50 = knots[2], knot_90 = knots[3], p_overall = overall[2, "Pr(>F)"], p_nonlinear = nonlinear[2, "Pr(>F)"])
}
save_csv(do.call(rbind, fat_spline_rows), "adiposity_spline_tests.csv")

parameters <- data.frame(parameter = c("administrative_end", "ihd_codes", "diabetes_code", "liver_codes", "cox_ties", "days_per_year", "non_liver_fat_scale"), value = c("2023-03-31", "I20-I25", "E11", "K76.0,K75.8", "efron", "365", "natural_log"))
save_csv(parameters, "analysis_parameters.csv")
writeLines(capture.output(utils::sessionInfo()), file.path(options_run$output, "session_info.txt"))
writeLines("Completed", file.path(options_run$output, "SUCCESS"))
