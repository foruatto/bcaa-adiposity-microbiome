argv <- commandArgs(trailingOnly = TRUE)

usage_text <- paste(
  "Usage:",
  "  Rscript microbiome_analysis.R --metadata <sampleMetadata.csv> --output-dir <directory>",
  "",
  "The metadata file must contain sample_id, study_name, subject_id, body_site, disease, BMI, country, non_westernized, number_reads, antibiotics_current_use, pregnant, gender, age, age_category, and days_from_first_collection.",
  "The workflow retrieves public relative abundance, pathway abundance, and gene-family resources through curatedMetagenomicData.",
  sep = "\n"
)

read_options <- function(x) {
  if (length(x) == 1L && x == "--help") {
    cat(usage_text, "\n")
    quit(save = "no", status = 0L)
  }
  keys <- c("--metadata", "--output-dir")
  if (length(x) != 4L || !all(keys %in% x[c(1L, 3L)])) stop(usage_text, call. = FALSE)
  values <- setNames(x[c(2L, 4L)], x[c(1L, 3L)])
  list(metadata = values[["--metadata"]], output = values[["--output-dir"]])
}

options_run <- read_options(argv)

required_packages <- c("curatedMetagenomicData", "SummarizedExperiment", "Maaslin2", "dplyr", "tidyr", "tibble", "ggplot2", "metafor")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1L), quietly = TRUE)]
if (length(missing_packages)) stop("Missing packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)

if (!file.exists(options_run$metadata)) stop("Metadata file does not exist.", call. = FALSE)
if (dir.exists(options_run$output) && length(list.files(options_run$output, all.files = TRUE, no.. = TRUE))) stop("Output directory must be empty.", call. = FALSE)
dir.create(options_run$output, recursive = TRUE, showWarnings = FALSE)

save_csv <- function(x, name) utils::write.csv(x, file.path(options_run$output, name), row.names = FALSE, na = "")

assert_fields <- function(x, fields, label) {
  absent <- setdiff(fields, names(x))
  if (length(absent)) stop(label, " is missing: ", paste(absent, collapse = ", "), call. = FALSE)
}

load_metadata <- function(path) {
  separator <- if (grepl("[.]tsv$", path, ignore.case = TRUE)) "\t" else ","
  utils::read.table(path, header = TRUE, sep = separator, quote = "\"", comment.char = "", check.names = FALSE, stringsAsFactors = FALSE, na.strings = c("", "NA"))
}

metadata_fields <- c("sample_id", "study_name", "subject_id", "body_site", "disease", "BMI", "country", "non_westernized", "number_reads", "antibiotics_current_use", "pregnant", "gender", "age", "age_category", "days_from_first_collection")
metadata <- load_metadata(options_run$metadata)
assert_fields(metadata, metadata_fields, "metadata")
metadata$sample_id <- as.character(metadata$sample_id)
metadata$study_name <- as.character(metadata$study_name)
metadata$subject_id <- as.character(metadata$subject_id)
if (anyNA(metadata$sample_id) || any(!nzchar(metadata$sample_id)) || anyDuplicated(metadata$sample_id)) stop("sample_id must be unique and nonmissing.", call. = FALSE)

metadata$BMI <- suppressWarnings(as.numeric(metadata$BMI))
metadata$age <- suppressWarnings(as.numeric(metadata$age))
metadata$number_reads <- suppressWarnings(as.numeric(metadata$number_reads))
metadata$days_from_first_collection <- suppressWarnings(as.numeric(metadata$days_from_first_collection))

eligible <- metadata[
  (is.na(metadata$age_category) | metadata$age_category == "adult") &
  metadata$body_site == "stool" &
  metadata$disease == "healthy" &
  is.finite(metadata$BMI) &
  !is.na(metadata$country) &
  !is.na(metadata$non_westernized) &
  is.finite(metadata$number_reads) & metadata$number_reads >= 5000000 &
  (is.na(metadata$antibiotics_current_use) | metadata$antibiotics_current_use != "yes") &
  (is.na(metadata$pregnant) | metadata$pregnant != "yes"),
  , drop = FALSE
]

eligible$collection_order <- ifelse(is.finite(eligible$days_from_first_collection), eligible$days_from_first_collection, 0)
eligible <- eligible[order(eligible$study_name, eligible$subject_id, eligible$collection_order, eligible$sample_id), , drop = FALSE]
subject_key <- ifelse(is.na(eligible$subject_id) | !nzchar(eligible$subject_id), paste0("missing_subject_", eligible$sample_id), paste(eligible$study_name, eligible$subject_id, sep = "::"))
eligible <- eligible[!duplicated(subject_key), , drop = FALSE]
if (nrow(eligible) < 20L) stop("Fewer than 20 eligible samples.", call. = FALSE)

get_assay <- function(sample_table, assay_name) {
  result <- curatedMetagenomicData::returnSamples(sample_table, assay_name, rownames = "short")
  if (!inherits(result, "SummarizedExperiment")) stop("Resource retrieval did not return a SummarizedExperiment.", call. = FALSE)
  if (!assay_name %in% SummarizedExperiment::assayNames(result)) stop("Requested assay is unavailable: ", assay_name, call. = FALSE)
  matrix <- as.matrix(SummarizedExperiment::assay(result, assay_name))
  storage.mode(matrix) <- "numeric"
  shared <- intersect(sample_table$sample_id, colnames(matrix))
  if (length(shared) < 20L) stop("Fewer than 20 samples were returned for ", assay_name, ".", call. = FALSE)
  list(metadata = sample_table[match(shared, sample_table$sample_id), , drop = FALSE], matrix = matrix[, shared, drop = FALSE])
}

abundance_resource <- get_assay(eligible, "relative_abundance")
analysis_metadata <- abundance_resource$metadata
abundance_matrix <- abundance_resource$matrix

pc_candidates <- grep("(^|[_. ])Prevotella([_. ])copri($|[_. ])", rownames(abundance_matrix), ignore.case = TRUE, value = TRUE)
if (!length(pc_candidates)) pc_candidates <- grep("Prevotella.*copri", rownames(abundance_matrix), ignore.case = TRUE, value = TRUE)
if (!length(pc_candidates)) stop("Prevotella copri was not found in relative abundance resources.", call. = FALSE)
pc_abundance <- colSums(abundance_matrix[pc_candidates, , drop = FALSE], na.rm = TRUE)

analysis_metadata$pc_abundance <- pc_abundance[match(analysis_metadata$sample_id, names(pc_abundance))]
analysis_metadata$pc_detected <- as.integer(analysis_metadata$pc_abundance > 0.000001)
analysis_metadata$log_pc <- log1p(pmax(analysis_metadata$pc_abundance, 0))
analysis_metadata$overweight <- factor(ifelse(analysis_metadata$BMI >= 25, "Yes", "No"), levels = c("No", "Yes"))
analysis_metadata$obesity <- factor(ifelse(analysis_metadata$BMI >= 30, "Yes", "No"), levels = c("No", "Yes"))
analysis_metadata$lifestyle <- factor(ifelse(analysis_metadata$non_westernized == "no", "Westernized", ifelse(analysis_metadata$non_westernized == "yes", "Non-westernized", NA)), levels = c("Westernized", "Non-westernized"))
analysis_metadata$gender <- factor(analysis_metadata$gender)

westernized_metadata <- analysis_metadata[analysis_metadata$lifestyle == "Westernized", , drop = FALSE]
study_counts <- aggregate(cbind(overweight_yes = as.integer(westernized_metadata$overweight == "Yes"), overweight_no = as.integer(westernized_metadata$overweight == "No"), pc_yes = westernized_metadata$pc_detected, pc_no = 1L - westernized_metadata$pc_detected) ~ study_name, westernized_metadata, sum)
study_counts$total <- as.integer(table(westernized_metadata$study_name)[study_counts$study_name])
accepted_studies <- study_counts$study_name[study_counts$overweight_yes >= 15 & study_counts$overweight_no >= 15 & study_counts$pc_yes >= 15 & study_counts$pc_no >= 15 & study_counts$total >= 40]
analysis_metadata <- westernized_metadata[westernized_metadata$study_name %in% accepted_studies, , drop = FALSE]
analysis_metadata <- droplevels(analysis_metadata)
if (nrow(analysis_metadata) < 40L || length(unique(analysis_metadata$study_name)) < 2L) stop("Fewer than two studies satisfy the prespecified group-size criteria.", call. = FALSE)

flow <- data.frame(step = c("Metadata input", "Eligible healthy stool samples", "Samples with abundance resources", "Final analytical samples"), samples = c(nrow(metadata), nrow(eligible), ncol(abundance_matrix), nrow(analysis_metadata)))
save_csv(flow, "sample_flow.csv")
save_csv(study_counts, "study_counts.csv")

model_row <- function(data, outcome, exposure, covariates, family, study) {
  variables <- unique(c(outcome, exposure, covariates))
  section <- droplevels(data[stats::complete.cases(data[variables]), variables, drop = FALSE])
  if (nrow(section) < 20L || length(unique(section[[outcome]])) < 2L || stats::sd(section[[exposure]]) == 0) return(NULL)
  fitted <- tryCatch(stats::glm(stats::reformulate(c(exposure, covariates), outcome), data = section, family = family), error = function(e) NULL)
  if (is.null(fitted) || !isTRUE(fitted$converged)) return(NULL)
  coefficient <- stats::coef(fitted)[exposure]
  error <- sqrt(diag(stats::vcov(fitted)))[exposure]
  if (!is.finite(coefficient) || !is.finite(error) || error <= 0) return(NULL)
  data.frame(study = study, exposure = exposure, outcome = outcome, coefficient = coefficient, standard_error = error, estimate = exp(coefficient), lower_95 = exp(coefficient - 1.96 * error), upper_95 = exp(coefficient + 1.96 * error), p_value = 2 * stats::pnorm(-abs(coefficient / error)), sample_size = nrow(section))
}

meta_analyze <- function(rows, label) {
  if (is.null(rows) || !nrow(rows)) return(data.frame())
  groups <- split(rows, interaction(rows$exposure, rows$outcome, drop = TRUE))
  pooled <- lapply(groups, function(section) {
    if (nrow(section) < 2L) return(NULL)
    fitted <- metafor::rma.uni(yi = section$coefficient, sei = section$standard_error, method = "DL")
    data.frame(analysis = label, exposure = section$exposure[1], outcome = section$outcome[1], studies = nrow(section), estimate = exp(as.numeric(fitted$b)), lower_95 = exp(fitted$ci.lb), upper_95 = exp(fitted$ci.ub), p_value = fitted$pval, heterogeneity_i2 = fitted$I2)
  })
  pooled <- Filter(Negate(is.null), pooled)
  if (length(pooled)) do.call(rbind, pooled) else data.frame()
}

study_models <- list()
for (study in unique(analysis_metadata$study_name)) {
  section <- analysis_metadata[analysis_metadata$study_name == study, , drop = FALSE]
  specifications <- list(
    c("overweight", "log_pc"),
    c("overweight", "pc_detected")
  )
  for (specification in specifications) {
    outcome <- specification[1]
    exposure <- specification[2]
    for (adjustment in list(character(), "age", "gender", c("age", "gender"))) {
      row <- model_row(section, outcome, exposure, adjustment, stats::binomial(), study)
      if (!is.null(row)) {
        row$adjustment <- if (length(adjustment)) paste(adjustment, collapse = "+") else "unadjusted"
        study_models[[length(study_models) + 1L]] <- row
      }
    }
  }
}
study_model_results <- if (length(study_models)) do.call(rbind, study_models) else data.frame()
save_csv(study_model_results, "study_logistic_models.csv")
if (nrow(study_model_results)) {
  meta_rows <- do.call(rbind, lapply(split(study_model_results, study_model_results$adjustment), function(section) meta_analyze(section, section$adjustment[1])))
  save_csv(meta_rows, "random_effects_meta_analysis.csv")
}

wilson_interval <- function(events, total, confidence = 0.95) {
  z <- stats::qnorm(1 - (1 - confidence) / 2)
  proportion <- events / total
  denominator <- 1 + z^2 / total
  center <- (proportion + z^2 / (2 * total)) / denominator
  half_width <- z * sqrt(proportion * (1 - proportion) / total + z^2 / (4 * total^2)) / denominator
  c(proportion = proportion, lower = max(0, center - half_width), upper = min(1, center + half_width))
}

detection_rows <- list()
for (study in unique(analysis_metadata$study_name)) for (group in levels(analysis_metadata$overweight)) {
  section <- analysis_metadata[analysis_metadata$study_name == study & analysis_metadata$overweight == group, , drop = FALSE]
  interval <- wilson_interval(sum(section$pc_detected), nrow(section))
  detection_rows[[length(detection_rows) + 1L]] <- data.frame(study = study, overweight = group, detected = sum(section$pc_detected), total = nrow(section), proportion = interval[1], lower_95 = interval[2], upper_95 = interval[3])
}
save_csv(do.call(rbind, detection_rows), "pc_detection_summary.csv")

comparison_rows <- list()
for (study in unique(analysis_metadata$study_name)) {
  section <- analysis_metadata[analysis_metadata$study_name == study, , drop = FALSE]
  contingency <- table(section$overweight, section$pc_detected)
  detection_test <- if (all(dim(contingency) == c(2, 2))) stats::fisher.test(contingency) else NULL
  positive <- section[section$pc_detected == 1L, , drop = FALSE]
  abundance_test <- if (length(unique(positive$overweight)) == 2L) stats::wilcox.test(log_pc ~ overweight, positive, exact = FALSE) else NULL
  comparison_rows[[length(comparison_rows) + 1L]] <- data.frame(study = study, detection_p = if (is.null(detection_test)) NA_real_ else detection_test$p.value, positive_abundance_p = if (is.null(abundance_test)) NA_real_ else abundance_test$p.value)
}
save_csv(do.call(rbind, comparison_rows), "pc_group_comparisons.csv")

pathway_resource <- get_assay(analysis_metadata, "pathway_abundance")
pathway_metadata <- pathway_resource$metadata
pathway_matrix <- pathway_resource$matrix
unstratified <- pathway_matrix[!grepl("\\|", rownames(pathway_matrix)), , drop = FALSE]
stratified <- pathway_matrix[grepl("\\|", rownames(pathway_matrix)) & !grepl("^(UNMAPPED|UNINTEGRATED)\\|", rownames(pathway_matrix)), , drop = FALSE]
target_pathways <- grep("branched.*amino.*acid.*biosynthesis|valine.*biosynthesis|leucine.*biosynthesis|isoleucine.*biosynthesis", rownames(unstratified), ignore.case = TRUE, value = TRUE)
if (!length(target_pathways)) stop("No BCAA biosynthesis pathways were found.", call. = FALSE)

correlation_one <- function(x, y, label, group = "All") {
  valid <- is.finite(x) & is.finite(y)
  if (sum(valid) < 10L || stats::sd(x[valid]) == 0 || stats::sd(y[valid]) == 0) return(data.frame(pathway = label, group = group, rho = NA_real_, p_value = NA_real_, sample_size = sum(valid)))
  tested <- suppressWarnings(stats::cor.test(x[valid], y[valid], method = "spearman", exact = FALSE))
  data.frame(pathway = label, group = group, rho = unname(tested$estimate), p_value = tested$p.value, sample_size = sum(valid))
}

pc_positive <- pathway_metadata$pc_detected == 1L
correlation_rows <- list()
for (pathway in target_pathways) {
  correlation_rows[[length(correlation_rows) + 1L]] <- correlation_one(pathway_metadata$log_pc[pc_positive], pathway_matrix[pathway, pc_positive], pathway, "All")
  for (gender in levels(pathway_metadata$gender)) {
    selected <- pc_positive & pathway_metadata$gender == gender
    correlation_rows[[length(correlation_rows) + 1L]] <- correlation_one(pathway_metadata$log_pc[selected], pathway_matrix[pathway, selected], pathway, as.character(gender))
  }
}
correlations <- do.call(rbind, correlation_rows)
correlations$q_value <- stats::p.adjust(correlations$p_value, method = "BH")
save_csv(correlations, "bcaa_pathway_correlations.csv")

split_feature <- strsplit(rownames(stratified), "\\|")
left <- vapply(split_feature, `[`, character(1L), 1L)
right <- vapply(split_feature, `[`, character(1L), 2L)
left_taxon <- grepl("g__|s__", left)
right_taxon <- grepl("g__|s__", right)
pathway_name <- ifelse(left_taxon & !right_taxon, right, ifelse(right_taxon & !left_taxon, left, NA_character_))
taxon_name <- ifelse(left_taxon & !right_taxon, left, ifelse(right_taxon & !left_taxon, right, NA_character_))
usable <- !is.na(pathway_name) & !is.na(taxon_name) & pathway_name %in% target_pathways
stratified <- stratified[usable, , drop = FALSE]
pathway_name <- pathway_name[usable]
taxon_name <- taxon_name[usable]

contribution_rows <- list()
for (pathway in unique(pathway_name)) {
  index <- which(pathway_name == pathway)
  total <- unstratified[pathway, ]
  stratified_total <- colSums(stratified[index, , drop = FALSE], na.rm = TRUE)
  for (taxon in unique(taxon_name[index])) {
    taxon_values <- colSums(stratified[index[taxon_name[index] == taxon], , drop = FALSE], na.rm = TRUE)
    contribution_rows[[length(contribution_rows) + 1L]] <- data.frame(pathway = pathway, taxon = taxon, summed_abundance = sum(taxon_values), mean_percent_of_unstratified = mean(100 * taxon_values / pmax(total, 1e-12)), mean_percent_of_stratified = mean(100 * taxon_values / pmax(stratified_total, 1e-12)), detected_samples = sum(taxon_values > 0))
  }
}
contribution_results <- do.call(rbind, contribution_rows)
contribution_results$taxon_group <- ifelse(grepl("Prevotella.*copri", contribution_results$taxon, ignore.case = TRUE), "Prevotella_copri", contribution_results$taxon)
save_csv(contribution_results, "pathway_contributions_by_taxon.csv")
merged_contributions <- aggregate(cbind(summed_abundance, mean_percent_of_unstratified, mean_percent_of_stratified, detected_samples) ~ pathway + taxon_group, contribution_results, sum)
save_csv(merged_contributions, "pathway_contributions_merged.csv")

target_species <- c("Bacteroides uniformis", "Gemmiger formicilis", "Prevotella copri", "Anaerostipes hadrus", "Blautia obeum", "Eubacterium rectale", "Blautia wexlerae", "Ruminococcus torques", "Bifidobacterium adolescentis")
species_rows <- list()
for (species in target_species) {
  matches <- grep(paste0("^", gsub(" ", "[ _.]", species), "$"), rownames(abundance_matrix), ignore.case = TRUE, value = TRUE)
  if (!length(matches)) next
  values <- colSums(abundance_matrix[matches, analysis_metadata$sample_id, drop = FALSE], na.rm = TRUE)
  for (group in levels(analysis_metadata$overweight)) {
    selected <- analysis_metadata$overweight == group
    positive_values <- values[selected & values > 0]
    species_rows[[length(species_rows) + 1L]] <- data.frame(species = species, overweight = group, samples = sum(selected), detected = sum(values[selected] > 0), detection_rate = mean(values[selected] > 0), median_positive_abundance = if (length(positive_values)) stats::median(positive_values) else NA_real_)
  }
}
species_summary <- if (length(species_rows)) do.call(rbind, species_rows) else data.frame()
save_csv(species_summary, "target_species_summary.csv")

gene_resource <- get_assay(analysis_metadata, "gene_families")
gene_matrix <- gene_resource$matrix
gene_matches <- grep("branched|ilvA|ilvB|ilvC|ilvD|ilvE|leuA|leuB|leuC|leuD|valine|leucine|isoleucine", rownames(gene_matrix), ignore.case = TRUE, value = TRUE)
gene_rows <- if (length(gene_matches)) data.frame(gene_family = gene_matches, mean_abundance = rowMeans(gene_matrix[gene_matches, , drop = FALSE], na.rm = TRUE), prevalence = rowMeans(gene_matrix[gene_matches, , drop = FALSE] > 0, na.rm = TRUE)) else data.frame(gene_family = character(), mean_abundance = numeric(), prevalence = numeric())
save_csv(gene_rows, "bcaa_gene_family_summary.csv")

align_maaslin <- function(matrix, metadata) {
  shared <- intersect(metadata$sample_id, colnames(matrix))
  features <- as.data.frame(t(matrix[, shared, drop = FALSE]), check.names = FALSE)
  rownames(features) <- shared
  aligned_metadata <- metadata[match(shared, metadata$sample_id), , drop = FALSE]
  rownames(aligned_metadata) <- shared
  list(features = features, metadata = aligned_metadata)
}

run_maaslin <- function(matrix, metadata, fixed_effects, random_effects, folder) {
  aligned <- align_maaslin(matrix, metadata)
  fixed_effects <- fixed_effects[fixed_effects %in% names(aligned$metadata)]
  random_effects <- random_effects[random_effects %in% names(aligned$metadata)]
  if (!length(fixed_effects)) stop("No valid fixed effects for MaAsLin2.", call. = FALSE)
  variable_features <- aligned$features[, vapply(aligned$features, function(z) stats::sd(z, na.rm = TRUE) > 0, logical(1L)), drop = FALSE]
  if (!ncol(variable_features)) stop("No variable features for MaAsLin2.", call. = FALSE)
  Maaslin2::Maaslin2(input_data = variable_features, input_metadata = aligned$metadata, output = file.path(options_run$output, folder), fixed_effects = fixed_effects, random_effects = random_effects, normalization = "NONE", transform = "LOG", analysis_method = "LM", standardize = FALSE, plot_scatter = FALSE, plot_heatmap = FALSE)
}

run_maaslin(unstratified, pathway_metadata, "overweight", "study_name", "maaslin_overweight")
run_maaslin(unstratified, pathway_metadata, "BMI", "study_name", "maaslin_continuous_bmi")
for (gender in levels(pathway_metadata$gender)) {
  selected <- !is.na(pathway_metadata$gender) & pathway_metadata$gender == gender
  if (sum(selected) >= 20L && length(unique(pathway_metadata$overweight[selected])) == 2L) run_maaslin(unstratified[, selected, drop = FALSE], droplevels(pathway_metadata[selected, , drop = FALSE]), "overweight", "study_name", paste0("maaslin_gender_", make.names(gender)))
}

parameters <- data.frame(parameter = c("minimum_reads", "pc_detection_threshold", "bmi_overweight_threshold", "bmi_obesity_threshold", "minimum_group_size", "minimum_study_size", "meta_method", "multiple_testing"), value = c("5000000", "0.000001", "25", "30", "15", "40", "DerSimonian-Laird", "Benjamini-Hochberg"))
save_csv(parameters, "analysis_parameters.csv")
writeLines(capture.output(utils::sessionInfo()), file.path(options_run$output, "session_info.txt"))
writeLines("Completed", file.path(options_run$output, "SUCCESS"))
