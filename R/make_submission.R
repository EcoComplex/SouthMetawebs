# Builds the double-anonymous main document (Main.Rmd) from Manuscript.Rmd
# for the Journal of Biogeography submission, and renders Main and Title.
# Run from the repository root: source("R/make_submission.R")
# Main.Rmd and Title.Rmd are not tracked by git (see .gitignore).

ms <- readLines("Manuscript.Rmd", encoding = "UTF-8")

yaml_end    <- which(ms == "---")[2]
abs_start   <- which(ms == "## Abstract")[1]
intro_start <- which(ms == "## Introduction")[1]
stopifnot(!is.na(yaml_end), !is.na(abs_start), !is.na(intro_start))

# Keywords: 6-10, in alphabetical order (journal requirement)
keywords <- paste0(
  "**Keywords:** Antarctica, ecosystem size, food web structure, ",
  "latitudinal gradient, macroecology, marine food webs, metaweb, ",
  "network complexity, Southwest Atlantic, trophic level"
)

# Drop the author and affiliation block (between YAML and Abstract)
main <- c(ms[1:yaml_end], "",
          ms[abs_start:(intro_start - 1)], keywords, "",
          ms[intro_start:length(ms)])

# Replace the identifying repository link (GitHub URL and Zenodo DOI)
repo_pattern <- "in a public GitHub repository \\(https://github\\.com/[^)]*\\), with archived code and data available on Zenodo \\([^)]*\\)"
repo_new <- "in a public repository (anonymized link for review: ANONYMOUS_REPO_URL)"
if (!any(grepl(repo_pattern, main)))
  warning("Repository sentence not found: check the code availability text in Main.Rmd")
main <- sub(repo_pattern, repo_new, main)

# Flag any remaining identifying text
id_pattern <- "github|EcoComplex|zenodo|CADIC|UNTDF|CESIMAR|IIMyC|Helsinki|lsaravia"
hits <- grep(id_pattern, main, ignore.case = TRUE)
if (length(hits) > 0)
  warning("Possible identifying text in Main.Rmd, lines: ", paste(hits, collapse = ", "))

writeLines(main, "Main.Rmd", useBytes = TRUE)

# Sync Title.Rmd with Manuscript.Rmd: title and author/affiliation block.
# Title.Rmd is edited by hand only for running title, corresponding author,
# acknowledgements, funding, etc.; authors and affiliations come from Manuscript.Rmd.
tt <- readLines("Title.Rmd", encoding = "UTF-8")
ms_title <- grep("^title:", ms[1:yaml_end], value = TRUE)[1]
tt[grep("^title:", tt)[1]] <- ms_title
run_line  <- grep("^\\*\\*Running title:", tt)[1]
corr_line <- grep("^\\*\\*Corresponding author:", tt)[1]
stopifnot(!is.na(run_line), !is.na(corr_line))
authors_block <- ms[(yaml_end + 1):(abs_start - 1)]
authors_block <- authors_block[cumsum(nzchar(trimws(authors_block))) > 0]  # drop leading blanks
tt <- c(tt[1:run_line], "", authors_block, tt[corr_line:length(tt)])
writeLines(tt, "Title.Rmd", useBytes = TRUE)

rmarkdown::render("Main.Rmd")
# The journal does not accept a PDF title page: render it as Word
rmarkdown::render("Title.Rmd", output_format = "word_document")
