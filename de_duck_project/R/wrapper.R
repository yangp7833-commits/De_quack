#!/usr/bin/env Rscript

# de_duck_project R wrapper for the local Python de_duck library.
# This wrapper uses reticulate to expose core.py and viz.py functions from R.

if (!requireNamespace("reticulate", quietly = TRUE)) {
  stop("Please install the 'reticulate' package to use de_duck from R.")
}
if (!requireNamespace("R6", quietly = TRUE)) {
  stop("Please install the 'R6' package to use the de_duck wrapper.")
}

resolve_de_duck_path <- function(python_path = NULL) {
  # 1. If the user explicitly provided a path, use it.
  if (!is.null(python_path)) {
    return(normalizePath(python_path, mustWork = TRUE))
  }
  
  # 2. Check if the library is installed in the current Python environment via pip
  if (reticulate::py_module_available("de_duck")) {
    py_env <- reticulate::import("de_duck")
    # This automatically finds the real pip installation directory!
    pip_path <- dirname(py_env$`__file__`)
    return(normalizePath(pip_path, mustWork = TRUE))
  }
  
  # 3. Fallback for local development or uninstalled source code
  local_candidates <- c(
    file.path(getwd(), "de_duck"),
    file.path(getwd(), "de_duck_project", "de_duck"),
    getwd()
  )
  
  for (candidate in local_candidates) {
    if (file.exists(file.path(candidate, "core.py"))) {
      return(normalizePath(candidate, mustWork = TRUE))
    }
  }
  
  stop("de_duck Python module not found. Please run 'pip install de_duck' or provide python_path.")
}

load_de_duck_modules <- function(python_path = NULL) {
  python_path <- resolve_de_duck_path(python_path)
  message("Loading de_duck Python modules from: ", python_path)
  core <- reticulate::import_from_path("core", path = python_path, convert = TRUE)
  viz <- reticulate::import_from_path("viz", path = python_path, convert = TRUE)
  list(core = core, viz = viz)
}

DeDuck <- R6::R6Class(
  "DeDuck",
  public = list(
    py_core = NULL,
    py_viz = NULL,
    py_db = NULL,
    python_path = NULL,
    db_path = NULL,

    initialize = function(python_path = NULL, db_path = "SQL.duckdb", use_condaenv = NULL, use_virtualenv = NULL) {
      self$python_path <- python_path
      self$db_path <- db_path
      if (!is.null(use_condaenv)) {
        reticulate::use_condaenv(use_condaenv, required = TRUE)
      }
      if (!is.null(use_virtualenv)) {
        reticulate::use_virtualenv(use_virtualenv, required = TRUE)
      }
      modules <- load_de_duck_modules(python_path)
      self$py_core <- modules$core
      self$py_viz <- modules$viz
      invisible(self)
    },

    connect = function(db_path = NULL) {
      target_path <- if (is.null(db_path)) self$db_path else db_path
      
      self$py_db <- self$py_core$de_duckling(db_path = target_path)
      self$py_db$connect()
      invisible(self)
    },

    close = function() {
      if (!is.null(self$py_db)) {
        self$py_db$close()
        self$py_db <- NULL
      }
      invisible(self)
    },

    insert_to_database = function(info, tool = NULL, date = NULL, file_path = NULL, experiment_name = NULL, comparison_label = NULL) {
      stopifnot(!is.null(self$py_db))
      self$py_db$insert_to_database(
        info,
        tool = if (is.null(tool)) NULL else tool,
        date = if (is.null(date)) NULL else date,
        file_path = if (is.null(file_path)) NULL else file_path,
        experiment_name = if (is.null(experiment_name)) NULL else experiment_name,
        comparison_label = if (is.null(comparison_label)) NULL else comparison_label
      )
      invisible(TRUE)
    },

    query = function(table, save_path = NULL, ...) {
      stopifnot(!is.null(self$py_db))
      filters <- list(...)
      args <- c(list(table = table, save_path = save_path), filters)
      result <- do.call(self$py_db$query, args)
      result
    },

    delete_gene_results = function(...) {
      stopifnot(!is.null(self$py_db))
      args <- list(...)
      do.call(self$py_db$delete_gene_results, args)
      invisible(TRUE)
    },

    delete_experiments = function(...) {
      stopifnot(!is.null(self$py_db))
      args <- list(...)
      do.call(self$py_db$delete_experiments, args)
      invisible(TRUE)
    },

    initialize_gene_table = function(species = "human") {
      stopifnot(!is.null(self$py_db))
      self$py_db$initialize_gene_table(species)
      invisible(TRUE)
    },

    execute_raw = function(query) {
      stopifnot(!is.null(self$py_db))
      self$py_db$execute_raw(query)
      invisible(TRUE)
    },

    volcano_plot = function(df, padj = 0.05, log2fc = 1, plot_file = NULL, pvalue = NULL, significant_up = "red", insignificant = "grey", significant_down = "blue", label = 5) {
      result <- self$py_viz$volcano_plot(
        df,
        padj = padj,
        log2fc = log2fc,
        plot_file = if (is.null(plot_file)) NULL else plot_file,
        pvalue = if (is.null(pvalue)) NULL else pvalue,
        significant_up = significant_up,
        insignificant = insignificant,
        significant_down = significant_down,
        label = label
      )
      result
    },
    
    # Garbage collection guard to unlock DuckDB automatically
    finalize = function() {
      self$close()
    }
    )
    )

# Convenience constructor for R users
de_duck <- function(python_path = NULL, db_path = "SQL.duckdb", use_condaenv = NULL, use_virtualenv = NULL) {
  DeDuck$new(
    python_path = python_path,
    db_path = db_path,
    use_condaenv = use_condaenv,
    use_virtualenv = use_virtualenv
  )
}