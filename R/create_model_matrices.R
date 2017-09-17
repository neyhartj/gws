#' Create a random effect incidence matrix
#'
#' @description
#' Creates an incidence matrix for the random effects of a mixed-model
#'
#' @param random A \code{formula} specifying the the random effect terms of the
#' linear mixed-model. Must be preceded by a \code{~} and can include the modifier
#' verbs \code{\link[at]{sommer}}, \code{\link[us]{sommer}}, and \code{\link[diag]{base}}
#' to specify different structures. If a random effect term has an accompanying
#' variance-covariance matrix (see \code{vcov}), the term should be nested using
#' the function \code{\link[g]{sommer}}
#' @param data A \code{data.frame} that includes the response and random effect
#' predictors. If not supplied, the function will look in the current running
#' environment.
#' @param vcov A named list of variance-covariance matrices that correspond to
#' the random effect terms that are specified with the \code{\link[g]{sommer}}
#' function. If not passed, an identity matrix is supplied.
#'
#' @return
#' A nested list of model matrices and variance-covariance matrices. For each
#' random effect term, the incidence matrix \code{Z} is returned along with its
#' accompanying variance-covariance matrix \code{K}.
#'
#' @import dplyr
#' @import stringr
#' @importFrom sommer at g us and
#'
#' @export
#'
ranef_model_matrix <- function(random, data, vcov) {

  # Input check
  if (!inherits(random, "formula"))
    stop("The input 'random' must be a model formula.")

  # If data is missing, use the current environment
  if (missing(data)) {
    data <- environment(random)

  } else {
    # Convert data to a data.frame
    data <- as.data.frame(data) %>%
      mutate_if(is.character, as.factor)

  }

  # If vcov is supplied, perform some error checking
  if (!missing(vcov)) {
    if (is.null(names(vcov)))
      stop("vcov must be a named list.")

    # Make sure each is has column names
    if (any(lapply(vcov, dimnames) %>% sapply(is.null)))
      stop("The matrices in vcov must have dimnames.")

    # Extract the names
    vcov_names <- names(vcov)

  }

  # Edit the formula to remove the intercept
  random1 <- as.formula(paste("~ -1 +", random)[-1])

  # Create the base model matrix
  model_frame_base <- model.frame(random1, data = data)

  # Determine the variables in the formula
  rand_vars <- str_split(string = as.character(random)[-1], pattern = "\\+", simplify = FALSE) %>%
    unlist() %>%
    str_trim()

  # Create a list
  Z_list <- list()

  # Cycle through the random effect terms and create a model.matrix
  for (term in rand_vars) {

    ## Detect the modifier verbs in the formula (i.e. at, diag, us)
    is_at <- str_detect(string = term, "at\\(")
    is_us <- str_detect(string = term, "us\\(")
    is_diag <- str_detect(string = term, "diag\\(")

    # Are any terms nested in the 'g' function?
    is_g <- str_detect(string = term, "g\\(")

    # Create a formula for that specific term
    term_formula <- as.formula(paste("~ -1 +", term))

    # Pull out the model matrix
    model_matrix <- model.matrix(object = term_formula, data = model_frame_base)

    # Find the regressor that with the modifier
    if (is_at) {
      at_var <- str_extract(string = term, pattern = "at\\([A-Za-z0-9]*\\)") %>%
        str_replace_all(pattern = "at\\(|\\)", replacement = "")

      # What is the variable to the right of the parentheses/colon?
      sec_var <- str_split(string = term, pattern = ":") %>%
        unlist() %>%
        tail(-1) %>%
        paste(collapse = "|") %>%
        str_replace_all(pattern = "g\\(", "g\\\\(") %>%
        str_replace_all(pattern = "\\)", "\\\\)")

      # Pull out the levels of that variable
      at_var_levels <- levels(data[[at_var]])

      # For each level, subset the model_matrix
      for (lev in at_var_levels) {
        # Subset the matrix
        mat <- model_matrix[,str_detect(string = colnames(model_matrix), pattern = lev)]

        # Adjust the columnnames
        colnames(mat) <- colnames(mat) %>%
          str_replace_all(pattern = lev, replacement = "") %>%
          str_replace_all(pattern = at_var, replacement = lev) %>%
          str_replace_all(pattern = sec_var, replacement = "")

        # Determine a list name
        lev_list_name <- str_replace(string = term, pattern = at_var, replacement = lev)

        # If the term has a specified covariance matrix, figure out what term it is
        if (is_g) {
          # Is the 'g' at the beginning of the term?
          is_start_g <- str_detect(string = term, pattern = "^g")

          g_term <- str_extract(string = term, pattern = "g\\(.*\\)") %>%
            str_replace_all(pattern = "g\\(|\\)", replacement = "")

          # Is the name in the vcov list?
          if (!g_term %in% vcov_names) {
            stop("A term was specified using g(.), but the term is not present
                 in the 'vcov' list.")

          } else {
            # Pull out the vcov matrix
            vcov_term <- vcov[[g_term]]

          }

        } else {
          is_start_g <- FALSE

          vcov_term <- diag(ncol(mat)) %>%
            structure(dimnames = list(colnames(mat), colnames(mat)))

        }


        # Add to the list
        Z_list[[lev_list_name]] <- list(Z = mat, K = vcov_term)

      }

    } else {

      # If so, figure out what term it is
      if (is_g) {
        # Is the 'g' at the beginning of the term?
        is_start_g <- str_detect(string = term, pattern = "^g")

        # The line '.*' matches all strings.
        g_term <- str_extract(string = term, pattern = "g\\(.*\\)") %>%
          str_replace_all(pattern = "g\\(|\\)", replacement = "")

        # Is the name in the vcov list?
        if (!g_term %in% vcov_names) {
          stop("A term was specified using g(.), but the term is not present in the 'vcov' list.")

        } else {
          # Pull out the vcov matrix
          vcov_term <- vcov[[g_term]]

        }

      } else {
        is_start_g <- FALSE

        vcov_term <- diag(ncol(model_matrix)) %>%
          structure(dimnames = list(colnames(model_matrix), colnames(model_matrix)))

      }

      # Remove the variable name from the matrix column names
      if (is_start_g) {
        colnames(model_matrix) <- colnames(model_matrix) %>%
          str_replace_all(pattern = "g\\(.*\\)", replacement = "")

      } else {
        colnames(model_matrix) <- colnames(model_matrix) %>%
          str_replace_all(pattern = term, replacement = "")

      }

      # Return the matrix
      Z_list[[term]] <- list(Z = model_matrix, K = vcov_term)

    }

  }

  ## Need to revise this for "un" and "diag" verbs



  # Return the list
  return(Z_list)

} # Close the function



#' Create a fixed effect incidence matrix
#'
#' @description
#' Creates an incidence matrix for the fixed effects of a mixed-model. This is a
#' very thin wrapper around \code{\link[model.matrix]{stats}}.
#'
#' @param fixed A \code{formula} specifying the the response variable and the
#' fixed effect terms of the linear mixed-model. Must take the form \code{response
#' ~ predictor(s)}.
#' @param data A \code{data.frame} that includes the response and fixed effect
#' predictors. If not supplied, the function will look in the current running
#' environment.
#'
#' @return
#' A model matrix.
#'
#' @export
#'
fixef_model_matrix <- function(fixed, data) {

  model.matrix(fixed, data = data)

} # Close the function



#' Create a residual effect covariance matrix
#'
#' @description
#' Creates an covariance matrix for the random effects of a mixed-model
#'
#' @param resid A \code{formula} specifying the the structure of the residual
#' effects in the linear mixed-model. Must be preceded by a \code{~} and can
#' include the modifier verbs \code{\link[at]{sommer}}, \code{\link[us]{sommer}},
#' and \code{\link[diag]{base}} to specify different structures. The formulae
#' can take the same form as those in \code{\link[ranef_model_matrix]{gws}},
#' except the term must be "units".
#' @param data A \code{data.frame} that includes the response and predictors.
#' If not supplied, the function will look in the current running environment.
#'
#' @return
#' A nested list of model matrices.
#'
#' @import dplyr
#' @import stringr
#' @importFrom sommer at g us and
#'
#' @export
#'
resid_model_matrix <- function(resid, data) {

  # Input check
  if (!inherits(resid, "formula"))
    stop("The input 'resid' must be a model formula.")

  # If data is missing, use the current environment
  if (missing(data)) {
    data <- environment(resid)

  } else {
    # Convert data to a data.frame
    data <- as.data.frame(data) %>%
      mutate_if(is.character, as.factor)

  }

  # Detect "units" variable in the model and in the specific format
  detect_format <- str_detect(string = as.character(resid)[2],
                              pattern = "at\\(.*\\):units|units")

  if (!detect_format) {
    stop("The 'resid' formula must take the form 'resid = ~ units' or 'resid = ~ at(.):units'.")

  } else {
    # Edit the formula to remove the intercept and the "units"
    resid1 <- paste("~ -1 +", as.character(resid)[2]) %>%
      str_replace_all(pattern = ":units|\\+ units", replacement = "") %>%
      str_trim() %>%
      as.formula()

  }

  # Create the base model frame
  model_frame_base <- model.frame(resid1, data = data)

  # If the model_frame is dimensionless (i.e. resid = ~ units), return the R list
  if (ncol(model_frame_base) == 0) {
    return(list(units = diag(nrow(model_frame_base))))

  } else {
    # Otherwise create a list of R matrices
    resid_vars <- as.character(resid1)[2] %>%
      str_replace(pattern = "-1 \\+ ", replacement = "")

    # Create a list
    R_list <- list()

    # Cycle through the random effect terms and create a model.matrix
    for (term in resid_vars) {

      ## Detect the modifier verbs in the formula (i.e. at, diag, us)
      is_at <- str_detect(string = term, "at\\(")
      is_us <- str_detect(string = term, "us\\(")
      is_diag <- str_detect(string = term, "diag\\(")

      # Create a formula for that specific term
      term_formula <- as.formula(paste("~ -1 +", term))

      # Pull out the model matrix
      model_matrix <- model.matrix(object = term_formula, data = model_frame_base)

      # Find the regressor that with the modifier
      if (is_at) {
        at_var <- str_extract(string = term, pattern = "at\\([A-Za-z0-9]*\\)") %>%
          str_replace_all(pattern = "at\\(|\\)", replacement = "")

        # Pull out the levels of that variable
        at_var_levels <- levels(data[[at_var]])

        # For each level, subset the model_matrix
        for (lev in at_var_levels) {
          # Subset the matrix - this will become the diagonal of the matrix
          diag_mat <- model_matrix[,str_detect(string = colnames(model_matrix), pattern = lev)]

          # Create an empty matrix
          mat <- matrix(data = 0, nrow = length(diag_mat), ncol = length(diag_mat))
          diag(mat) <- diag_mat

          # Determine a list name
          lev_list_name <- str_replace(string = term, pattern = at_var, replacement = lev)

          # Add to the list
          R_list[[lev_list_name]] <- mat

        }
      }
    }
  }
  ## Need to revise this for "un" and "diag" verbs

  # Return the list
  return(R_list)

} # Close the function









