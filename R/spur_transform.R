#' Transform Variables Using Spatial Transformations
#'
#' This function applies spatial transformations to variables, similar to the Stata spurtransform command.
#' It supports different transformation methods: "lbmgls" (default), "iso", "cluster", and "nn".
#'
#' @param data A data frame containing the variables and spatial coordinates.
#' @param vars Character vector of variable names to transform.
#' @param prefix String prefix to add to transformed variable names.
#' @param transformation Type of transformation: "lbmgls" (default), "iso", "cluster", or "nn".
#' @param radius Numeric radius parameter for "iso" transformation.
#' @param clustvar Name of clustering variable for "cluster" transformation.
#' @param latlong Logical; if TRUE, coordinates are latitude/longitude.
#' @param lat Name of the latitude column (or y-coordinate if latlong=FALSE).
#' @param lon Name of the longitude column (or x-coordinate if latlong=FALSE).
#' @param replace Logical; if TRUE, replace existing variables with the same names.
#' @param separately Logical; if TRUE, handle missing values in each variable separately.
#'
#' @return A data frame with the original and transformed variables.
#' @export
#' 

# Stata equivalent:
# program spurtransform, sortpreserve 
#   version 14
#   syntax varlist(numeric) [if] [in] , PREfix(string) [ transformation(string) radius(real -1) clustvar(varname) latlong Replace separately] 
options(digits = 8)
spur_transform <- function(data, vars, prefix = "tr", transformation = "lbmgls", 
                          radius = NULL, clustvar = NULL, latlong = TRUE,
                          lat = "lat", lon = "lon", 
                          replace = FALSE, separately = FALSE) {
  
    # Check inputs - keep this validation logic
    
    # Stata equivalent:
    # if "`transformation'"=="" {
    #   local transformation "lbmgls"
    # }
    # Default is already handled in R function signature
    
    # Stata equivalent:
    # if "`transformation'"!="iso" & `radius'!=-1 {
    #   di as error "Option radius only allowed with transformation(iso)."
    #   exit
    # }
    if (transformation == "iso" && is.null(radius)) {
        stop("Radius required for iso transformation")
    }
    
    # Stata equivalent:
    # if "`transformation'"=="iso" & `radius'<=0 {
    #   di as error "Radius must be positive."
    #   exit
    # }
    if (transformation == "iso" && radius <= 0) {
        stop("Radius must be positive")
    }
    
    # Stata equivalent:
    # if "`transformation'"=="cluster" & "`clustvar'"=="" {
    #   di as error "Clustvar missing."
    #   exit
    # }
    if (transformation == "cluster" && is.null(clustvar)) {
        stop("Clustvar required for cluster transformation")
    }
    
    # Stata equivalent:
    # if "`transformation'"!="cluster" & "`clustvar'"!="" {
    #   di as error "Option clustvar only allowed with transformation(cluster)."
    #   exit
    # }
    if (transformation != "cluster" && !is.null(clustvar)) {
        stop("Clustvar only allowed with transformation(cluster)")
    }
    
    # Check for coordinates - no direct Stata equivalent
    # Stata expects specific columns while R allows specification
    if (!all(c(lat, lon) %in% names(data))) {
        stop(paste("Data frame must contain columns:", lat, "and", lon))
    }
    
    base_sample <- rep(TRUE, nrow(data))
    
    # Check for required coordinates
    if (!all(c(lat, lon) %in% names(data))) {
        stop(paste("Data frame must contain columns:", lat, "and", lon))
    }
    
    # Extract coordinates and standardize names
    coords <- data[, c(lat, lon)]
    names(coords) <- c("lat", "lon")  # Rename to standard names for internal use
    
    # Convert cluster variable if needed
    if (transformation == "cluster") {
        # Convert cluster variable to consecutive numbers (like egen group)
        cluster <- as.numeric(factor(data[[clustvar]]))
    } else {
        cluster <- NULL
    }
    
    # Prepare result dataframe
    result <- data
    
    # Stata equivalent:
    # foreach name of local varlist {
    for (var in vars) {
        # Handle missing values
        # Stata equivalent:
        # qui gen `touse2' = `touse'
        # if "`separately'"!="" {
        #   qui replace `touse2' = 0 if missing(`name')
        # }
        # if "`transformation'"=="cluster" {
        #   qui replace `touse2' = 0 if missing(`clustervar')
        # }
        # Start with base sample for each variable (equivalent to qui gen `touse2' = `touse')
        use_rows <- base_sample
        
        # Add variable-specific filtering
        if (separately) {
            use_rows <- use_rows & !is.na(data[[var]])
        }
        if (transformation == "cluster") {
            use_rows <- use_rows & !is.na(data[[clustvar]])
        }
        
        # Skip if no observations to use
        if (sum(use_rows) == 0) {
            warning(paste("No valid observations for variable:", var))
            next
        }

        # Create coordinates specific to this variable's non-missing pattern
        # Stata equivalent:
        # mata: s = get_s_matrix("`touse2'", "`latlong'")
        var_coords <- get_s_matrix(coords[use_rows, ], latlong) # identical
        
        if (transformation == "cluster") {
          # Get cluster values for the observations to use
          cluster_matrix <- get_cluster_matrix(data[[clustvar]], use_rows) # identical
        } else {
          cluster_matrix <- 0
        }
        
        H <- make_transform(var_coords, transformation, radius, cluster_matrix, latlong = latlong) # minimally different
        
        # Apply transformation with the variable-specific matrix
        # Stata equivalent:
        # mata: hy = transform("`name'", H, "`touse2'", "`transformation'")
        
        transformed <- transform(data, var, H, transformation, use_rows)  # identical
        
        # Insert results
        # Stata equivalent:
        # quietly mata: st_addvar("double", "`prefix'`name'")
        # mata: st_store(., "`prefix'`name'", "`touse2'", hy)
        result_vec <- rep(NA, nrow(data))
        result_vec[use_rows] <- transformed
        result[[paste0(prefix, var)]] <- result_vec
    }
    
    return(result)
}