# ToDO:

# multiply ha area by actual shrub/treecover before
# estimating SAL costs

# refine SAL costs to reflect slope



BCN_target_df = function(
    data,
    original_col,
    mask_col,
    replacement_col,
    target_col = "random", # Optional: if not provided, random assignment is used
    threshold = NULL,
    budget = NULL,
    relative = TRUE,
    target_decreasing = FALSE,
    year = NULL,
    verbose = FALSE
) {
  # --- Input Validation ---
  # Check if the input is a data.frame
  if (!inherits(data, "data.frame")) {
    stop("Input 'data' must be a data.frame object.")
  }
  
  # Define mandatory column arguments
  mandatory_cols_args = c("original_col", "mask_col", "replacement_col")
  # Check if all mandatory column arguments are provided
  for (arg in mandatory_cols_args) {
    if (is.null(get(arg))) {
      stop(paste("Argument '", arg, "' must be provided.", sep = ""))
    }
  }
  
  # Check if all specified column names exist in the dataframe
  # Collect all column names that are expected to be in the dataframe
  all_specified_cols = c(original_col, mask_col, replacement_col)
  if (target_col != "random") {
    all_specified_cols = c(all_specified_cols, target_col)
  }
  
  # Check if all collected column names are actually present in the dataframe
  if (!all(all_specified_cols %in% names(data))) {
    missing_cols = setdiff(all_specified_cols, names(data))
    stop(paste("The following specified columns are missing from the dataframe:", paste(missing_cols, collapse = ", ")))
  }
  
  # Check if either 'threshold' or 'budget' is provided
  if (is.null(threshold) && is.null(budget)) {
    stop("Either 'threshold' or 'budget' must be provided.")
  }
  
  # Check if 'relative' is a boolean value
  if (!is.logical(relative)) {
    stop("'relative' must be a boolean value (TRUE or FALSE).")
  }
  
  # Ensure 'mask_col' column contains only 0s, 1s, or NAs
  if (!all(data[[mask_col]] %in% c(0, 1, NA))) {
    stop(paste("Column '", mask_col, "' must contain only 0s, 1s, or NAs.", sep = ""))
  }
  
  # --- Calculate Intervention Cells ---
  # Identify rows where the mask column is 1 (eligible for intervention)
  # and exclude NA values from the mask.
  original_intervention_indices = which(data[[mask_col]] == 1 & !is.na(data[[mask_col]]))
  intervention_cells = length(original_intervention_indices)
  
  # If no cells are marked for intervention, return the original data with a message
  if (intervention_cells == 0) {
    if (verbose) {
      message("No cells found with 'intervention_mask' = 1. Returning original data.")
    }
    # Determine column name suffix for the new columns
    col_suffix = if (relative) {
      if (!is.null(threshold)) paste0("_", target_col, "_t", as.character(threshold * 100), "_y", year) else ""
    } else {
      if (!is.null(budget)) paste0("_", target_col, "_b", as.character(budget), "_y", year) else ""
    }
    data[[paste0("new_mask", col_suffix)]] = 0 # All zeros for new mask
    data[[paste0("modified_values", col_suffix)]] = data[[original_col]] # Original values
    return(data)
  }
  
  # --- Determine Target Cell Count ---
  # Calculate the absolute number of target cells based on 'relative' and 'threshold'/'budget'
  if (relative) {
    # If 'relative' is TRUE, 'threshold' must be between 0 and 1
    if (is.null(threshold) || (threshold <= 0 || threshold > 1)) {
      stop("If 'relative' is TRUE, 'threshold' (between 0 and 1) must be provided.")
    }
    target_cells = round(intervention_cells * threshold)
  } else {
    # If 'relative' is FALSE, 'budget' must be provided
    if (is.null(budget)) {
      stop("If 'relative' is FALSE, 'budget' must be provided.")
    }
    target_cells = budget
  }
  
  # Adjust 'target_cells' to be within valid bounds (0 to 'intervention_cells')
  if (target_cells > intervention_cells) {
    warning(paste("Target cells (", target_cells, ") exceed total intervention cells (", intervention_cells, "). Setting target_cells to intervention_cells."))
    target_cells = intervention_cells
  }
  if (target_cells < 0) {
    warning(paste("Target cells (", target_cells, ") is negative. Setting target_cells to 0."))
    target_cells = 0
  }
  
  # --- Generate Column Names for Output ---
  # Create a suffix for the new column names based on threshold or budget
  col_suffix = if (relative) {
    paste0("_", target_col, "_t", as.character(threshold * 100), "_y", year)
  } else {
    paste0("_", target_col, "_b", as.character(budget), "_y", year)
  }
  new_mask_col_name = paste0("new_mask", col_suffix)
  modified_values_col_name = paste0("modified_values", col_suffix)
  
  
  # --- Initialize New Columns ---
  # Create the new mask column, initially all zeros
  data[[new_mask_col_name]] = 0
  
  # Create the modified values column, initially a copy of original_col
  data[[modified_values_col_name]] = data[[original_col]]
  
  # --- Select Target Cells for Modification ---
  target_cells_global_indices = integer(0) # Initialize an empty vector for selected indices
  
  # Check if 'target_col' is provided and has non-NA values for intervention cells
  # This determines whether to use random selection or ordered selection
  has_valid_target_surface = (target_col != "random") &&
    any(!is.na(data[[target_col]][original_intervention_indices]))
  
  if (!has_valid_target_surface) {
    # If no valid 'target_col' values, randomly select target cells from intervention mask
    if (target_cells > 0) {
      target_cells_global_indices = sample(
        original_intervention_indices, # Sample from the indices of eligible intervention cells
        size = target_cells,
        replace = FALSE # Do not sample with replacement
      )
    }
  } else {
    # Create a temporary data structure to link original indices with their target_surface values
    temp_selection_data = data.frame(
      original_index = original_intervention_indices,
      target_value = data[[target_col]][original_intervention_indices]
    )
    
    # Remove rows where 'target_value' is NA, as they cannot be ordered
    temp_selection_data = temp_selection_data[!is.na(temp_selection_data$target_value), ]
    
    if (nrow(temp_selection_data) == 0) {
      # If all 'target_col' values for intervention cells are NA after filtering,
      # fall back to random selection.
      if (verbose) {
        message("No valid 'target_col' values within intervention mask for ordering. Randomly selecting target cells.")
      }
      if (target_cells > 0) {
        target_cells_global_indices = sample(
          original_intervention_indices,
          size = target_cells,
          replace = FALSE
        )
      }
    } else {
      # Order the temporary data by 'target_value' in decreasing order
      ordered_temp_data = temp_selection_data[order(temp_selection_data$target_value, decreasing = target_decreasing), ]
      
      # Select the 'original_index' for the top 'target_cells' rows
      if (target_cells > 0) {
        target_cells_global_indices = ordered_temp_data$original_index[1:min(target_cells, nrow(ordered_temp_data))]
      }
    }
  }
  
  # --- Apply Replacement Values and Update New Mask ---
  # If any target cells were selected, update their 'modified_values' with 'replacement_values'
  # and set the new mask column to 1 for these cells.
  if (length(target_cells_global_indices) > 0) {
    data[[modified_values_col_name]][target_cells_global_indices] = data[[replacement_col]][target_cells_global_indices]
    data[[new_mask_col_name]][target_cells_global_indices] = 1
  }
  
  # --- Verbose Output ---
  if (verbose) {
    message(paste("Number of target cells selected for modification:", length(target_cells_global_indices)))
    message(paste("Total number of intervention-eligible cells considered:", intervention_cells))
    message(paste("Calculated target cells threshold/budget:", target_cells))
  }
  
  # Return the dataframe with the new columns
  return(data)
}
