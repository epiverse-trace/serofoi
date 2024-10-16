library(ggplot2)
library(rlang)
library(purrr)


# Dear Code Reviewer: If I'm reinventing the wheel with these functions, 
# please let me know in the PR comments :)



# Utility functions to extract relevant data from the plot, for comparison in each test
extract_plot_data <- function(plot) {
  list(
    classes = class(plot),
    layers = extract_layers(plot$layers),
    coordinates = list(
      limits = plot$coordinates$limits
    ),
    labels = plot$labels
  )
}

extract_layers <- function(layers) {
  map(layers, function(layer) {
    # Extract mappings
    list(
      mapping =
        if (is.null(layer$mapping)) {
          NULL
        } else {
          map(layer$mapping, \(x) {
            quo_name(get_expr(x))
          })
        },
      geom_params =
        if (is.null(layer$geom_params)) {
          NULL
        } else {
          layer$geom_params[!names(layer$geom_params) %in% c("grob")]
        }
    )
  })
}

# Custom function to compare lists with tolerance for numeric values
compare_lists_with_tolerance <- function(list1, list2, tol, path = "") {
  
  # Check if both are lists
  if (is.list(list1) && is.list(list2)) {
    # Compare lengths
    if (length(list1) != length(list2)) {
      stop(sprintf("Different lengths at %s: length %d vs %d", path, length(list1), length(list2)))
    }
    
    # Compare each element recursively
    for (name in names(list1)) {
      new_path <- paste0(path, "$", name)
      
      # Check if name exists in both lists
      if (!name %in% names(list2)) {
        stop(sprintf("Name '%s' missing in second list at %s", name, new_path))
      }
      
      # Recursively compare elements
      if (!compare_lists_with_tolerance(list1[[name]], list2[[name]], tol = tol, path = new_path)) {
        return(FALSE)
      }
    }
    return(TRUE)
  }
  
  # If numeric, compare with tolerance
  if (is.numeric(list1) && is.numeric(list2)) {
    if (!all(abs(list1 - list2) <= tol)) {
      stop(sprintf("Numeric mismatch at %s: %s vs %s with tolerance %e", path, list1, list2, tol))
    }
    return(TRUE)
  }
  
  # Otherwise, compare directly
  if (!identical(list1, list2)) {
    stop(sprintf("Mismatch at %s: %s vs %s", path, list1, list2))
  }
  return(TRUE)
}

# Using expect_true to assert comparison
expect_lists_equal_with_tolerance <- function(list1, list2, tol = 1e-2) {
  tryCatch({
    expect_true(compare_lists_with_tolerance(list1, list2, tol = tol))
  }, error = function(e) {
    fail(e$message)
  })
}