lut_path <- "C:/Users/Nelson/Documents/Projects/map_generator"
lut_filename <- "feature_lut.csv"
feature_lut <- read.csv(paste0(lut_path, "/", lut_filename),
                        stringsAsFactors = FALSE)

features <- unique(feature_lut$feature)

x_dim <- 10
y_dim <- 10

initialize_map <- function(x_dim,
                           y_dim,
                           features){
  setNames(object = lapply(X = features,
                           x_dim= x_dim,
                           y_dim = y_dim,
                           FUN = function(X, x_dim, y_dim) {
                             matrix(data = X,
                                    nrow = y_dim,
                                    ncol = x_dim)
                           }),
           nm = features)
}

map_potentials_list <- initialize_map(x_dim = x_dim,
                                      y_dim = y_dim,
                                      features = features)

# Make the map
map <- matrix(data = NA,
              nrow = y_dim,
              ncol = x_dim)

while (any(is.na(map))) {
  # Get the coordinates we'll be centering on
  #!!!!!!! Use vectors of coords that correspond to NA cells in the map
  unresolved_cell_coords <- which(x = is.na(map),
                                  arr.ind = TRUE)
  random_coord_index <- sample(x = seq_len(nrow(unresolved_cell_coords)),
                               size = 1)
  current_x <- as.vector(unresolved_cell_coords[random_coord_index, "col"])
  current_y <- as.vector(unresolved_cell_coords[random_coord_index, "row"])
  
  # First we pick a value
  cell_potentials <- sapply(X = map_potentials_list,
                            x = current_x,
                            y = current_y,
                            FUN = function(X, x, y) {
                              X[y, x]
                            })
  cell_potentials <- cell_potentials[!is.na(cell_potentials)]
  
  # Then we write it into the map
  if (length(cell_potentials) > 0) {
    new_cell_value <- sample(x = cell_potentials,
                             size = 1)
    # Just to de-emphasize the blanks, we'll only accept it as a value if we draw
    # it twice in a row. There'll be plenty of blanks surrounding lines
    if (new_cell_value == "blank") {
      sample(x = cell_potentials,
             size = 1)
    }
    map[current_y, current_x] <- new_cell_value
  } else {
    # In case we end up with impossible cells that have no potential
    map[current_y, current_x] <- "blank"
  }
  
  current_cell_value <- map[current_y, current_x]
  
  # And now we collapse the surrounding cells' potentials
  # All diagonals will be set to "blank"
  diagonal_coords <- list("left_down" = c("x" = current_x - 1,
                                          "y" = current_y - 1),
                          "right_down" = c("x" = current_x + 1,
                                           "y" = current_y - 1),
                          "left_up" = c("x" = current_x - 1,
                                        "y" = current_y + 1),
                          "right_up" = c("x" = current_x + 1,
                                         "y" = current_y + 1))
  
  # Make sure the diagonals are all in-bounds
  diagonal_coords <- lapply(X = diagonal_coords,
                            x_max = x_dim,
                            y_max = y_dim,
                            FUN = function(X, x_max, y_max) {
                              x_valid <- 0 < X[["x"]] & X[["x"]] < x_max
                              y_valid <- 0 < X[["y"]] & X[["y"]] < y_max
                              if (x_valid & y_valid) {
                                c("x" = X[["x"]],
                                  "y" = X[["y"]])
                              } else {
                                NULL
                              }
                            })
  diagonal_coords <- diagonal_coords[!sapply(diagonal_coords, is.null)]
  
  for (coords in diagonal_coords) {
    map_potentials_list <- lapply(X = features,
                                  potentials_list = map_potentials_list,
                                  x = coords[["x"]],
                                  y = coords[["y"]],
                                  FUN = function(X, potentials_list, x, y) {
                                    current_matrix <- potentials_list[[X]]
                                    if (X != "blank") {
                                      current_matrix[y, x] <- NA
                                    }
                                    current_matrix
                                  })
    names(map_potentials_list) <- features
  }
  
  orthogonal_coords <- list("left" = c("x" = current_x - 1,
                                       "y" = current_y),
                            "right" = c("x" = current_x + 1,
                                        "y" = current_y),
                            "up" = c("x" = current_x,
                                     "y" = current_y + 1),
                            "down" = c("x" = current_x,
                                       "y" = current_y - 1))
  
  orthogonal_coords <- lapply(X = orthogonal_coords,
                              x_max = x_dim,
                              y_max = y_dim,
                              FUN = function(X, x_max, y_max) {
                                x_valid <- 0 < X[["x"]] & X[["x"]] < x_max
                                y_valid <- 0 < X[["y"]] & X[["y"]] < y_max
                                if (x_valid & y_valid) {
                                  c("x" = X[["x"]],
                                    "y" = X[["y"]])
                                } else {
                                  NULL
                                }
                              })
  orthogonal_coords <- orthogonal_coords[!sapply(orthogonal_coords, is.null)]
  
  # For each direction orthogonally
  for (direction in names(orthogonal_coords)) {
    # Going in that direction from our current cell, is there an outgoing connection?
    current_outgoing_direction_valid <- feature_lut[["outgoing_connection"]][feature_lut$feature == current_cell_value & feature_lut$direction == direction]
    # Which features don't match the above answer?
    !identical(feature_lut$incoming_connection, current_outgoing_direction_valid)
    current_invalid_feature_indices <- sapply(X = feature_lut$incoming_connection,
                                              value = current_outgoing_direction_valid,
                                              FUN = function(X, value) {
                                                !identical(X, value)
                                              })
    current_invalid_features <- feature_lut$feature[current_invalid_feature_indices]
    
    # For each feature that's invalid, 
    for (feature in current_invalid_features) {
      map_potentials_list[[feature]][orthogonal_coords[[direction]][["y"]], orthogonal_coords[[direction]][["x"]]] <- NA
    }
  }
}

plotting_df <- do.call(rbind,
                       lapply(X = unique(as.vector(map)),
                              map = map,
                              FUN = function(X, map) {
                                output <- as.data.frame(x = which(map == X, arr.ind = TRUE))
                                output$feature <- X
                                names(output) <- c("y", "x", "feature")
                                output
                              }))

plotting_df$color <- "gray90"
plotting_df$color[plotting_df$feature == "blank"] <- "white"

ggplot() +
  geom_raster(data = plotting_df,
              aes(x = x,
                  y = y,
                  fill = color))
