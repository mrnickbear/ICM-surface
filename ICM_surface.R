#Goal - convert a Civil 3D surface into layer cake bands that don't overlap, for ICM ground modification
#
#test
#Civil 3D
#1 - turn on border and contours with desired interval
#2 - Select surface, Surface Tools panel, Extract from Surface, Border and Major Contour, save drawing.


library(sf)
library(dplyr)
library(rmapshaper) #for ms_simplify which maintains topology.  sf::st_simplify would create gaps and overlaps.   #https://cran.r-project.org/web/packages/rmapshaper/vignettes/rmapshaper.html
library(mapview)

#For the Z value extraction
library(tidyr)
library(purrr)
library(glue)

# 1. Load the DXF file
# Note: Ensure the 'Layer' or 'Elevation' attribute name matches your file

# Path that can be copy/pasted (escaping removed)
# S:\HYD\Projects - studies\STUDIES\R_Ochoa\1_Project\1.0_Project-By-Name\2817J_6thStreet Reroute\4th&King-Yard\

dxf_path <- paste0("S:\\HYD\\Projects - studies\\STUDIES\\R_Ochoa\\1_Project\\1.0_Project-By-Name\\2817J_6thStreet Reroute\\4th&King-Yard\\",
                   "References\\Public Works Coordination 20251222\\4th-King-Yard_Grading-Drainage-Files\\",
                   "XCGRADA_nb_linesOnly_worldCS.dxf")

raw_data <- st_read(dxf_path) %>% st_set_crs(2227)

# dput(
# raw_data[2,]
# )

# %>% st_zm(raw_data) #st_zm() to drop Z values


# # 2. Separate Contours and Boundary
# # Assuming standard DXF attributes; adjust 'Layer' names as needed
# contours <- raw_data %>% filter(Linetype == "Continuous") 
# boundary_3d <- raw_data %>% filter(Linetype == "HIDDENX2")
# 
# 
# # 3. Simplify Contours
# # 'keep' defines the proportion of points to retain. 
# # sys = TRUE uses the system mapshaper if installed (faster for large files)
# contours_simple <- ms_simplify(contours, keep = 0.05, keep_shapes = TRUE)



# =================================================================
# STEP 0: Calculate Z Stats & Isolate Boundary
# =================================================================
raw_data <- raw_data %>%
  mutate(z_stats = map(st_geometry(.), function(geom) {
    coords <- st_coordinates(geom)
    if ("Z" %in% colnames(coords)) {
      list(
        avg_z = mean(coords[, "Z"], na.rm = TRUE),
        min_z = min(coords[, "Z"], na.rm = TRUE),
        max_z = max(coords[, "Z"], na.rm = TRUE)
      )
    } else {
      list(avg_z = NA_real_, min_z = NA_real_, max_z = NA_real_)
    }
  })) %>%
  unnest_wider(z_stats) %>%
  st_as_sf() %>%
  mutate(z_variance = max_z - min_z) %>%
  st_zm(drop = TRUE, what = "ZM")

boundary_data <- raw_data %>% filter(z_variance > 0)
contour_data  <- raw_data %>% filter(z_variance == 0 | is.na(z_variance))

# =================================================================
# STEP 1 & 2: Slice Boundary into Cells using Core SF
# =================================================================
# Extract both as raw line segments
boundary_line <- st_geometry(boundary_data) %>% st_cast("LINESTRING")
contour_lines <- st_cast(contour_data, "LINESTRING")

extend_contours_to_boundary <- function(contour_lines, boundary_line) {
  boundary_geom <- st_union(boundary_line)
  boundary_bbox <- st_bbox(boundary_geom)
  boundary_diagonal <- sqrt(
    (boundary_bbox[["xmax"]] - boundary_bbox[["xmin"]]) ^ 2 +
      (boundary_bbox[["ymax"]] - boundary_bbox[["ymin"]]) ^ 2
  )
  boundary_tolerance <- boundary_diagonal * 0.0001
  boundary_extension <- boundary_diagonal * 0.001
  contour_geom <- st_geometry(contour_lines)
  crs <- st_crs(contour_geom)
  
  extended_lines <- lapply(seq_along(contour_geom), function(i) {
    coords <- st_coordinates(contour_geom[[i]])[, c("X", "Y"), drop = FALSE]
    
    if (nrow(coords) >= 2) {
      start_point <- st_sfc(st_point(coords[1, ]), crs = crs)
      end_point <- st_sfc(st_point(coords[nrow(coords), ]), crs = crs)
      
      if (as.numeric(st_distance(start_point, boundary_geom)) <= boundary_tolerance) {
        direction <- coords[1, ] - coords[2, ]
        direction_length <- sqrt(sum(direction ^ 2))
        if (direction_length > 0) {
          coords[1, ] <- coords[1, ] + boundary_extension * direction / direction_length
        }
      }
      
      if (as.numeric(st_distance(end_point, boundary_geom)) <= boundary_tolerance) {
        direction <- coords[nrow(coords), ] - coords[nrow(coords) - 1, ]
        direction_length <- sqrt(sum(direction ^ 2))
        if (direction_length > 0) {
          coords[nrow(coords), ] <- coords[nrow(coords), ] + boundary_extension * direction / direction_length
        }
      }
    }
    
    st_linestring(coords)
  })
  
  do.call(st_sfc, c(extended_lines, list(crs = crs)))
}

contour_lines_for_polygonizing <- extend_contours_to_boundary(contour_lines, boundary_line)

# Combine the boundary outer line and the open contour lines.
# Extending boundary-touching contour ends just past the boundary lets polygonize
# create closed boundary/contour bands for contours that do not close on themselves.
line_network <- st_union(c(contour_lines_for_polygonizing, boundary_line))

# Polygonize the combined network into independent jigsaw cells
sliced_surface <- line_network %>% 
  st_polygonize() %>% 
  st_collection_extract("POLYGON") %>% 
  st_as_sf() %>% 
  mutate(cell_id = row_number())

# =================================================================
# STEP 3: Assign Elevations Based on Bordering Lines
# =================================================================
# Because the lines are open, a cell sits *between* two contours.
# We will find contour lines that touch each cell's edge.
cell_edges <- sliced_surface %>%
  st_set_geometry(st_boundary(st_geometry(.)))

cell_line_match <- st_join(cell_edges, contour_data, join = st_intersects)

cell_elevations <- cell_line_match %>% 
  st_drop_geometry() %>% 
  group_by(cell_id) %>% 
  summarize(
    # Assign the cell to the lower bounding contour line it touches
    z_layer = min(avg_z, na.rm = TRUE),
    .groups = "drop"
  ) %>% 
  # If a cell only touches the outer boundary and no contour lines, 
  # it's the absolute lowest base cell.
  mutate(z_layer = ifelse(is.infinite(z_layer), min(contour_data$avg_z, na.rm = TRUE), z_layer))

# Join the true elevations back to our sliced pieces
terrain_cells_valued <- sliced_surface %>% 
  left_join(cell_elevations, by = "cell_id")

# =================================================================
# STEP 4: Prepare and Hollow Out the Cake Layers
# =================================================================
outer_polygon_area <- function(geometry) {
  polygons <- st_cast(st_geometry(geometry), "POLYGON")
  crs <- st_crs(polygons)
  
  sum(vapply(seq_along(polygons), function(i) {
    polygon <- polygons[[i]]
    as.numeric(st_area(st_sfc(st_polygon(list(polygon[[1]])), crs = crs)))
  }, numeric(1)))
}

cast_layer_polygons <- function(layers) {
  polygon_layers <- lapply(seq_len(nrow(layers)), function(i) {
    layer_geometry <- st_geometry(layers[i, ])
    geometry_type <- unique(as.character(st_geometry_type(layer_geometry)))
    
    if (all(!geometry_type %in% c("POLYGON", "MULTIPOLYGON", "GEOMETRYCOLLECTION"))) {
      return(NULL)
    }
    
    polygons <- tryCatch(
      suppressWarnings(st_collection_extract(layer_geometry, "POLYGON")),
      error = function(e) layer_geometry
    )
    polygons <- st_cast(polygons, "POLYGON", warn = FALSE)
    polygons <- polygons[!st_is_empty(polygons)]
    
    if (length(polygons) == 0) {
      return(NULL)
    }
    
    st_sf(
      z_layer = rep(layers[i, ]$z_layer, length(polygons)),
      geometry = polygons
    )
  })
  
  polygon_layers <- polygon_layers[!vapply(polygon_layers, is.null, logical(1))]
  
  if (length(polygon_layers) == 0) {
    return(st_sf(z_layer = numeric(), geometry = st_sfc(crs = st_crs(layers))))
  }
  
  do.call(rbind, polygon_layers)
}

remove_layer_overlaps <- function(layers) {
  if (nrow(layers) == 0) {
    return(layers)
  }
  
  layers <- cast_layer_polygons(layers)
  
  if (nrow(layers) == 0) {
    return(layers)
  }
  
  layers <- layers %>%
    mutate(area = as.numeric(st_area(geometry))) %>%
    arrange(area)
  
  clean_layers <- list()
  occupied_geometry <- st_sfc(crs = st_crs(layers))
  
  for (i in seq_len(nrow(layers))) {
    current_layer <- layers[i, ]
    current_geom <- st_make_valid(st_geometry(current_layer))
    
    if (length(occupied_geometry) > 0) {
      overlaps <- st_intersects(occupied_geometry, current_geom, sparse = FALSE)[, 1]
      
      if (any(overlaps)) {
        overlap_mask <- st_make_valid(st_union(occupied_geometry[overlaps]))
        current_geom <- st_difference(current_geom, overlap_mask)
        current_geom <- current_geom[!st_is_empty(current_geom)]
        
        if (length(current_geom) > 0) {
          current_geom <- st_cast(current_geom, "POLYGON")
          current_geom <- current_geom[!st_is_empty(current_geom)]
        }
      }
    }
    
    if (length(current_geom) > 0) {
      current_layer <- st_sf(
        z_layer = rep(layers[i, ]$z_layer, length(current_geom)),
        geometry = current_geom
      )
      clean_layers[[length(clean_layers) + 1]] <- current_layer
      occupied_geometry <- c(occupied_geometry, st_geometry(current_layer))
    }
  }
  
  do.call(rbind, clean_layers) %>%
    mutate(area = as.numeric(st_area(geometry)))
}

clean_simplified_layers <- function(simplified_layers, original_layers, min_area, snap_tolerance) {
  simplified_layers <- cast_layer_polygons(st_make_valid(simplified_layers)) %>%
    mutate(area = as.numeric(st_area(geometry))) %>%
    filter(area > min_area)
  
  if (nrow(simplified_layers) == 0) {
    return(simplified_layers)
  }
  
  snapped_geometry <- st_snap(
    st_geometry(simplified_layers),
    st_union(st_geometry(simplified_layers)),
    tolerance = snap_tolerance
  )
  
  simplified_layers <- simplified_layers %>%
    st_set_geometry(st_make_valid(snapped_geometry)) %>%
    cast_layer_polygons() %>%
    mutate(area = as.numeric(st_area(geometry))) %>%
    filter(area > min_area)
  
  simplified_layers <- remove_layer_overlaps(simplified_layers)
  
  original_footprint <- st_union(st_geometry(original_layers))
  simplified_footprint <- st_union(st_geometry(simplified_layers))
  gaps <- st_difference(original_footprint, simplified_footprint)
  gaps <- gaps[!st_is_empty(gaps)]
  
  if (length(gaps) > 0) {
    gap_parts <- st_sf(geometry = suppressWarnings(st_collection_extract(gaps, "POLYGON"))) %>%
      st_cast("POLYGON", warn = FALSE) %>%
      mutate(area = as.numeric(st_area(geometry))) %>%
      filter(area > 0)
    
    if (nrow(gap_parts) > 0) {
      gap_buffers <- st_buffer(st_geometry(gap_parts), snap_tolerance)
      neighbor_matrix <- st_intersects(gap_buffers, simplified_layers, sparse = FALSE)
      
      for (gap_index in seq_len(nrow(gap_parts))) {
        neighbors <- which(neighbor_matrix[gap_index, ])
        
        if (length(neighbors) >= 2) {
          gap_boundary <- st_boundary(st_geometry(gap_parts[gap_index, ]))
          neighbor_boundaries <- st_boundary(st_geometry(simplified_layers[neighbors, ]))
          shared_lengths <- as.numeric(st_length(st_intersection(gap_boundary, neighbor_boundaries)))
          
          if (all(shared_lengths == 0)) {
            shared_lengths <- as.numeric(st_area(st_intersection(gap_buffers[gap_index], simplified_layers[neighbors, ])))
          }
          
          target_index <- neighbors[which.max(shared_lengths)]
          st_geometry(simplified_layers)[target_index] <- st_union(
            st_geometry(simplified_layers[target_index, ]),
            st_geometry(gap_parts[gap_index, ])
          )
        }
      }
    }
  }
  
  remove_layer_overlaps(simplified_layers) %>%
    filter(area > min_area)
}

# Keep each polygon piece separate so multipart layers are not collapsed.
solid_layers <- terrain_cells_valued %>% 
  select(z_layer) %>%
  rename(geometry = x) %>%
  st_cast("POLYGON") %>%
  mutate(
    area = as.numeric(st_area(geometry)),
    outer_area = vapply(seq_along(geometry), function(i) {
      outer_polygon_area(geometry[i])
    }, numeric(1))
  ) %>%
  filter(area>10) %>%
  arrange(outer_area)

# Punch out previously kept smaller pieces to create true non-overlapping rings/ribbons
hollow_layers <- list()
occupied_geometry <- st_sfc(crs = st_crs(solid_layers))

for (i in seq_len(nrow(solid_layers))) {
  
  #i <- 26
  cat(glue("punching out i=",i,"\n\n"))
  current_layer <- solid_layers[i, ]
  current_geom <- st_make_valid(st_geometry(current_layer))
  
  if (length(occupied_geometry) > 0) {
    overlaps <- st_intersects(occupied_geometry, current_geom, sparse = FALSE)[, 1]
    
    if (any(overlaps)) {
      internal_mask <- st_make_valid(st_union(occupied_geometry[overlaps]))
      difference <- st_difference(current_geom, internal_mask)
      difference <- difference[!st_is_empty(difference)]
      
      if (length(difference) > 0) {
        difference <- st_cast(difference, "POLYGON")
        difference <- difference[!st_is_empty(difference)]
      }
      
      if (length(difference) > 0) {
        current_layer <- st_sf(
          z_layer = rep(solid_layers[i, ]$z_layer, length(difference)),
          geometry = difference
        )
      } else {
        current_layer <- NULL
      }
    }
  }
  
  if (!is.null(current_layer)) {
    current_layer <- current_layer[, c("z_layer", "geometry")]
    hollow_layers[[length(hollow_layers) + 1]] <- current_layer
    occupied_geometry <- c(occupied_geometry, st_geometry(current_layer))
  }
  
  # mapview(current_layer)
  # mapview(difference)
  # mapview(internal_mask)+mapview( current_geom)
}

# 1. Rename all geometry columns to 'geometry'
hollow_layers_clean <- lapply(hollow_layers, function(df) {
  st_geometry(df) <- "geometry"
  return(df)
})

# 2. Select only common columns (e.g., z_layer and geometry)
hollow_layers_clean <- lapply(hollow_layers_clean, function(df) {
  df[, c("z_layer", "geometry")]
})

# Re-combine and calculate the true individual ring areas
final_layer_cake <- do.call(rbind, hollow_layers_clean) %>% 
  # st_as_sf() %>% 
  mutate(area = as.numeric(st_area(geometry))) %>% filter(area > 10)



# =================================================================
# STEP 5: Verification Print & Plot
# =================================================================
print(final_layer_cake %>% st_drop_geometry() %>% select(z_layer, area))

# plot(final_layer_cake["z_layer"], 
#      pal = terrain.colors(nrow(final_layer_cake)), 
#      main = "True Layer-Cake Bands (Open Contours Accounted For)")

# mapview(final_layer_cake[4,])
mapview(final_layer_cake, z = "z_layer")


simplify_buffer_distance <- 1
simplify_min_area <- 10
simplify_snap_tolerance <- simplify_buffer_distance * 2

simp_final <- ms_simplify(final_layer_cake %>% 
                            # Dual buffers to remove thin parts of ring
                            st_buffer(dist = -simplify_buffer_distance) %>%
                            st_buffer(dist = simplify_buffer_distance),
                          keep = 0.1) %>%
  clean_simplified_layers(
    original_layers = final_layer_cake,
    min_area = simplify_min_area,
    snap_tolerance = simplify_snap_tolerance
  )

mapview(simp_final, z = "z_layer")


# mapview(raw_data)
# mapview(contours_simple)

# mapview(final_layer_cake, z = "z_layer")
# 
# mapview(cake_layers, z = "z_layer")
# 
# mapview( arrange(desc(solid_layers)), z = "z_layer")
# 
# 
# mapview( 
#   filter(solid_layers, z_layer  < 10)
#   )

# flc_areas <- solid_layers %>%
#   mutate(
#     # Calculate area and strip the 'units' attribute for a clean numeric column
#     area = as.numeric(st_area(geometry))
#   )
