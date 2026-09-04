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

# 1. Load the DXF file
# Note: Ensure the 'Layer' or 'Elevation' attribute name matches your file

# Path that can be copy/pasted (escaping removed)
# S:\HYD\Projects - studies\STUDIES\R_Ochoa\1_Project\1.0_Project-By-Name\2817J_6thStreet Reroute\4th&King-Yard\

dxf_path <- paste0("S:\\HYD\\Projects - studies\\STUDIES\\R_Ochoa\\1_Project\\1.0_Project-By-Name\\2817J_6thStreet Reroute\\4th&King-Yard\\",
                   "References\\Public Works Coordination 20251222\\4th-King-Yard_Grading-Drainage-Files\\",
                   "XCGRADA_nb_linesOnly_worldCS.dxf")

test_dxf <- "XCGRADA_nb_linesOnly_worldCS.dxf"

raw_data <- st_read(test_dxf) %>% st_set_crs(2227)

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
  mutate(z_variance = max_z - min_z)

boundary_data <- raw_data %>% filter(z_variance > 0)
contour_data  <- raw_data %>% filter(z_variance == 0 | is.na(z_variance))

# =================================================================
# STEP 1 & 2: Slice Boundary into Cells using Core SF
# =================================================================
# Extract both as raw line segments
boundary_line <- st_geometry(boundary_data) %>% st_cast("LINESTRING")
contour_lines <- st_cast(contour_data, "LINESTRING")

# Combine the boundary outer line and the open contour lines.
# This forces the open lines to touch the boundary line, sealing the edges.
line_network <- st_union(st_geometry(contour_lines), boundary_line)

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
# We will find all contour lines that touch each cell's edge.
cell_line_match <- st_join(sliced_surface, contour_data, join = st_intersects)

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
# STEP 4: Dissolve and Hollow Out the Cake Layers
# =================================================================
# Now, grouping by z_layer will yield distinct, staggered solid shapes!
solid_layers <- terrain_cells_valued %>% 
  group_by(z_layer) %>% 
  summarize(geometry = st_union(x), .groups = "drop") %>% 
  arrange(z_layer)

# Punch out the nested pieces to create true non-overlapping rings/ribbons
hollow_layers <- lapply(1:nrow(solid_layers), function(i) {
  
  i <- 5
  current_layer <- solid_layers[i, ]
  
  # Find other layers whose centroids are inside the current layer
  other_layers <- solid_layers[-i, ]
  centroids <- st_centroid(other_layers)
  is_inside <- st_intersects(centroids, current_layer, sparse = FALSE)[, 1]
  
  if (any(is_inside)) {
    internal_mask <- st_union(other_layers[is_inside, ])
    current_layer$geometry <- st_difference(
      st_make_valid(current_layer$geometry), 
      st_make_valid(internal_mask)
    )
  }
  return(current_layer)
})

# Re-combine and calculate the true individual ring areas
final_layer_cake <- do_call(rbind, hollow_layers) %>% 
  st_as_sf() %>% 
  mutate(area = as.numeric(st_area(geometry)))

# =================================================================
# STEP 5: Verification Print & Plot
# =================================================================
print(final_layer_cake %>% st_drop_geometry() %>% select(z_layer, area))

plot(final_layer_cake["z_layer"], 
     pal = terrain.colors(nrow(final_layer_cake)), 
     main = "True Layer-Cake Bands (Open Contours Accounted For)")




# mapview(raw_data)
# mapview(contours_simple)

mapview(final_layer_cake, z = "z_layer")

mapview(cake_layers, z = "z_layer")

mapview( arrange(desc(solid_layers)), z = "z_layer")


mapview( 
  filter(solid_layers, z_layer  < 10)
  )

flc_areas <- solid_layers %>%
  mutate(
    # Calculate area and strip the 'units' attribute for a clean numeric column
    area = as.numeric(st_area(geometry))
  )
