###############################################################
### 0 - Set graphical parameters for output plots

ptype=".pdf"
pscale=0.9
pheight=11.25
pwidth=20
punits="cm"
pbg="transparent" #background color

## Set adjusted ggplot light theme
default_theme <- theme_get()
theme_light_adj <- theme_light()

# Set transparent background for printing on color
theme_light_adj$panel.background = element_rect(fill = "transparent") # bg of the panel
theme_light_adj$plot.background = element_rect(fill = "transparent", color = NA) #bg of the plot
theme_light_adj$legend.background = element_rect(fill = "transparent", color= NA) # get rid of legend bg
theme_light_adj$legend.box.background = element_rect(fill = "transparent", color= NA) # get rid of legend panel bg
# Set theme
theme_set(theme_light_adj)

# Define new ggplot color theme
# scales::show_col(tableau_color_pal('Classic 10')(10)) 
scale_colour_discrete <- function(...) scale_color_tableau(palette="Classic 10")
scale_fill_discrete <- function(...) scale_fill_tableau(palette="Classic 10")
#scale_colour_discrete <- function(...) scale_color_tableau(palette="Classic Color Blind")
#scale_fill_discrete <- function(...) scale_fill_tableau(palette="Classic Color Blind")
#scale_colour_gradient2 <- function(...) scale_colour_gradient2_tableau(palette="Classic Orange-White-Blue")
#scale_fill_gradient2 <- function(...) scale_fill_gradient2_tableau(palette="Classic Orange-White-Blue")
scale_colour_gradient2 <- function(...) scale_colour_gradient2_tableau(palette="Classic Orange-Blue")
scale_fill_gradient2 <- function(...) scale_fill_gradient2_tableau(palette="Classic Orange-Blue")


# Save colors for heatmaps
palettes <- ggthemes_data[["tableau"]][["color-palettes"]][["ordered-diverging"]]
divcolors <- palettes$`Classic Orange-White-Blue`$value
divcolors <- palettes$`Classic Orange-Blue`$value
#scales::show_col(divcolors)
tableau_div <- colorRampPalette(rev(divcolors), bias = 2)
#scales::show_col(tableau_div(100))