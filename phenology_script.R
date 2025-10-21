# 1. CHUNK OF CODE
# I will use the phenocamr package which 
# interfaces with the phenocam network API
# to download time series of vegetation 
# greenness and derived phenology metrics
library(phenocamr)

# download greenness time series,
# calculate phenology (phenophases),
# amend with DAYMET data
phenocamr::download_phenocam(
  site = "harvard$",
  veg_type = "DB",
  roi_id = "1000",
  daymet = TRUE,
  phenophase = TRUE,
  trim = 2022,
  out_dir = "data-raw/" #to indicate in which folder to save the downloaded data
)

harvard_phenocam_data <- readr::read_csv(         #reads the data into data frame
  file.path("data-raw/harvard_DB_1000_3day.csv"), #points to where to find that data to read in
  comment = "#"                                   #ignore lines that start w/ #
)

# reading in harvard phenology only retaining
# spring (rising) phenology for the GCC 90th
# percentile time series (the default)
harvard_phenology <- readr::read_csv(
  file.path("data-raw/harvard_DB_1000_3day_transition_dates.csv"),
  comment = "#"
) |>
  dplyr::filter(
    direction == "rising", #keeps only rising phenophase transitions aka only the spring data
    gcc_value == "gcc_90"  #keeps only the column type gcc_90 aka only for the 90th percentile ( value below which 90& of data fall)
  )

#CODE FOR FIGURE 6.1 (not part of what we need to correct I assume)
library(ggplot2)
ggplot(harvard_phenocam_data) +
  geom_line(
    aes(
      as.Date(date),
      smooth_gcc_90
    ),
    colour = "grey25"
  ) +
  geom_point(
    data = harvard_phenology,
    aes(
      as.Date(transition_25),
      threshold_25
    )
  ) +
  labs(
    x = "",
    y = "GCC"
  ) +
  theme_bw() +
  theme(
    legend.position = "none"
  )

#N 2. CHUNK OF CODE
# return mean daily temperature as well
# as formal dates (for plotting)
library(dplyr)
harvard_temp <- harvard_phenocam_data |> #|> puts output of one fct in next fct
  group_by(year) |> #to reset GDD calculation at beginning of each year
  dplyr::mutate(    #creates new column with mean daily temp
    tmean = (tmax..deg.c. + tmin..deg.c.)/2
  ) |> 
  dplyr::mutate(
    date = as.Date(date),
    gdd = cumsum(ifelse(tmean >= 5, tmean - 5, 0)) #calculates the daily GDD contribution using the formula from text
  ) |>
  dplyr::select( #selects / only keeps the columns you need for further analysis
    date,
    year,
    tmean,
    gdd
  ) |>
  ungroup() #removes the grouping setting from beginning 

# convert the harvard phenology data and only
# retain required data
harvard_phenology <- harvard_phenology |> #takes pheno dataset and pipes it forward
  mutate(
    doy = as.numeric(format(as.Date(transition_25),"%j")), #access the trans_25 vector w/ the leaf out date & returns it into day of the year in numeric format  
    year = as.numeric(format(as.Date(transition_25),"%Y")) #does same for year of this date -> later needed to plot GDD
  ) |>
  select( #keeps only the relevant columns
    year,
    doy,
    transition_25,
    threshold_25
  ) #now gives a table w/ transition dates (=leaf out dates) & their date of the year & year


#CODE TO PLOT FIG 6.2
library(patchwork) #to be able to combine ggplot's at the end
# grab only the 2010 value of spring phenology
harvard_phenology_2010 <- harvard_phenology |>
  dplyr::filter(
    year == 2010
  )

harvard_gdd_value <- harvard_temp |>
  dplyr::filter(
    date == harvard_phenology_2010$transition_25
  )

p <- ggplot(harvard_temp) +
  geom_line(
    aes(
      date,
      tmean
    )
  ) +
  geom_point(
    aes(
      date,
      tmean,
      colour = tmean > 5,
      group = 1
    )
  ) +
  geom_vline(
    data = harvard_phenology_2010,
    aes(
      xintercept = as.Date(transition_25)
    )
  ) +
  scale_colour_discrete(
    type = c(
      "blue",
      "red"
    )
  ) +
  labs(
    x = "",
    y = "Temperature (deg. C)"
  ) +
  xlim(
    c(
      as.Date("2010-01-01"),
      as.Date("2010-06-30")
    )
  ) +
  theme_bw() +
  theme(
    legend.position = "none"
  )

p2 <- ggplot(harvard_temp) +
  geom_line(
    aes(
      date,
      gdd
    )
  ) +
  geom_point(
    aes(
      date,
      gdd,
      colour = tmean > 5,
      group = 1
    )
  ) +
  scale_colour_discrete(
    type = c(
      "blue",
      "red"
    )
  ) +
  geom_vline(
    data = harvard_phenology_2010,
    aes(
      xintercept = as.Date(transition_25)
    )
  ) +
  geom_hline(
    data = harvard_gdd_value,
    aes(
      yintercept = gdd
    ),
    lty = 2
  ) +
  labs(
    x = "",
    y = "GDD (deg. C)"
  ) +
  xlim(
    c(
      as.Date("2010-01-01"),
      as.Date("2010-06-30")
    )
  ) +
  ylim(c(0, 1000)) +
  theme_bw()  +
  theme(
    legend.position = "none"
  )

# compositing
p + p2 + 
  plot_layout(ncol = 1) + 
  plot_annotation(
    tag_levels = "a",
    tag_prefix = "(",
    tag_suffix = ")"
  )

# 3. CHUNK OF CODE | 6.3
gdd_model <- function(temp, par) { #takes a temp series -> takes 2 parameters
  # split out parameters from a simple vector of parameter values
  temp_threshold <- par[1] #takes the first element of the parameter vector & assigns it to temp_threshold. represents base temp
  gdd_crit <- par[2] #takes second param & assigns it. critical value at which pheno event occ
  
  # accumulate growing degree days for temperature data
  gdd <- cumsum(ifelse(temp > temp_threshold, temp - temp_threshold, 0)) 
  
  # figure out when the number of growing degree days exceeds the minimum value 
  #required for leaf development, only return the first value
  doy <- unlist(which(gdd >= gdd_crit)[1]) #get 1st day when GDD reaches critical value
  
  return(doy)
}

# 4. CHUNK OF CODE
# confirm that the model function returns expected results (i.e. DOY 114)
# (we filter out the year 2010, but removing the filter would run the model for all years!)
prediction <- harvard_temp |>
  dplyr::filter(
    year == 2010
  ) |>
  group_by(year) |>
  summarize(
    pred = gdd_model(
      temp = tmean,
      par = c(5, 130.44)
    )  
  )

print(prediction)





