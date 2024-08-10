# author: Sergio

# 


# Define a custom labeling function
my_labeller <- function(variable, value) {
  # Define custom labels for each facet
  custom_labels <- c(
    ws100 = "Wind Speed at 100m",
    ws10 = "Wind Speed at 10m"
    # Add more custom labels if needed
  )
  
  # Return the custom label for the given facet
  return(custom_labels[value])
}



convert_era5_time <- function(era5.time){
  # Convert time variable from nc file to POSIXct
  epoch <- as.POSIXct("1900-01-01 00:00:00", tz = "UTC")
  
  #convert to numeric. This represents hours after epoch
  time_numeric <- as.numeric(era5.time)
  
  # updating. This needs to convert time to seconds and then adding to reference point.
  time_posixct <- epoch + (time_numeric * 3600)
  
  return(time_posixct)
}


inter_grid_to_points <- function(loc.df, spatiotemp.array, lon, lat, time){
  # bilinear interpolation
  interp_series <- apply(spatiotemp.array, MARGIN = 3, 
                         FUN = \(wind.array) {
                           control_dat <- list(x=lon, y=lat, z=wind.array)
                           interp_var <- interp.surface(
                             control_dat,
                             desired.loc %>% dplyr::select(Longitude,Latitude))
                         },
                         simplify = T)
  
  # converting time
  time_posixct <- convert_era5_time(time)
  
  # merging everything in a single long data frame
  era.interp <- cbind(desired.loc,interp_series) %>%
    pivot_longer(!c(src_id:Latitude),values_to = "eraws10") %>% 
    mutate(name=as.numeric(name)) %>% 
    left_join(
      data.frame(name=1:length(time),time=time_posixct),
      by = "name"
    ) %>% mutate(time=as.character(time))
  
  return(era.interp)
}




# nearest_ind <- find_nearest_gridpoint(desired.loc, grid_points)
inter_grid_to_points <- function(loc.df, spatiotemp.array, lon, lat, time){
  
  # find nearest point and extract series
  temp <- mapply(
    FUN = \(x,y) {
      # get nearest lon and lat index
      lon_index <- which.min(abs(x - lon))
      lat_index <- which.min(abs(y - lat))
      # extract values
      variable_data[["ws10"]][lon_index, lat_index, ]
    },
    desired.loc$Longitude,
    desired.loc$Latitude,
    SIMPLIFY = F
  )
  
  # Convert time to POSIXct
  time_posixct <- convert_era5_time(time)
  
  # Create a data frame with the extracted value and time
  era.interp <- data.frame(
    Longitude = lon,
    Latitude = lat,
    eraws10 = nearest_value,
    time = as.character(time_posixct)
  )
  
  return(era.interp)
}


get_nearest_era5 <- function(loc_df, spatiotemp.array, lon, lat, time){
  
  # extract series of nearest point to each item in loc_df
  nearest.series <- mapply(
    FUN = \(x,y) {
      # get nearest lon and lat index
      lon_index <- which.min(abs(x - lon))
      lat_index <- which.min(abs(y - lat))
      # extract values
      spatiotemp.array[lon_index, lat_index, ]
    },
    desired.loc$Longitude,
    desired.loc$Latitude,
    SIMPLIFY = F
  )
  
  # gather result in a single table
  nearest.series <- nearest.series %>% 
    do.call(cbind,.) %>%
    as.data.frame() %>% 
    setNames(desired.loc$src_id) %>% 
    mutate(OB_TIME=time) %>% 
    pivot_longer(cols = -OB_TIME, names_to = "SRC_ID", values_to = "eraws10")
  
  nearest.series
}



myposixct <- function(timevar){
  as.POSIXct(
    case_when(
      nchar(timevar) == 10 ~ paste0(timevar, " 00:00:00"),
      nchar(timevar) == 16 ~ paste0(timevar, ":00"),
      TRUE ~ as.character(timevar)),
    format = "%Y-%m-%d %H:%M:%S", tz="UTC")
}


one.station.scatter <- function(
    data,
    station.id = NULL, station.name = NULL,
    station.table= station.table.inscope,
    p.col = "darkblue", p.alpha = 0.2, l.col = "red",
    ...){
  
  dat.yr = year(data$OB_TIME[1])
  if(!is.null(station.id)){
    stat.pos <- which(station.table$src_id==station.id)
    # get name
    stat.name <- station.table$Name[stat.pos]
    # get Area
    stat.area <- station.table$Area[stat.pos]
    # create subtitle
    my.sub <- sprintf("Area: %s, year: %04d",stat.area, dat.yr)
    
    # One station
    fig <- data %>% 
      filter(SRC_ID==station.id) %>% 
      # slice(sample(1:nrow(midas_clean),size = 100000,replace = F)) %>% 
      ggplot(aes(eraws10,WIND_SPEED))+
      geom_point(col=p.col,alpha=p.alpha)+
      geom_abline(slope = 1, intercept = 0, color = l.col) +  # Add 45-degree line
      # geom_smooth(method = "gam",formula = y ~s(x, bs ="cs"))+
      xlab("ERA5")+ylab("MIDAS")+labs(title = stat.name, subtitle = my.sub)+
      coord_fixed(ratio=1)+
      theme_bw()
  }
  print(fig)
}


one.station.ts <- function(
    data,
    station.id = NULL, station.name = NULL,
    station.table = station.table.inscope,
    roll.window.hr = 12,
    mypal2 = ggsci::pal_lancet()(2),
    group.var = FALSE,
    pal2col = ggsci::pal_lancet()(2),
    ...){
  
  dat.yr = year(data$OB_TIME[1])
  if(!is.null(station.id)){
    stat.pos <- which(station.table$src_id==station.id)
    # get name
    stat.name <- station.table$Name[stat.pos]
    # get Area
    stat.area <- station.table$Area[stat.pos]
    # create subtitle
    my.sub <- sprintf("Area: %s",stat.area)
    
    fig <- midas_clean %>% 
      filter(SRC_ID==station.id) %>% 
      mutate(trimestre = lubridate::quarter(OB_TIME)) %>% 
      {if (group.var) group_by(., trimestre) else .} %>%  # Conditionally group_by
      mutate(across(c(eraws10,WIND_SPEED),~ zoo::rollmean(., roll.window.hr, fill = c("extend", NA, "extend")))) %>% 
      ungroup() %>% 
      pivot_longer(cols = c(eraws10,WIND_SPEED)) %>% 
      ggplot(aes(OB_TIME,value, col=name))+
      geom_line()+
      {if (group.var) facet_wrap(
        ~trimestre, scales = "free",
        labeller=labeller(trimestre = \(x) sprintf("%04dQ%s",dat.yr,x))) } +
      xlab("time")+ylab("Wind speed (m/s)")+
      labs(title = sprintf("%s. %02d Hours Rolling average",stat.name,roll.window.hr),
           subtitle = my.sub,
           color= "source")+
      theme(legend.position = "bottom")+
      scale_color_manual(values = c("eraws10" = mypal2[1], "WIND_SPEED" = mypal2[2]),
                         labels = c("eraws10" = "ERA5", "WIND_SPEED" = "MIDAS"))
  }
  print(fig)
}

one.time.scatter <- function(
    data,
    one.time,
    plotly.convert = FALSE,
    p.col = "darkblue", l.col = "red",
    ...){
  # Create the ggplot object
  gg <- data %>% 
    rename(midas = WIND_SPEED, era5=eraws10) %>% 
    filter(OB_TIME == one.time) %>% 
    ggplot(aes(era5, midas,
               text = paste(Name, " - ", Area))) +
    geom_point(col = p.col) +
    geom_abline(slope = 1, intercept = 0, color = l.col) +
    xlab("ERA5") + ylab("MIDAS") +
    labs(title = paste0(format(one.time, "%d %b %y at %H:%M"))) +
    # coord_fixed(ratio = 1) +
    theme_bw()
  
  if(plotly.convert){
    ggplotly(gg) 
  }else{
    print(gg)
  }
  
}


# Define the function
process_data <- function(data, time_func) {
  data %>%
    mutate(tgroup = time_func(OB_TIME)) %>%
    group_by(SRC_ID, Name,Area, tgroup) %>%
    summarise(
      abs.bias = mean(abs(bias)),
      across(c(WIND_SPEED,WIND_DIRECTION, eraws10, bias), ~mean(., na.rm = TRUE)),
      .groups = "drop"
    )
}



# centroids function

get_polygon <- function(data, x.metric = "WIND_SPEED",  y.metric = "abs.bias",...){
  
  # get convex hull
  centroid_convex_hull <- data %>%
    st_as_sf(coords = c(x.metric, y.metric)) %>%
    st_set_crs(4326) %>%
    group_by(SRC_ID) %>%
    summarise(geometry = st_union(geometry)) %>%
    mutate(centroid = st_centroid(geometry),
           hull = st_convex_hull(geometry))
  
  
  centroid_convex_hull
}


plot_mean_polygons <- function(data, time_func = month,...){
  
  # group data and get means
  result <- process_data(data, time_func) %>% 
    left_join(stations.welev %>% 
                dplyr::select(src_id, elevation),
              by = c("SRC_ID" = "src_id"))
  
  # get polygons
  convex_hull <- get_polygon(result,...)
  
  # Plot
  
  fig <- ggplot() +
    geom_sf(data = convex_hull, aes(geometry = hull), fill = "#B19CD9", alpha = 0.5) +
    geom_sf(data = convex_hull, aes(geometry = centroid), color = "darkblue", size = 2) +
    labs(title = "Mean absolute bias by station. 2023",...) +
    scale_x_continuous(labels = scales::number_format())+
    scale_y_continuous(labels = scales::number_format())+
    xlab("MIDAS wind speed (m/s)") + ylab("Mean absolute bias")
  print(fig)
  invisible(fig)
}





my_mesh <- function(
    points, 
    n.points = nrow(points)*3,
    max.edge = c(0.9,5),
    cutoff = 0.1,
    min.angle = 21,
    figure = TRUE,
    ...){
  # Generate the INLA mesh
  loc.mesh <-  fmesher::fm_mesh_2d_inla(
    points,
    max.edge = max.edge, 
    cutoff = cutoff,
    min.angle = min.angle,
    max.n.strict = c(n.points,60),
    # plot.delay = TRUE,
    offset = c(-.05,-.05)
  )
  if(figure) plot(loc.mesh)
  
  invisible(loc.mesh)
}

my_mesh <- function(
    points, 
    n.points = nrow(points)*3,
    figure = TRUE,
    ...){
  # Generate the INLA mesh
  loc.mesh <-  fmesher::fm_mesh_2d_inla(
    points,
    max.n.strict = c(n.points,60),
    # plot.delay = TRUE,
    offset = c(-.05,-.05),
    ...
  )
  if(figure) plot(loc.mesh)
  
  invisible(loc.mesh)
}


premodel <- function(
    mesh, points, data, prior.range = c(1, 0.5),
    prior.sigma = c(1, 0.5), alpha = 2, ...){
  # Creating A matrix (Observation/prediction matrix)
  loc.A <- inla.spde.make.A(mesh, loc = points)
  #Creating Matern SPDE model object with PC prior
  loc.spde = inla.spde2.pcmatern(mesh = mesh,
                                 prior.range = prior.range,
                                 prior.sigma = prior.sigma,
                                 alpha = alpha)
  #Generating the SPDE model index vector
  loc.w <- inla.spde.make.index('w', n.spde = loc.spde$n.spde)
  
  # projgrid <- inla.mesh.projector(loc.mesh, xlim = range(Locations[, 1]), 
  #   ylim = range(Locations[, 2]), dims = nxy)
  
  projgrid <- inla.mesh.projector(mesh, xlim = range(points[, 1]), 
                                  ylim = range(points[, 2]))
  
  #First we make the model matrix using the model formula,
  #but without response and intercept.
  X0 <- model.matrix( ~ 0 + WIND_SPEED + elevation, data = data)
  
  X <- as.data.frame(X0) # convert to a data frame.
  # Making the stack
  N <- nrow(data) #Saving the number of rows in the data
  
  stack <- inla.stack(
    # specify the response variable
    data = list(y = data$bias),
    # Vector of Multiplication factors for fixed effects
    A = list(1, 1, loc.A),
    #Specify the fixed and random effects
    effects = list(
      # specify the manual intercept!
      Intercept = rep(1, N),
      # attach the model matrix
      X = X,
      # attach the w
      w = loc.w) )
  
  list(loc.spde = loc.spde, stack = stack, projgrid=projgrid)
}


run.model <- function(loc.spde, stack, projgrid){
  inlamodel <- inla(y ~ 0 + Intercept + WIND_SPEED + elevation +
                      f(w, model = loc.spde),
                    family = "Gaussian",
                    data = inla.stack.data(stack),
                    control.compute = list(cpo=T,dic = T),
                    control.predictor = list(A = inla.stack.A(stack)))
  
  inlamodel
}

summary.info <- function(inlamodel){
  model.summary <- summary(inlamodel)
  
  print(model.summary)
  model.nlscpo=-sum(log(inlamodel$cpo$cpo)) # Negative Log Sum CPO
  cat("NLSCPO of INLA model 2:",model.nlscpo,"\n")
  dic = inlamodel$dic$dic
  cat("DIC of INLA model 2:",dic,"\n") # Deviance Information Criteria
  sd.res = sd(df.sample$bias-inlamodel$summary.fitted.values$mean[1:N])
  cat("Standard deviation of mean residuals for INLA model:",
      sd.res,"\n") # 
  
  rmse = sqrt( mean((df.sample$bias-inlamodel$summary.fitted.values$mean[1:N])^2))
  cat("RMSE for INLA model:",
      rmse,"\n")
  
  invisible(list(
    summary=model.summary,
    nlscpo=model.nlscpo,
    dic=dic,
    sd.res,
    rmse=rmse))
}

main.model <- function(points,n.points,data,...){
  # mesh
  loc.mesh <- my_mesh(points, n.points,...)
  
  # premodel
  inputs <- premodel(loc.mesh, points, data,...)
  
  # model
  model <- run.model(inputs$loc.spde,inputs$stack,inputs$projgrid)
  
  
  list(model=model, mesh=loc.mesh)
}



fractional.hours <- function(time){
  # Convert to POSIXlt to extract components
  time_lt <- as.POSIXlt(time)
  
  # Extract components
  hours <- time_lt$hour
  minutes <- time_lt$min
  seconds <- time_lt$sec
  # Calculate fractional hours
  fractional_hours <- hours + (minutes / 60) + (seconds / 3600)
  fractional_hours
}

fractional.months <- function(time){
  
  # Convert to POSIXlt to extract components
  time_lt <- as.POSIXlt(time)
  
  # Extract components
  year_start <- as.POSIXct(paste(format(time_lt, "%Y"), "-01-01 00:00:00", sep=""))
  
  # Calculate the difference in seconds since the start of the year
  time_difference <- difftime(time, year_start, units = "secs")
  
  # Get the number of days in the current year to account for leap years
  days_in_year <- ifelse(leap_year(year_start), 366, 365)
  
  # Calculate fractional months
  fractional_months <- as.numeric(time_difference) / (days_in_year / 12 * 24 * 3600)
  
  # Print result
  fractional_months
}

coords_in_km <- function(sf_data){
  
  # Extract the coordinates
  coords_meters <- st_coordinates(sf_data)
  
  # Convert meters to kilometers
  coords_kilometers <- coords_meters / 1000
  
  # Recreate the geometry with the converted coordinates
  # Note: st_sfc needs a list of geometries, so we apply st_point to each row
  geom_km <- st_sfc(
    lapply(1:nrow(coords_kilometers), 
           function(i) st_point(coords_kilometers[i, 1:2])), 
    crs = st_crs(sf_data))
  
  
  # Update the geometry in the sf object
  sf_data_utm_km <- sf_data
  st_geometry(sf_data_utm_km) <- geom_km
  
  sf_data_utm_km
}


coords_in_km2 <- function(sf_data){
  
  
  # Transform coordinates to kilometers
  sf_data_km <- st_geometry(sf_data) / 1000

# Create a new sf object with the transformed coordinates
  sf_data_km_sf <- st_set_geometry(sf_data, sf_data_km)
 
  sf_data_km_sf
}

book.plot.field <- function(field, mesh, projector, xlim, ylim, 
                            dims=c(300,300), poly, asp = 1, 
                            axes = FALSE, xlab = '', ylab = '', 
                            col = book.color.c(), ...){
  ## you can supply field as a matrix vector or like a named list with 'x', 'y' and 'z' as for image
  ## when field is a vector, it will project it using projector, assuming projector will create a matrix 
  ## when mesh is supplied and projector not, projector will be created and used to project field
  if (missing(mesh)) {
    if (missing(projector)) {
      if (missing(xlim) | missing(ylim)) {
        image.plot(field, asp = asp, axes = axes, 
                   xlab = xlab, ylab = ylab, col = col, ...)
      } else {
        image.plot(field, xlim = xlim, ylim = ylim, asp = asp, 
                   axes = axes, xlab = xlab, ylab = ylab, col = col, ...)
      }
    } else {
      if (missing(xlim)) xlim <- range(projector$x)
      if (missing(ylim)) ylim <- range(projector$y)
      field.proj <- inla.mesh.project(projector, field)
      image.plot(x = projector$x, y = projector$y, z = field.proj, 
                 asp=asp, axes=axes, xlab = xlab, ylab = ylab, 
                 col=col, xlim=xlim, ylim=ylim, ...)
    }
  } else {
    if (missing(xlim)) xlim <- range(mesh$loc[,1])
    if (missing(ylim)) ylim <- range(mesh$loc[,2])
    projector <- inla.mesh.projector(mesh, xlim = xlim,
                                     ylim = ylim, dims=dims)
    field.proj <- inla.mesh.project(projector, field)
    image.plot(x = projector$x, y = projector$y, z = field.proj, 
               asp=asp, axes=axes, xlab = xlab, ylab = ylab, col=col, ...)
  }
  if (!missing(poly)) 
    plot(poly, add = TRUE, col = 'grey')
}
book.color.c = function(n = 201) {
  return(viridis(n))
}
book.color.c2 = function(n = 201) {
  return(magma(n))
}

get_country_polygon <- function(
    country=c("ireland","united kingdom"),
    crs.code = 32629, # Transforming to UTM Zone 29N (EPSG:32629)
    units = "km"
){
  shape <- ne_countries(
    scale = "medium", 
    country = c("ireland","united kingdom"),
    returnclass = "sf") %>% 
    st_transform(., crs.code) 
  
  
  if (units == "km") shape <- shape[1] %>% coords_in_km2() 
  return(shape) # return country polygon
}


mygfplot <- function(
    model,
    mesh,
    obs.points,
    add.country = TRUE,
    plot.obs = FALSE,
    color_palette1 = colorRampPalette(brewer.pal(11, "RdBu"))(200),
    color_palette2 = book.color.c2(),
    stepsize = 4,
    ...
){
  
  # projector on a fine grid
  x.range <- diff(range(obs.points[, 1]))
  y.range <- diff(range(obs.points[, 2]))
  nxy <- round(c(x.range, y.range) / stepsize)
  
  pgrid0 <- inla.mesh.projector(
    loc.mesh, 
    xlim = range(Locations[, 1])*c(0.8,1.2), 
    ylim = range(Locations[, 2])*c(0.975,1.025),
    dims = nxy)
  
  # extract projections
  prd0.m <- inla.mesh.project(
    pgrid0,
    model$summary.random$w$mean)                          
  prd0.s <- inla.mesh.project(
    pgrid0,
    model$summary.random$w$sd)
  
  book.plot.field(
    list(x = pgrid0$x, y = pgrid0$y, z = prd0.m),
    col = color_palette1,
    # poly = loc.mesh,
    # axes = TRUE,
    ...
  )
  if(add.country) plot(get_country_polygon(), add=TRUE, col=NA)
  if(plot.obs) points(obs.points, pch=19, cex = 0.5)
  book.plot.field(
    list(x = pgrid0$x, y = pgrid0$y, z = prd0.s),
    col = color_palette2,
    # poly = loc.mesh,
    # axes = TRUE,
    ...
  )
  if(add.country) plot(get_country_polygon(), border = "lightgray", 
                       add=TRUE, col="NA")
  if(plot.obs) points(obs.points, pch=19, cex = 0.5,col="lightgray")
}

field_plot <- function(
    field.atribute,
    mesh,
    obs.points,
    add.country = TRUE,
    plot.obs = FALSE,
    color_palette1 = colorRampPalette(brewer.pal(11, "RdBu"))(200),
    color_palette2 = book.color.c2(),
    border.country = "black",
    point.col = "black",
    stepsize = 4,
    ...
){
  
  # projector on a fine grid
  x.range <- diff(range(obs.points[, 1]))
  y.range <- diff(range(obs.points[, 2]))
  nxy <- round(c(x.range, y.range) / stepsize)
  
  pgrid0 <- inla.mesh.projector(
    loc.mesh, 
    xlim = range(Locations[, 1])*c(0.8,1.2), 
    ylim = range(Locations[, 2])*c(0.975,1.025),
    dims = nxy)
  
  # extract projections
  prd0.m <- inla.mesh.project(
    pgrid0,
    field.atribute)                          
  
  book.plot.field(
    list(x = pgrid0$x, y = pgrid0$y, z = prd0.m),
    col = color_palette1,
    # poly = loc.mesh,
    # axes = TRUE,
    ...
  )
  if(add.country) plot(get_country_polygon(), add=TRUE, col=NA, 
                       border = border.country , ...)
  if(plot.obs) points(obs.points, pch=19, cex = 0.5,
                      col = point.col, ...)
  invisible(list(grid = pgrid0,nxy = nxy))
}