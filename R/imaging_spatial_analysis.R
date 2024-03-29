get_coor_from_string <- function(string){
  temp <- gsub('^.|.$', '', string) # remove the first and last string
  temp <- strsplit(temp, split=', ')
  return(c(as.numeric(temp[[1]][1]), as.numeric(temp[[1]][2])))
}

get_count_from_xy_coor <- function(file, topleft, bottomright, image=NULL, plot=FALSE){
  
  res <- read.table(file, header=TRUE, sep=',')
  xmin <- topleft[1]; ymin <- topleft[2]
  xmax <- bottomright[1]; ymax <- bottomright[2]

  cell_count <- c()
  for(i in 1:nrow(res)){
    coor_temp <- res[i, "xy_coordinate"]
    coor_temp <- get_coor_from_string(string=coor_temp)
    x_coor <- coor_temp[1]
    y_coor <- coor_temp[2]
    x_inside <- y_inside <- FALSE
    if((xmin < x_coor) & (x_coor < xmax)){x_inside <- TRUE}
    if((ymin < y_coor) & (y_coor < ymax)){y_inside <- TRUE}
    if(x_inside & y_inside){
      cell_count[i] <- TRUE
    } else {
      cell_count[i] <- FALSE
    }
  }
  
  if(plot==TRUE){
    if(is.null(image)) {stop("Image should be input to plotting")}
    img <- magick::image_read(image)
    info <- magick::image_info(img)
    plot(img)
    ymin_plot <- info$height-ymin
    ymax_plot <- info$height-ymax
    rect(xleft=xmin, xright=xmax, ybottom=ymin_plot, ytop=ymax_plot, border='yellow')
  }
  
  return(sum(cell_count))
}

#get_count_from_xy_coor(file='/Users/shan/Desktop/M926910_Position1_CD3-BUV395_sizes_coordinates.txt', topleft=c(500,0), bottomright=c(1392,500), image='/Users/shan/Desktop/M926910_Position1_CD3-BUV395.tiff', plot=TRUE)
