get_coor_from_string <- function(string){
  temp <- gsub('^.|.$', '', string) # remove the first and last string
  temp <- strsplit(temp, split=', ')
  return(c(as.numeric(temp[[1]][1]), as.numeric(temp[[1]][2])))
}

get_count_from_xy_coor <- function(file, topleft, bottomright){
  
  res <- read.table(file, header=T, sep=',')
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
  
  sum(cell_count)
}

#get_count_from_xy_coor(file='/Users/shan/Desktop/M926910_Position1_CD3-BUV395_sizes_coordinates.txt', topleft=c(0,0), bottomright=c(1392,600))
