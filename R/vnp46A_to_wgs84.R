

#'vnp46A_to_wgs84
#'@description project VNP46A black marble products to WGS84. Input can be multi-band. Can use in conjunction with
#'[purrr::map] to reproject all bands of all rasters contained in folder
#'@param path \code{character} path to hd5 black marble product (VNP46A1 or VNP46A2)
#'@return VNP46A raster/stack in WGS84
#'@examples \dontrun{
#'    library(terra)
#'    vnp46a1_wgs84 <- vnp46A_to_wgs84(path = "../../raster/VNP46A1.A2022043.h21v04.001.2022044083115.h5")
#'    # DNB band is 5th band so we can just write that out or we could combine bands however we want
#'    vnp46a1_dnb <- vnp46a1_wgs84[[5]]
#'    writeRaster(vnp46a1_dnb,"vnp46a1_dnb.tiff")
#'}

vnp46A_to_wgs84 <-  function(path){

  terra_ras <- suppressWarnings(terra::rast(path))
  wgs84_extent <- get_vnp46A_extent(path)
  terra::ext(terra_ras) <- wgs84_extent
  terra::crs(terra_ras) <- "epsg:4326"
  cat(crayon::green("returning terra SpatRaster with CRS = WGS84\n"))
  return(terra_ras)
}


#' get_vnp46A_extent
#' @param path \code{character} path to hd5 black marble product (VNP46A1 or VNP46A2)
#' @return A \code{numeric} vector containing new extent in correct extent format for `SpatRaster` class ([terra] package)

get_vnp46A_extent <- function(path){
  file_meta<- terra::describe(path,meta=T,parse=T)
  tile_number_meta<- stringr::str_subset(file_meta,"^HorizontalTileNumber.+|^VerticalTileNumber.+")
  tile_numbers<- readr::parse_number(tile_number_meta)
  horizontal_tile_number <- tile_numbers[1]
  vertical_tile_number <- tile_numbers[2]
  west_bound_coord <-  (10*horizontal_tile_number) - 180
  north_bound_coord <-  90-(10*vertical_tile_number)
  east_bound_coord <-  west_bound_coord + 10
  south_bound_coord = north_bound_coord - 10

  return(c(west_bound_coord,east_bound_coord,south_bound_coord,north_bound_coord))
}




#' group_vnp46a_images_by_doy
#'
#' @param ras_list named list of `SpatRaster`s. Names must follow NASA VNP46A file naming convention
#'
#' @return nested list where each level contains all images for doy.
#' @export
#'
#' @examples \dontrun{
#'  library(terra)
#'
#'  h5_filenames <- list.files("raster/h5/")
#'  h5_paths <- list.files("raster/h5/", full.names = T)

#'  black_marble <- h5_paths |>
#'     map(~vnp46A_to_wgs84(path =.x)) |>
#'      set_names(h5_filenames)
#'
#'  black_marble_by_doy<- group_vnp46a_images_by_doy(vnp46a1_wgs84)
#'
#'
#' }
group_vnp46a_images_by_doy<- function(ras_list){
  ras_names<- names(ras_list)
  doys<- stringr::str_split(ras_names,"\\.") |>
    map(2) |>
    map(~str_sub(.x,start = 6,end = 8)) |>
    unlist()
  ras_name_doy_lookup <- data.frame(ras_names,doys)
  u_doys <- unique(ras_name_doy_lookup$doys)
  doy_group_list <-  list()
  for(i in 1:length(u_doys)){
    doy_temp <- u_doys[i]
    ras_name_temp<- ras_name_doy_lookup |>
      filter(doys==doy_temp) |>
      pull(ras_names)

    doy_group_list[[doy_temp]] <- terra::sprc(ras_list[ras_name_temp])

  }
  return(doy_group_list)


}

#' mosaic_by_doy
#'
#' @param ras_list named list of `SpatRaster`s. Names must follow NASA VNP46A file naming convention
#'
#' @return nested list where each level contains all images for doy.
#' @export
#'
#' @examples \dontrun{
#'  library(terra)
#'
#'  h5_filenames <- list.files("raster/h5/")
#'  h5_paths <- list.files("raster/h5/", full.names = T)

#'  black_marble <- h5_paths |>
#'     map(~vnp46A_to_wgs84(path =.x)) |>
#'      set_names(h5_filenames)
#'
#'  black_marble_by_doy<- group_vnp46a_images_by_doy(vnp46a1_wgs84)
#'
#'
#' }
mosaic_vnp46a_by_doy<- function(ras_list,output=NULL){
  ras_names<- names(ras_list)
  doys<- stringr::str_split(ras_names,"\\.") |>
    map(2) |>
    map(~str_sub(.x,start = 6,end = 8)) |>
    unlist()
  ras_name_doy_lookup <- data.frame(ras_names,doys)
  u_doys <- unique(ras_name_doy_lookup$doys)
  doy_group_list <-  list()
  for(i in 1:length(u_doys)){
    doy_temp <- u_doys[i]
    ras_name_temp<- ras_name_doy_lookup |>
      filter(doys==doy_temp) |>
      pull(ras_names)

    terra_collection <- terra::sprc(ras_list[ras_name_temp])
    image_mosaic <- terra::mosaic(terra_collection)
    doy_group_list[[doy_temp]] <- image_mosaic
    if(!is.null(output)){
      cat(crayon::green("writing output to ",output,"\n"))
      cal_date <- as.Date(as.numeric(doy_temp), origin = '2021-12-31')|> stringr::str_remove_all(pattern = "-")

      ras_name_temp |> length()
      ras_name_temp[1]
      file_name_prefix <- ras_name_temp |>
        stringr::str_split(pattern = "\\.") |>
        map(1) |> unique() |> unlist()
      output_path_name <-  paste0(output,"/",file_name_prefix,".",cal_date,".wgs84.tiff")
      cat(crayon::green("writing ", output_path_name,"\n"))
      terra::writeRaster(image_mosaic,output_path_name)
    }


  }
  cat(crayon::green("returning SpatRaster mosaics as list"))
  return(doy_group_list)



}
# vnp_doy_paths_to_mosaic_paths <-  function
