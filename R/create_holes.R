#' Generate holes in a spatial simulation object
#'
#' This function generates holes (regions of low probability) in a spatial simulation
#' object based on user-defined parameters. The function uses a kernel density
#' estimate to simulate holes, and returns a modified version of the input object
#' with the holes added. The function also has options to compute a density heatmap
#' and to assign points within the holes to be dropped or kept based on a scaled
#' probability value.
#'
#' @param sim_object A spatial simulation object of class `SpatSimObj`
#' @param xmin Minimum x-coordinate for the holes (default: NA)
#' @param xmax Maximum x-coordinate for the holes (default: NA)
#' @param ymin Minimum y-coordinate for the holes (default: NA)
#' @param ymax Maximum y-coordinate for the holes (default: NA)
#' @param sdmin Minimum standard deviation for the kernels (default: 1/2)
#' @param sdmax Maximum standard deviation for the kernels (default: 2)
#' @param hole_prob A vector of length 2 with the minimum and maximum probabilities of
#'   a point being within a hole (default: c(0.2, 0.35))
#' @param force Logical; if TRUE, forces the function to simulate outside the window
#'   boundaries (default: FALSE)
#' @param density_heatmap Logical; if TRUE, computes a density heatmap (default: FALSE)
#' @param step_size The step size for the grid (default: 1)
#' @param cores The number of cores to use for parallel processing (default: 1)
#' @param overwrite boolean to replace holes if they have been simulated previously
#' @param use_window boolean whether to use the simulation window to set x and y limits
#'
#' @details The function first checks that the input object is of the correct class,
#' and that no parameters are NULL. If any parameters are NULL, the function stops
#' with an error message. If the x- and y-ranges for the holes extend beyond the
#' boundaries of the simulation window, the function also stops with an error message,
#' unless the force parameter is set to TRUE. The function then produces kernel
#' parameter lists for each simulated pattern, and generates a grid based on the user-defined
#' step size. If density_heatmap is set to TRUE, the function computes a density heatmap
#' using the CalculateGrid function. Finally, the function computes hole probabilities
#' for each simulated pattern, assigns each point to be dropped or kept based on a
#' scaled probability value, and returns the modified simulation object.
#'
#' @return A modified spatial simulation object with holes added
#' @export
#'
#' @examples
#' sim_object <- CreateSimulationObject()
#'
#' #simulate points
#' sim_object <- GenerateSpatialPattern(sim_object, lambda = 20)
#'
#' # Generate tissue with default parameters
#' sim_object <- GenerateTissue(sim_object)
#'
#' # Generate holes in the simulation object
#' sim_object <- GenerateHoles(sim_object, hole_prob = c(0.1, 0.3), force = TRUE)
#'
GenerateHoles = function(sim_object, xmin = NA, xmax = NA, ymin = NA, ymax = NA,
                         sdmin = 1/2, sdmax = 2,
                         hole_prob = c(0.2, 0.35),
                         force = FALSE,
                         density_heatmap = FALSE, step_size = 1, cores = 1, overwrite = FALSE,
                         use_window = FALSE){
  if(!methods::is(sim_object, "SpatSimObj")) stop("`sim_object` must be of class 'SpatSimObj'")
  if(any(is.null(c(xmin, xmax, ymin, ymax, sdmin, sdmax, hole_prob)))) stop("Cannot have `NULL` parameters")

  if(!is.empty(sim_object@Holes, "Simulated Kernels") & overwrite == FALSE) stop("Already have hole kernels and `overwrite == FALSE`")
  if(!is.empty(sim_object@Holes, "Simulated Kernels") & overwrite == TRUE){
    message("Overwriting existing hole kernels")
    # #tissue
    # sim_object@Tissue@`Simulated Kernels` = list()
    # sim_object@Tissue@`Density Grids` = list()
    #holes
    message("Resetting Hole slots")
    sim_object@Holes@`Simulated Kernels` = list()
    sim_object@Holes@`Density Grids` = list()
    # #cells
    # for(i in seq(sim_object@Cells)){
    #   sim_object@Cells[[i]]@`Simulated Kernels` = list()
    #   sim_object@Cells[[i]]@`Density Grids` = list()
    # }
    #spatial files
    sim_object@`Spatial Files` = lapply(sim_object@`Spatial Files`, function(spat){
      spat %>%
        dplyr::select(-dplyr::contains("Hole"))
    })
    #letting know finished
    message("Reset...Continuing.")
  }

  #check if using window is TRUE
  if(use_window){
    xmin = sim_object@Window$xrange[1]
    xmax = sim_object@Window$xrange[2]
    ymin = sim_object@Window$yrange[1]
    ymax = sim_object@Window$yrange[2]
  }

  #create parameter vector
  params = list(xmin = xmin, xmax = xmax,
             ymin = ymin, ymax = ymax,
             sdmin = sdmin, sdmax = sdmax,
             hole_prob = hole_prob)
  #if no parameters are input then use the initialized
  params = mapply(replace_na, sim_object@Holes@Parameters, params, SIMPLIFY = FALSE)
  #update initialized parameters with custom input from user
  sim_object@Holes@Parameters <- params
  #get the window size
  win_limits = c(sim_object@Window$xrange, sim_object@Window$yrange)
  #check whether the parameters would simulate outside window
  if(any((unlist(params[c(1, 3)]) < win_limits[c(1,3)]) |
         (unlist(params[c(2, 4)]) > win_limits[c(2,4)])) & force == FALSE){
    stop("x and y range outside simulation window limits")
  }
  #inform user parameter window inside simulation window
  if(any(c(unlist(params[c(1, 3)]) > win_limits[c(1,3)],
           unlist(params[c(2, 4)]) < win_limits[c(2,4)]))){
    message("x and y range inside window boundary")
  }
  #produce kernel parameter list for k clusters in each simulated pattern
  sim_object@Holes@`Simulated Kernels` = lapply(seq(sim_object@Sims), function(hld){
    do.call(generate_holes, params)
  })
  #make gric based on user step size for their window
  grid = expand.grid(x = seq(win_limits[1], win_limits[2], step_size),
                     y = seq(win_limits[3], win_limits[4], step_size))

  if(density_heatmap){
    message("Computing density heatmap")
    sim_object@Holes@`Density Grids` = pbmcapply::pbmclapply(sim_object@Holes@`Simulated Kernels`, function(gauss_tab){
      cbind(grid,
            prob = CalculateGrid(grid, gauss_tab, cores = cores))
    }, mc.cores = 1)
  }

  if(is.empty(sim_object, "Spatial Files")){
    sim_object@`Spatial Files` = lapply(sim_object@Patterns, data.frame)
  }

  message("Computing hole probability")
  sim_object@`Spatial Files` = pbmcapply::pbmclapply(seq(sim_object@`Spatial Files`), function(spat_num){
    df = cbind(sim_object@`Spatial Files`[[spat_num]],
               #this is for making sharp hole edges
               # CalculateGridHoles(sim_object@`Spatial Files`[[spat_num]],
               #                                      sim_object@Holes@`Simulated Kernels`[[spat_num]], cores = cores)
               `Hole Probability` = CalculateGrid(sim_object@`Spatial Files`[[spat_num]],
                                                  sim_object@Holes@`Simulated Kernels`[[spat_num]], cores = cores)
               )
    #scale
    df$`Hole Probability Scaled` = sqrt(df$`Hole Probability`)
    # df$`Hole Assignment` = ifelse(df$`Hole Z` < sim_object@Holes@`Simulated Kernels`[[spat_num]][df$`Closest Hole`,"max_dist"],
    #                               "Drop", "Keep")
    df$`Hole Assignment` = ifelse(stats::rbinom(nrow(df), size = 1, prob = df$`Hole Probability Scaled`) == 1, "Drop", "Keep")
    return(df)
  }, mc.cores = 1)

  return(sim_object)
}
