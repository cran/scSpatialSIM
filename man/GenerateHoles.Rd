% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_holes.R
\name{GenerateHoles}
\alias{GenerateHoles}
\title{Generate holes in a spatial simulation object}
\usage{
GenerateHoles(
  sim_object,
  xmin = NA,
  xmax = NA,
  ymin = NA,
  ymax = NA,
  sdmin = 1/2,
  sdmax = 2,
  hole_prob = c(0.2, 0.35),
  force = FALSE,
  density_heatmap = FALSE,
  step_size = 1,
  cores = 1,
  overwrite = FALSE,
  use_window = FALSE
)
}
\arguments{
\item{sim_object}{A spatial simulation object of class \code{SpatSimObj}}

\item{xmin}{Minimum x-coordinate for the holes (default: NA)}

\item{xmax}{Maximum x-coordinate for the holes (default: NA)}

\item{ymin}{Minimum y-coordinate for the holes (default: NA)}

\item{ymax}{Maximum y-coordinate for the holes (default: NA)}

\item{sdmin}{Minimum standard deviation for the kernels (default: 1/2)}

\item{sdmax}{Maximum standard deviation for the kernels (default: 2)}

\item{hole_prob}{A vector of length 2 with the minimum and maximum probabilities of
a point being within a hole (default: c(0.2, 0.35))}

\item{force}{Logical; if TRUE, forces the function to simulate outside the window
boundaries (default: FALSE)}

\item{density_heatmap}{Logical; if TRUE, computes a density heatmap (default: FALSE)}

\item{step_size}{The step size for the grid (default: 1)}

\item{cores}{The number of cores to use for parallel processing (default: 1)}

\item{overwrite}{boolean to replace holes if they have been simulated previously}

\item{use_window}{boolean whether to use the simulation window to set x and y limits}
}
\value{
A modified spatial simulation object with holes added
}
\description{
This function generates holes (regions of low probability) in a spatial simulation
object based on user-defined parameters. The function uses a kernel density
estimate to simulate holes, and returns a modified version of the input object
with the holes added. The function also has options to compute a density heatmap
and to assign points within the holes to be dropped or kept based on a scaled
probability value.
}
\details{
The function first checks that the input object is of the correct class,
and that no parameters are NULL. If any parameters are NULL, the function stops
with an error message. If the x- and y-ranges for the holes extend beyond the
boundaries of the simulation window, the function also stops with an error message,
unless the force parameter is set to TRUE. The function then produces kernel
parameter lists for each simulated pattern, and generates a grid based on the user-defined
step size. If density_heatmap is set to TRUE, the function computes a density heatmap
using the CalculateGrid function. Finally, the function computes hole probabilities
for each simulated pattern, assigns each point to be dropped or kept based on a
scaled probability value, and returns the modified simulation object.
}
\examples{
sim_object <- CreateSimulationObject()

#simulate points
sim_object <- GenerateSpatialPattern(sim_object, lambda = 20)

# Generate tissue with default parameters
sim_object <- GenerateTissue(sim_object)

# Generate holes in the simulation object
sim_object <- GenerateHoles(sim_object, hole_prob = c(0.1, 0.3), force = TRUE)

}
