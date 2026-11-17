#' Estimate the Mean by Iterative Interval Shrinking
#'
#' Computes the sample mean by numerically minimizing the squared-error loss
#' function using an iterative interval-shrinking search algorithm. The
#' algorithm starts with an interval bounded by the minimum and maximum
#' observations and repeatedly reduces the interval until convergence.
#'
#' @param data A numeric vector containing the observations.
#' @param tol A positive numeric tolerance controlling the stopping criterion.
#'   The algorithm stops when the absolute change in the objective function
#'   value between iterations is less than or equal to `tol`.
#'
#' @return An object of class `"my_mean"` containing:
#' \describe{
#'   \item{sol}{The estimated sample mean.}
#'   \item{data}{The original data vector.}
#'   \item{iter}{The number of iterations required for convergence.}
#'   \item{tol}{The tolerance used by the algorithm.}
#' }
#'
#' @details
#' The function minimizes
#' \deqn{
#' \sum_{i=1}^n (x_i - \mu)^2
#' }
#' using a one-dimensional interval-shrinking search. Because the squared-error
#' loss is convex, the procedure converges to the sample mean.
#'
#' @examples
#' x <- c(2, 4, 5, 7, 12)
#'
#' fit <- my_mean(x)
#'
#' fit$sol
#' mean(x)
#'
#' @export
my_mean <- function(
    data,
    tol = 10^(-12)
) {

  ###########################################################################
  # CHECK THE INPUT
  ###########################################################################

  # The function is designed for numerical data.

  if (!is.numeric(data)) {

    stop(
      "data must be numeric."
    )

  }


  # We also need at least one observation.

  if (length(data) == 0) {

    stop(
      "data must contain at least one observation."
    )

  }


  # This simple version does not allow missing observations.

  if (any(is.na(data))) {

    stop(
      "data cannot contain missing values."
    )

  }


  # The tolerance should be positive.

  if (tol <= 0) {

    stop(
      "tol must be positive."
    )

  }



  ###########################################################################
  # STEP 1: INITIAL SEARCH INTERVAL
  ###########################################################################

  # The sample mean must lie between:
  #
  #   the minimum observation
  #
  # and
  #
  #   the maximum observation.
  #
  # We therefore use these values as the endpoints of our
  # initial search interval.



  # Left endpoint:

  mu1 <- min(data)


  # Right endpoint:

  mu2 <- max(data)


  # Our initial interval is:
  #
  #   [mu1, mu2]



  ###########################################################################
  # STEP 2: INITIAL MIDPOINT
  ###########################################################################

  # Start at the midpoint of the interval.

  mid_mu <- (mu1 + mu2) / 2


  # Calculate the squared-error loss at this midpoint.

  mid_loss <- sum(
    (data - mid_mu)^2
  )



  ###########################################################################
  # STEP 3: INITIALIZE THE STOPPING RULE
  ###########################################################################

  # At each iteration, we will compare:
  #
  #   the loss at the OLD midpoint
  #
  # with:
  #
  #   the loss at the NEW midpoint.
  #
  #
  # Define:
  #
  #   diff_loss =
  #
  #       |old loss - new loss|
  #
  #
  # The algorithm stops when:
  #
  #   diff_loss <= tol



  # Start with an infinite value so that the while loop
  # is guaranteed to begin.

  diff_loss <- Inf



  ###########################################################################
  # STEP 4: INITIALIZE AN ITERATION COUNTER
  ###########################################################################

  # We will also count how many times the algorithm updates
  # the search interval.

  iteration <- 0



  ###########################################################################
  # STEP 5: ITERATIVE INTERVAL-SHRINKING SEARCH
  ###########################################################################

  while (diff_loss > tol) {


    # Increase the iteration counter.

    iteration <- iteration + 1



    # -----------------------------------------------------------------------
    # CREATE TWO NEW CANDIDATE POINTS
    # -----------------------------------------------------------------------

    # left_mid_mu lies halfway between:
    #
    #   mu1 and mid_mu

    left_mid_mu <- (mu1 + mid_mu) / 2


    # right_mid_mu lies halfway between:
    #
    #   mid_mu and mu2

    right_mid_mu <- (mid_mu + mu2) / 2



    # -----------------------------------------------------------------------
    # CALCULATE THE LOSS AT EACH CANDIDATE
    # -----------------------------------------------------------------------

    left_mid_loss <- sum(
      (data - left_mid_mu)^2
    )


    right_mid_loss <- sum(
      (data - right_mid_mu)^2
    )



    # -----------------------------------------------------------------------
    # DECIDE WHICH PART OF THE INTERVAL TO DISCARD
    # -----------------------------------------------------------------------

    # The squared-error loss is CONVEX.
    #
    # In one dimension, it has a U shape.
    #
    # If:
    #
    #   left_mid_loss < right_mid_loss
    #
    # the left candidate gives a better value of the objective.
    #
    # Therefore, we shrink the interval from the RIGHT.


    if (left_mid_loss < right_mid_loss) {

      mu2 <- right_mid_mu


    } else {

      # Otherwise, the right-side candidate is at least as good,
      # so shrink the interval from the LEFT.

      mu1 <- left_mid_mu

    }



    # -----------------------------------------------------------------------
    # CALCULATE THE NEW MIDPOINT
    # -----------------------------------------------------------------------

    mid_mu <- (mu1 + mu2) / 2



    # -----------------------------------------------------------------------
    # CALCULATE THE LOSS AT THE NEW MIDPOINT
    # -----------------------------------------------------------------------

    new_mid_loss <- sum(
      (data - mid_mu)^2
    )



    # -----------------------------------------------------------------------
    # CHECK FOR CONVERGENCE
    # -----------------------------------------------------------------------

    diff_loss <- abs(
      mid_loss - new_mid_loss
    )



    # -----------------------------------------------------------------------
    # UPDATE THE REFERENCE LOSS
    # -----------------------------------------------------------------------

    # The current loss becomes the "old" loss for the next iteration.

    mid_loss <- new_mid_loss

  }



  ###########################################################################
  # STEP 6: CREATE THE OUTPUT
  ###########################################################################

  # We want to return MORE than just the estimated mean.
  #
  # We also want to preserve:
  #
  #   - the original data;
  #   - the number of iterations;
  #   - the numerical tolerance.
  #
  # Since these are several different objects, a LIST is convenient.


  output <- list(

    # Numerical solution:

    sol = mid_mu,


    # Original observations:

    data = data,


    # Number of iterations:

    iter = iteration,


    # Tolerance used:

    tol = tol

  )



  ###########################################################################
  # STEP 7: GIVE THE OUTPUT A CUSTOM CLASS
  ###########################################################################

  # output is fundamentally still a LIST.
  #
  # However, we also assign it the custom class:
  #
  #   "my_mean"
  #
  # This will allow generic functions such as plot() to recognize
  # our object and use methods that we define specifically for it.

  class(output) <- c(
    "my_mean",
    "list"
  )



  ###########################################################################
  # STEP 8: RETURN THE OBJECT
  ###########################################################################

  # Since output is the final expression in the function, R
  # automatically returns it.
  #
  # We could equivalently write:
  #
  #   return(output)

  output

}


#' Plot the Squared-Error Loss Function for a my_mean Object
#'
#' Produces a visualization of the squared-error loss function
#' associated with a fitted `"my_mean"` object. The plot shows the
#' objective function over a grid of candidate values for \eqn{\mu},
#' together with vertical reference lines indicating the numerical
#' solution returned by the optimization algorithm and the exact
#' sample mean.
#'
#' @param x An object of class `"my_mean"` returned by [my_mean()].
#' @param ... Additional graphical arguments passed to [plot()].
#'
#' @return
#' Called for its side effect of producing a plot. Invisibly returns
#' `NULL`.
#'
#' @details
#' The function evaluates the squared-error loss
#' \deqn{
#' \sum_{i=1}^n (x_i - \mu)^2
#' }
#' over a fine grid of candidate values of \eqn{\mu}. The resulting
#' curve illustrates the convex nature of the objective function.
#'
#' A dashed vertical line identifies the numerical solution obtained
#' by [my_mean()], while a second reference line identifies the exact
#' sample mean computed with [mean()]. Since the sample mean minimizes
#' the squared-error loss, the two lines should coincide up to
#' numerical tolerance.
#'
#' @examples
#' x <- c(2, 4, 5, 7, 12)
#'
#' fit <- my_mean(x)
#'
#' plot(fit)
#'
#' # Pass additional graphical arguments
#' plot(fit, col = "steelblue", lwd = 2)
#'
#' @seealso [my_mean()]
#'
#' @export
plot.my_mean <- function(x, ...) {

  ###########################################################################
  # EXTRACT INFORMATION FROM THE OBJECT
  ###########################################################################

  # x is the object returned by my_mean().
  #
  # Therefore:
  #
  #   x$data
  #
  # contains the original observations, and:
  #
  #   x$sol
  #
  # contains our numerical solution.



  ###########################################################################
  # CREATE A GRID OF POSSIBLE mu VALUES
  ###########################################################################

  # To visualize the loss function, construct a fine sequence from
  # the minimum observation to the maximum observation.

  mu_grid <- seq(
    from = min(x$data),
    to = max(x$data),
    length.out = 500
  )



  ###########################################################################
  # CREATE STORAGE FOR THE LOSS
  ###########################################################################

  # We need one loss value for every possible value in mu_grid.

  loss_grid <- rep(
    NA,
    length(mu_grid)
  )



  ###########################################################################
  # CALCULATE THE LOSS OVER THE GRID
  ###########################################################################

  for (i in seq_along(mu_grid)) {

    loss_grid[i] <- sum(
      (x$data - mu_grid[i])^2
    )

  }



  ###########################################################################
  # PLOT THE LOSS FUNCTION
  ###########################################################################

  plot(
    mu_grid,
    loss_grid,
    type = "l",
    xlab = expression(mu),
    ylab = "Squared-error loss",
    main = "Optimization of the Squared-Error Loss",
    ...
  )


  # Notice the use of:
  #
  #   ...
  #
  # This allows the user to pass additional graphical arguments
  # to plot().
  #
  # For example:
  #
  #   plot(result, lwd = 2)



  ###########################################################################
  # MARK THE NUMERICAL SOLUTION
  ###########################################################################

  abline(
    v = x$sol,
    lty = 2
  )



  ###########################################################################
  # MARK THE EXACT SAMPLE MEAN
  ###########################################################################

  abline(
    v = mean(x$data),
    lty = 3
  )



  ###########################################################################
  # ADD A LEGEND
  ###########################################################################

  legend(
    "topright",
    legend = c(
      "Numerical solution",
      "Sample mean"
    ),
    lty = c(
      2,
      3
    ),
    bty = "n"
  )

}
