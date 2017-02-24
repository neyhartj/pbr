#' Plot AMMI Output
#'
#' @description
#' More flexible plotting function for the output of the \code{AMMI} function
#' in the \code{agricolae} package.
#'
#' @param AMMI.out The output from the \code{\link[AMMI]{agricolae}} function.
#' @param x The variable for the x-axis. Must be an \code{integer} specifying
#' the interaction principal component to plot. If \code{0}, the phenotypic
#' mean is plotted.
#' @param y The variable for the y-axis. Must be an \code{integer} specifying
#' the interaction principal component to plot.
#' @param show.gen \code{Logical}. Should the genotypes be highlighted with
#' arrows and text?
#' @param show.env \code{Logical}. Should the environments be highlighted with
#' arrows and text?
#'
#' @import dplyr
#' @import ggplot2
#'
#' @export
#'
plot_AMMI <- function(AMMI.out, x = 0, y = 1, show.gen = TRUE, show.env = TRUE) {

  # Error handling
  ## AMMI.out must be of class "AMMI"
  if (class(AMMI.out != "AMMI"))
    stop("'AMMI.out' must be of class 'AMMI.'")

  # Extract the useful information from the output
  analysis <- AMMI.out$analysis
  biplot <- AMMI.out$biplot

  # How many PCs are in the output?
  maxPC <- nrow(analysis)

  # Neither x nor y can be greater than maxPC
  if (x > maxPC)
    stop(str_c("'x' cannot be greater than the number of PCs: ", maxPC, "."))
  if (y > maxPC)
    stop(str_c("'y' cannot be greater than the number of PCs: ", maxPC, "."))
  # x cannot be less than 0
  if (x < 0)
    stop("'x' cannot be less than 0.")
  # y cannot be less than 1
  if (x < 1)
    stop("'y' cannot be less than 1.")

  # Tidy the biplot data
  biplot1 <- data.frame(
    name = row.names(biplot),
    biplot
  ) %>%
    tbl_df()

  # Change the third column to phenotypic mean
  names(biplot1)[3] <- "phenotypic_mean"

  # Assign x and y axis names
  x.name <- names(biplot1)[3 + x]
  y.name <- names(biplot1)[3 + y]
  # Designator for genotype/environment
  desig <- "type"
  name <- "name"

  # Pull out the proportion of variance explained by the PCs, if applicable
  if (x.name == "phenotypic_mean") {
    x.prop <- ""
    x.lab <- str_c("Phenotypic Mean ", x.prop)
  } else {
    x.prop <- str_c("(", analysis[x.name,"percent"], "%)")
    x.lab <- x.name
  }

  y.prop <- str_c("(", analysis[y.name,"percent"], "%)")
  y.lab <- str_c(y.name, y.prop, sep = " ")

  # Calculate x and y intercepts for straight lines
  x.inter <- biplot1 %>% select_(x.name) %>% unlist() %>% mean()
  y.inter <- biplot1 %>% select_(y.name) %>% unlist() %>% mean()


  # Start plotting
  gp <- biplot1 %>%
    ggplot(aes_string(x = x.name, y = y.name)) +
    geom_vline(xintercept = x.inter) +
    geom_hline(yintercept = y.inter) +
    geom_point(aes(col = type, shape = type), size = 2) +
    xlab(x.lab) +
    ylab(y.lab) +
    scale_color_discrete(guide = guide_legend(title = "Type")) +
    scale_shape_discrete(guide = guide_legend(title = "Type"))

  # Add arrows/text if requested
  if (show.gen) {
    # Create vectors of x1, x2, y1, y2
    x1 <- x.inter
    y1 <- y.inter
    x2 <- ifelse(biplot1$type == "GEN", biplot1[[x.name]], x1)
    y2 <- ifelse(biplot1$type == "GEN", biplot1[[y.name]], y1)

    gp1 <- gp +
      # Add arrows
      geom_segment(aes(x = x1, xend = x2, y = y1, yend = y2),
                   arrow = arrow(length = unit(0.01, "npc"))) +
      # Add genotype names
      geom_text(data = subset(biplot1, type == "GEN"),
                aes_string(x.name, y.name, label = name),
                size = 3, col = "red", check_overlap = T, nudge_y = -0.1)

  } else {
    # Reassign
    gp1 <- gp
  }

  if (show.env) {
    # Create vectors of x1, x2, y1, y2
    x1 <- x.inter
    y1 <- y.inter
    x2 <- ifelse(biplot1$type == "ENV", biplot1[[x.name]], x1)
    y2 <- ifelse(biplot1$type == "ENV", biplot1[[y.name]], y1)

    gp2 <- gp1 +
      # Add arrows
      geom_segment(aes(x = x1, xend = x2, y = y1, yend = y2),
                   arrow = arrow(length = unit(0.01, "npc"))) +
      # Add genotype names
      geom_text(data = subset(biplot1, type == "ENV"),
                aes_string(x.name, y.name, label = name),
                size = 3, col = "blue", check_overlap = T, nudge_y = -0.1)

  } else {
    # Reassign
    gp2 <- gp1

  }

  # Plot!
  print(gp2)

} # Close the function








