library(lubridate)
library(ggplot2)

source("pav_algorithms.R")
source("generate_synthetic_data.R")

compare_runtime <- function(pava_list, gen_data_list, n_fnc, plot_name, B = 1) {
  # - pava_list:
  # list of pava implementations following format used in pav_algorithms.R
  # - gen_data_list:
  # list of data generating functions following format used in generate_synthetic_data.R
  # - n_fnc:
  # function getting name of algorithm from pava_list and name of data generating
  # process from gen_data_list and returning sizes of data to test algorithm on
  # - plot_name:
  # name of the plot comparing runtimes
  # - my_colors:
  # named vector assigning algorithms in pava_list colors
  # - my_shapes:
  # named vector assigning algorithms in pava_list shape values
  # - B
  # number of times used to average runtimes
  collect_results <- data.frame()

  for (a in names(pava_list)) {
    cat("\n", a, ":\n", sep = "")
    for (dgp in names(gen_data_list)) {
      cat("-", dgp, " ", sep = "")

      n_vec <- n_fnc(a, dgp)

      for (n in n_vec) {
        data <- gen_data_list[[dgp]](n)
        time <- 0
        # average runtimes
        for (b in 1:B) {
          start <- now()
          res <- pava_list[[a]](data$x, data$y)
          end <- now()
          time <- time + as.numeric(seconds(end - start))
        }

        collect_results <- rbind(
          collect_results,
          data.frame(n = n, T = time / B, Data = dgp, Algorithm = a)
        )
      }
    }
  }
  collect_results$T[collect_results$T == 0] <- 10^-5   # for log transform
  ggplot(collect_results, aes(x = n, y = T, color = Algorithm, linetype = Data)) +
    geom_line(size = 0.5, aes(group = paste0(Algorithm, ":", Data))) +
    xlab("n") +
    scale_x_log10() +
    ylab("runtime [s]") +
    scale_y_log10() +
    ggtitle("Runtime of PAV Implementations on Synthetic Data") +
    theme_bw() +
    theme(legend.justification = c(1, 0), legend.position = c(0.99, 0.01),
          legend.box = "horizontal", legend.background = element_blank(),
          legend.box.background = element_rect(color = "gray"),
          panel.grid.major = element_line(size = 0.05),
          panel.grid.minor = element_line(size = 0.05),
          strip.background = element_blank())
  ggsave(plot_name, width = 250, height = 150, unit = "mm")
  invisible(collect_results)
}


# data generating processes
l_dgp <- list(IncreasingNormal = gen_inc_norm,
              UShapedNormal = gen_u_norm,
              InvUShapedNormal = gen_inv_u_norm,
              IncreasingNormalTies = gen_ties_norm)
l_algo <- list(GPAVA = wrap_gpava,
               MY_R = my_r_isoreg,
               ISOREG = wrap_isoreg,
               MY_CPP = my_cpp_isoreg)
num_avg <- 10

# n         10^7      10^8
# memory    160 MB    1.6 GB

# function that determines for each algorithm and data generating process
# sizes of data to measure runtime on
get_data_sizes <- function(name_algo, name_dgp) {
  n_min <- 2
  n_matrix <- matrix(c(4, 4, 4, 4,
                       5, 5, 5, 6,
                       6, 7, 6, 6,
                       7, 7, 7, 7), byrow = T, nrow = 4, ncol = 4)

  rownames(n_matrix) <- names(l_algo)
  colnames(n_matrix) <- names(l_dgp)
  return(10^(n_min:n_matrix[name_algo, name_dgp]))
}

set.seed(999)

compare_runtime(l_algo, l_dgp, get_data_sizes, "figures/joint_analysis.pdf",
                B = num_avg)
