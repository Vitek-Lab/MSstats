#' @keywords internal
anomaly_weights_z_vec <- function(s, p = 2, eps = 1e-8,
                                  s_floor = 0.1,
                                  w_min = 0.05, w_max = 5,
                                  normalize_sum_to_n = TRUE) {
  ok <- !is.na(s)
  w  <- rep(NA_real_, length(s))
  
  if (sum(ok) < 3L) {
    w[ok] <- 1
    return(w)
  }
  
  med <- stats::median(s[ok])
  mad <- stats::mad(s[ok], center = med, constant = 1.4826)
  
  scale <- max(mad, s_floor)
  
  z <- (s - med) / (scale + eps)
  z_pos <- pmax(0, z)
  
  w_star <- 1 / (1 + (z_pos^p))
  
  w_star <- pmax(w_star, w_min)
  w_star <- pmin(w_star, w_max)
  
  if (normalize_sum_to_n) {
    w_star[ok] <- w_star[ok] * sum(ok) / sum(w_star[ok])
  }
  
  w[ok] <- w_star[ok]
  w
}