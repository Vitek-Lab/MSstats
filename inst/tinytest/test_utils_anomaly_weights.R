# Test anomaly_weights_z_vec ------------------------------------------------

# Test 1: fewer than 3 non-NA values returns all weights of 1 for non-NA entries
s_few = c(1, NA, 3)
w_few = MSstats:::anomaly_weights_z_vec(s_few)
expect_equal(w_few, c(1, NA, 1))

# Test 2: all NA
s_all_na = c(NA_real_, NA_real_)
w_all_na = MSstats:::anomaly_weights_z_vec(s_all_na)
expect_true(all(is.na(w_all_na)))

# Test 3: no anomalies - roughly homogeneous values should get similar weights near 1
s_homog = c(10, 10.1, 9.9, 10.05, 9.95)
w_homog = MSstats:::anomaly_weights_z_vec(s_homog)
expect_equal(length(w_homog), length(s_homog))
expect_true(all(!is.na(w_homog)))
expect_true(all(w_homog >= 0.05 & w_homog <= 5))

# Test 4: a clear anomaly should receive a lower weight than the rest
s_anomaly = c(10, 10.1, 9.9, 10.05, 100)
w_anomaly = MSstats:::anomaly_weights_z_vec(s_anomaly)
expect_true(w_anomaly[5] < min(w_anomaly[1:4]))

# Test 5: normalize_sum_to_n = FALSE skips renormalization step (weights unbounded by n)
w_no_norm = MSstats:::anomaly_weights_z_vec(s_anomaly, normalize_sum_to_n = FALSE)
w_norm = MSstats:::anomaly_weights_z_vec(s_anomaly, normalize_sum_to_n = TRUE)
expect_true(!isTRUE(all.equal(w_no_norm, w_norm)))

# Test 6: NA values are preserved as NA in output, non-NA positions filled
s_mixed = c(1, 2, NA, 4, 5, 100)
w_mixed = MSstats:::anomaly_weights_z_vec(s_mixed)
expect_true(is.na(w_mixed[3]))
expect_true(all(!is.na(w_mixed[-3])))
