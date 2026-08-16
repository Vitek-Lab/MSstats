#' Stable Givens rotation
#'
#' Given two numbers \code{a} and \code{b}, finds a rotation (a cosine and
#' a sine) that zeroes out \code{b}: it returns \code{cos_val}, \code{sin_val}
#' and \code{radius} such that
#' \code{cos_val * a + sin_val * b == radius} and
#' \code{-sin_val * a + cos_val * b == 0}. Used by \code{.lsqr} at every
#' iteration to eliminate one entry of the (otherwise upper-bidiagonal)
#' matrix built up during Golub-Kahan bidiagonalization. Computing the
#' rotation this way (dividing by whichever of \code{a}/\code{b} is larger
#' in magnitude) avoids the overflow/precision loss a naive
#' \code{a / sqrt(a^2 + b^2)} formula can suffer. Ported from
#' \code{_sym_ortho} in \code{scipy/sparse/linalg/_isolve/lsqr.py}.
#'
#' @param a,b numeric scalars
#' @return list with elements \code{cos_val}, \code{sin_val}, \code{radius}
#' @keywords internal
.symOrtho = function(a, b) {
    if (b == 0) {
        list(cos_val = sign(a), sin_val = 0, radius = abs(a))
    } else if (a == 0) {
        list(cos_val = 0, sin_val = sign(b), radius = abs(b))
    } else if (abs(b) > abs(a)) {
        ratio = a / b
        sin_val = sign(b) / sqrt(1 + ratio ^ 2)
        cos_val = sin_val * ratio
        list(cos_val = cos_val, sin_val = sin_val, radius = b / sin_val)
    } else {
        ratio = b / a
        cos_val = sign(a) / sqrt(1 + ratio ^ 2)
        sin_val = cos_val * ratio
        list(cos_val = cos_val, sin_val = sin_val, radius = a / cos_val)
    }
}


#' Solve a linear least-squares system with LSQR
#'
#' Finds the (least-squares) solution of \code{A x = b} using LSQR
#' (C. C. Paige and M. A. Saunders, "LSQR: An algorithm for sparse linear
#' equations and sparse least squares", ACM TOMS 8(1), 43-71, 1982), an
#' iterative method built on Golub-Kahan bidiagonalization. Ported from
#' \code{lsqr()} in \code{scipy/sparse/linalg/_isolve/lsqr.py}.
#'
#' The key property this function relies on: at no point does the
#' algorithm need \code{A} itself, or any factorization of it — every
#' step only ever multiplies a vector by \code{A} or by its transpose.
#' Those two operations are passed in as plain R functions
#' (\code{matvec}/\code{rmatvec}) rather than \code{A} being a matrix
#' argument, mirroring scipy's \code{LinearOperator} abstraction. This
#' is what lets a caller solve a system whose matrix is sparse and
#' structured (like the one-hot design matrices built by
#' \code{.buildDesignOperator}) without ever allocating the matrix in
#' memory: \code{matvec}/\code{rmatvec} can be written to exploit that
#' structure directly.
#'
#' Algorithm sketch: each iteration extends an orthogonal
#' bidiagonalization of \code{A} by one step (alternately calling
#' \code{matvec} and \code{rmatvec}), then applies one more Givens
#' rotation to fold the newly-revealed matrix entry into a running
#' upper-bidiagonal system, which immediately gives the next update to
#' the solution \code{x} and to the running estimates of the residual
#' norm, \code{norm(A)}, and \code{cond(A)} used for the convergence
#' tests below.
#'
#' @param matvec function(v) computing \code{A \%*\% v} for a length-n
#' vector \code{v}, returning a length-m vector
#' @param rmatvec function(u) computing \code{t(A) \%*\% u} for a length-m
#' vector \code{u}, returning a length-n vector
#' @param b numeric vector, length m (the right-hand side)
#' @param n number of unknowns (i.e. number of columns of the implicit A)
#' @param damp damping coefficient (Tikhonov regularization), default 0
#' @param atol,btol stopping tolerances on the backward error, default 1e-10
#' @param conlim stopping tolerance on the estimated condition number of A,
#' default 1e8
#' @param iter_lim maximum number of iterations, default \code{4 * n}
#'
#' @return list with elements \code{x} (solution vector), \code{istop}
#' (integer code for why iteration stopped: 1/2 = converged, 3/6 =
#' condition number too large, 4/5 = numerical stagnation, 7 = hit
#' \code{iter_lim}) and \code{itn} (number of iterations actually performed)
#' @keywords internal
.lsqr = function(matvec, rmatvec, b, n, damp = 0, atol = 1e-10, btol = 1e-10,
                 conlim = 1e8, iter_lim = NULL) {
    if (is.null(iter_lim)) {
        iter_lim = 4L * n
    }
    eps = .Machine$double.eps
    condition_tolerance = if (conlim > 0) 1 / conlim else 0
    damp_squared = damp ^ 2

    # Running estimates, updated incrementally every iteration rather than
    # recomputed from scratch (that's the whole point of the method: no
    # step ever looks back at A or at previous vectors).
    a_norm_estimate = 0            # estimate of norm(A) (or norm(A) padded with damp*I)
    a_cond_estimate = 0            # estimate of cond(A)
    direction_norm_sq = 0          # accumulates norm(A'A)^{-1}-related term used for a_cond_estimate
    residual_tail_sq = 0           # accumulates the part of the residual explained by damping
    x_norm = 0
    x_norm_sq = 0
    z = 0
    prev_cos = -1
    prev_sin = 0

    solution = rep(0, n)

    # Set up the first bidiagonalization vectors: beta*u = b (starting
    # from x = 0), alpha*v = A' u. These get extended by one more (beta,
    # u, alpha, v) each iteration below.
    u = b
    b_norm = sqrt(sum(b ^ 2))
    beta = b_norm

    if (beta > 0) {
        u = u / beta
        v = rmatvec(u)
        alpha = sqrt(sum(v ^ 2))
    } else {
        v = rep(0, n)
        alpha = 0
    }
    if (alpha > 0) {
        v = v / alpha
    }
    search_direction = v

    rho_bar = alpha       # diagonal entry (not yet finalized by a rotation) of the bidiagonal system
    phi_bar = beta        # remaining unexplained part of norm(b), before this iteration's rotation
    stop_code = 0
    iteration = 0

    a_t_r_norm = alpha * beta   # estimate of norm(A' * residual); 0 means b = 0, x = 0 solves it exactly
    if (a_t_r_norm == 0) {
        return(list(x = solution, istop = 0, itn = 0))
    }

    while (iteration < iter_lim) {
        iteration = iteration + 1

        # --- Continue the Golub-Kahan bidiagonalization by one step ---
        # These recurrences produce the next orthonormal pair (u, v) and
        # scalars (beta, alpha) such that, if A were fully expanded,
        # A * v_k = alpha_k * v_k... no wait, this is the two-term
        # recurrence beta*u_{k+1} = A*v_k - alpha_k*u_k and
        # alpha_{k+1}*v_{k+1} = A'*u_{k+1} - beta_{k+1}*v_k.
        # Note `v` here is the bidiagonalization direction, kept distinct
        # from `search_direction` (the vector actually used to update the
        # solution below) even though they start out equal — they diverge
        # after the first iteration.
        u = matvec(v) - alpha * u
        beta = sqrt(sum(u ^ 2))
        if (beta > 0) {
            u = u / beta
            a_norm_estimate = sqrt(a_norm_estimate ^ 2 + alpha ^ 2 + beta ^ 2 + damp_squared)
            v = rmatvec(u) - beta * v
            alpha = sqrt(sum(v ^ 2))
            if (alpha > 0) {
                v = v / alpha
            }
        }

        # --- Fold in the damping term (Tikhonov regularization), if any ---
        # A rotation that absorbs damp*I into the current diagonal entry
        # before the main elimination step below.
        if (damp > 0) {
            damped_rho_bar = sqrt(rho_bar ^ 2 + damp_squared)
            damp_cos = rho_bar / damped_rho_bar
            damp_sin = damp / damped_rho_bar
            residual_tail = damp_sin * phi_bar
            phi_bar = damp_cos * phi_bar
        } else {
            damped_rho_bar = rho_bar
            residual_tail = 0
        }

        # --- Eliminate the newly-revealed sub-diagonal entry (beta) ---
        # This is the step that turns the lower-bidiagonal system built
        # up so far into an upper-bidiagonal one, one row at a time.
        rotation = .symOrtho(damped_rho_bar, beta)
        rotation_cos = rotation$cos_val
        rotation_sin = rotation$sin_val
        rho = rotation$radius

        super_diagonal = rotation_sin * alpha
        rho_bar = -rotation_cos * alpha
        phi = rotation_cos * phi_bar
        phi_bar = rotation_sin * phi_bar
        phi_alpha_product = rotation_sin * phi

        # --- Update the solution and the search direction ---
        step_length = phi / rho
        direction_coef = -super_diagonal / rho
        scaled_direction = search_direction / rho

        solution = solution + step_length * search_direction
        search_direction = v + direction_coef * search_direction
        direction_norm_sq = direction_norm_sq + sum(scaled_direction ^ 2)

        # --- Update running estimate of norm(x) via one more rotation ---
        rotated_rho = prev_sin * rho
        neg_rotated_gamma_bar = -prev_cos * rho
        partial_z = phi - rotated_rho * z
        z_bar = partial_z / neg_rotated_gamma_bar
        x_norm = sqrt(x_norm_sq + z_bar ^ 2)
        gamma = sqrt(neg_rotated_gamma_bar ^ 2 + super_diagonal ^ 2)
        prev_cos = neg_rotated_gamma_bar / gamma
        prev_sin = super_diagonal / gamma
        z = partial_z / gamma
        x_norm_sq = x_norm_sq + z ^ 2

        # --- Update norm/condition-number estimates and test convergence ---
        a_cond_estimate = a_norm_estimate * sqrt(direction_norm_sq)
        residual_head_sq = phi_bar ^ 2
        residual_tail_sq = residual_tail_sq + residual_tail ^ 2
        residual_norm = sqrt(residual_head_sq + residual_tail_sq)
        a_t_r_norm = alpha * abs(phi_alpha_product)

        compatible_test = residual_norm / b_norm
        least_squares_test = a_t_r_norm / (a_norm_estimate * residual_norm + eps)
        condition_test = 1 / (a_cond_estimate + eps)
        machine_precision_test = compatible_test / (1 + a_norm_estimate * x_norm / b_norm)
        residual_tolerance = btol + atol * a_norm_estimate * x_norm / b_norm

        if (iteration >= iter_lim) stop_code = 7
        if (1 + condition_test <= 1) stop_code = 6
        if (1 + least_squares_test <= 1) stop_code = 5
        if (1 + machine_precision_test <= 1) stop_code = 4
        if (condition_test <= condition_tolerance) stop_code = 3
        if (least_squares_test <= atol) stop_code = 2
        if (compatible_test <= residual_tolerance) stop_code = 1

        if (stop_code != 0) break
    }

    list(x = solution, istop = stop_code, itn = iteration)
}


#' Build a matrix-free operator for a "1 + run + feature" design
#'
#' The model \code{log2inty ~ run + feature} (treatment-contrast coded,
#' i.e. the first level of each factor is a "reference" level that gets no
#' dummy column of its own) has a design matrix where every row has at
#' most 3 nonzero entries: the intercept (always 1), at most one `run`
#' dummy, and at most one `feature` dummy. That means the matrix-vector
#' products \code{.lsqr} needs (\code{A \%*\% v} and \code{t(A) \%*\% u}) can
#' be computed directly from two integer group-membership vectors,
#' without ever building the matrix \code{A} itself:
#' \itemize{
#'   \item \code{A \%*\% v} just looks up, for every row, which run/feature
#'     coefficient (if any) applies and adds it to the intercept.
#'   \item \code{t(A) \%*\% u} sums \code{u} within each run group and
#'     within each feature group (plus an overall sum for the intercept).
#' }
#' \code{run_levels}/\code{feature_levels} let the caller fix the full set
#' of levels (and which one is the reference) independently of which rows
#' are actually being fit — matching how \code{stats::model.matrix} derives
#' factor levels from a column before any rows are dropped for missing
#' values, so that a level with zero surviving observations still gets an
#' (all-zero, coefficient pinned near 0 by \code{.lsqr}'s minimum-norm
#' behavior) column rather than silently changing the model's
#' parameterization.
#'
#' @param run,feature vectors (same length), the two grouping variables
#' @param run_levels,feature_levels optional explicit level sets (else
#' taken from \code{sort(unique(run))} / \code{sort(unique(feature))} via
#' \code{factor()}'s default ordering, matching \code{model.matrix}'s
#' default treatment-contrast reference level)
#'
#' @return list with elements \code{matvec}, \code{rmatvec} (functions,
#' see \code{.lsqr}) and \code{p} (number of coefficients)
#' @keywords internal
.buildDesignOperator = function(run, feature, run_levels = NULL, feature_levels = NULL) {
    run_factor = if (is.null(run_levels)) factor(run) else factor(run, levels = run_levels)
    feature_factor = if (is.null(feature_levels)) factor(feature) else factor(feature, levels = feature_levels)

    n = length(run)
    n_run = nlevels(run_factor)
    n_feature = nlevels(feature_factor)

    # 0 = reference level (dropped by treatment contrasts, contributes
    # nothing beyond the intercept); 1..(n_run-1) indexes the non-reference
    # run coefficients, and likewise for feature.
    run_idx = as.integer(run_factor) - 1L
    feature_idx = as.integer(feature_factor) - 1L

    run_offset = 1L
    feature_offset = 1L + (n_run - 1L)
    p = feature_offset + (n_feature - 1L)

    matvec = function(coef) {
        fitted = rep(coef[1], n)
        has_run = run_idx > 0L
        fitted[has_run] = fitted[has_run] + coef[run_offset + run_idx[has_run]]
        has_feature = feature_idx > 0L
        fitted[has_feature] = fitted[has_feature] + coef[feature_offset + feature_idx[has_feature]]
        fitted
    }

    rmatvec = function(resid) {
        grad = numeric(p)
        grad[1] = sum(resid)
        if (n_run > 1L) {
            run_sums = rowsum(resid, run_idx)
            group_labels = as.integer(rownames(run_sums))
            nonref = group_labels > 0L
            grad[run_offset + group_labels[nonref]] = run_sums[nonref, 1]
        }
        if (n_feature > 1L) {
            feature_sums = rowsum(resid, feature_idx)
            group_labels = as.integer(rownames(feature_sums))
            nonref = group_labels > 0L
            grad[feature_offset + group_labels[nonref]] = feature_sums[nonref, 1]
        }
        grad
    }

    list(matvec = matvec, rmatvec = rmatvec, p = p)
}


#' Apply IRLS weights to a design operator
#'
#' Represents \code{diag(sqrt(w)) \%*\% A} without ever forming that
#' product: \code{matvec} scales the operator's output by \code{sqrt(w)},
#' \code{rmatvec} scales its input by \code{sqrt(w)} first. Rebuilding
#' just this thin wrapper each IRLS iteration (rather than
#' \code{.buildDesignOperator} itself) means the \code{run}/\code{feature}
#' group-index bookkeeping is done exactly once per protein, not once per
#' iteration.
#'
#' @param op a design operator (list with \code{matvec}/\code{rmatvec}/\code{p},
#' see \code{.buildDesignOperator})
#' @param w numeric vector of IRLS weights, one per observation
#' @return operator (list with \code{matvec}/\code{rmatvec}/\code{p}) for
#' the reweighted system
#' @keywords internal
.weightOperator = function(op, w) {
    sqrt_w = sqrt(w)
    list(
        matvec = function(coef) sqrt_w * op$matvec(coef),
        rmatvec = function(resid) op$rmatvec(sqrt_w * resid),
        p = op$p
    )
}


#' Fit a Huber M-estimator via IRLS, solving each weighted least-squares
#' step with LSQR against a matrix-free design operator
#'
#' Reimplements the "M" estimation / "Huber" scale-estimate branch of
#' \code{MASS::rlm.default} (see \code{MASS/R/rlm.R}) for the specific
#' model \code{log2inty ~ run + feature}, replacing every call to
#' \code{stats::lm.wfit(x, y, w, method = "qr")} with a weighted
#' \code{.lsqr()} solve against \code{.buildDesignOperator}'s matrix-free
#' operator. Only the code path used by \code{.fitHuber} is reproduced:
#' \code{init = "ls"}, \code{psi = MASS::psi.huber}, \code{scale.est =
#' "Huber"}, no case/prior weights.
#'
#' @param log2inty numeric response vector
#' @param run,feature grouping vectors, same length as \code{log2inty}
#' @param run_levels,feature_levels optional explicit level sets, see
#' \code{.buildDesignOperator}
#' @param k2 tuning constant for Huber's Proposal 2 scale estimate
#' @param maxit maximum number of IRLS iterations
#' @param acc IRLS convergence tolerance on the relative change in residuals
#'
#' @return list with elements \code{coefficients}, \code{residuals},
#' \code{fitted.values}, \code{scale}, \code{df.residual}, \code{rank} and
#' \code{converged}
#' @importFrom MASS psi.huber
#' @importFrom stats mad pnorm dnorm
#' @keywords internal
.rlmLSQR = function(log2inty, run, feature, run_levels = NULL, feature_levels = NULL,
                    k2 = 1.345, maxit = 20, acc = 1e-4) {
    y = log2inty
    n = length(y)
    op = .buildDesignOperator(run, feature, run_levels, feature_levels)
    p = op$p
    n1 = n - p

    coef = .lsqr(op$matvec, op$rmatvec, y, p)$x
    fitted = op$matvec(coef)
    resid = y - fitted
    scale = stats::mad(resid, center = 0)

    theta = 2 * stats::pnorm(k2) - 1
    gamma_const = theta + k2 ^ 2 * (1 - theta) - 2 * k2 * stats::dnorm(k2)

    done = FALSE
    for (iiter in seq_len(maxit)) {
        prev_resid = resid

        scale = sqrt(sum(pmin(resid ^ 2, (k2 * scale) ^ 2)) / (n1 * gamma_const))
        if (scale == 0) {
            done = TRUE
            break
        }

        w = MASS::psi.huber(resid / scale)
        weighted_op = .weightOperator(op, w)
        coef = .lsqr(weighted_op$matvec, weighted_op$rmatvec, sqrt(w) * y, p)$x
        fitted = op$matvec(coef)
        resid = y - fitted

        convi = sqrt(sum((prev_resid - resid) ^ 2) / max(1e-20, sum(prev_resid ^ 2)))
        done = convi <= acc
        if (done) break
    }

    list(
        coefficients = coef,
        residuals = resid,
        fitted.values = fitted,
        scale = scale,
        df.residual = n1,
        rank = p,
        converged = done
    )
}
