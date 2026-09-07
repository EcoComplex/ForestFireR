#' Render an animation of a spatial simulation
#'
#' Runs the same Gillespie dynamics as [simulate_spatial()], but instead of
#' returning only a final state (or the throttled fire-activity snapshots
#' from `capture_fire_snapshots`), captures the full `L`x`L` grid at
#' `n_frames` evenly spaced points in SIMULATED TIME across `[0, T]` and
#' encodes them into a video (default) or an animated GIF, returning the
#' path to the rendered file.
#'
#' Internally, this chains calls to [simulate_spatial_from_grid()], each one
#' asking to advance exactly to the next unfilled checkpoint time
#' `k * T / n_frames`, starting from the previous call's `final_grid` -- the
#' Gillespie process is Markovian, so this produces a trajectory
#' statistically identical to one continuous run.
#'
#' A single Gillespie step cannot be split: [simulate_spatial_from_grid()]'s
#' engine always applies whatever reaction it draws next in full, even when
#' that reaction's waiting time overshoots the requested horizon (there is
#' no way to ask it "how much time until the next event, without applying
#' it" -- that would need a change to the underlying C++, which this
#' function deliberately avoids). So a single call can overshoot its
#' requested checkpoint by a lot, especially early on or with slow rate
#' constants (no reaction is possible faster than the model's own smallest
#' active rate allows) -- e.g. requesting `T = 10` sliced into 30 frames
#' asks for one reaction roughly every `0.33` time units, but if the
#' landscape's next actual reaction is `12` time units away (very plausible
#' with, say, `default_rates()`'s regrowth/replacement rates, which the
#' paper itself only shows producing visible change over `T` on the order
#' of `100`-`300`), that one call alone jumps past `12/0.33 ~ 36` requested
#' checkpoints at once. This function handles that correctly: whenever a
#' call's actual elapsed time (`time_sim`) reaches or passes one or more
#' pending checkpoints, ALL of them are filled with that call's resulting
#' grid (correctly showing no change happened yet), and the next call
#' starts from the next still-unfilled checkpoint -- so the total number of
#' engine calls actually made is bounded by the number of *real* reactions
#' needed to cover `[0, T]`, which can be far fewer than `n_frames` for slow
#' dynamics (efficient) or up to `n_frames` for fast ones (e.g. an active
#' fire), never more. The practical takeaway: if you see very few distinct
#' frames (or a `verbose` log with very few calls before reaching `T`),
#' that's the model genuinely not doing much in `[0, T]` at these rates, not
#' a bug -- pick a larger `T` (see [default_rates()]'s documentation and
#' the paper's own figures for typical horizons) rather than a smaller
#' `n_frames`.
#'
#' Each engine call gets its own derived seed (`seed + <call index>`, when
#' `seed` is supplied) rather than reusing one seed repeatedly, which would
#' otherwise restart the underlying C RNG from the same state on every
#' call.
#'
#' If the native population goes extinct (or the process reaches an
#' absorbing state, `sumA == 0`) before `T`, matching [simulate_spatial()]'s
#' own `check_extinction` early stop, all remaining checkpoints are filled
#' with that final grid, so the rendered video still runs the full
#' requested length, ending on a frozen final state, instead of silently
#' being shorter than asked.
#'
#' @inheritParams simulate_spatial
#' @param initial_grid An `L`x`L` integer matrix using [FF_STATE]'s codes
#'   (typically from [generate_landscape_layers()]), used instead of
#'   building one internally from `density2`/`p` -- e.g. to seed an active
#'   fire directly instead of waiting on spontaneous ignition (`Lig_23`).
#'   When supplied, `L`, `density2` and `p` are ignored. Remember to also
#'   set `check_extinction = FALSE` if this landscape has no native
#'   vegetation by design (see [simulate_spatial_from_grid()]'s own docs
#'   for this same recipe).
#' @param n_frames Number of checkpoints to capture across `[0, T]` (the
#'   grid is also captured at `t = 0`, so the rendered animation has
#'   `n_frames + 1` frames). See Details for what happens when the model's
#'   actual dynamics are slower or faster than one reaction per
#'   `T / n_frames`.
#' @param fps Frames per second in the rendered output.
#' @param file Output path. The extension picks the renderer: `.gif` uses
#'   the `gifski` package, anything else (`.mp4`, `.mov`, `.webm`, ...) uses
#'   the `av` package. Defaults to a temporary `.mp4` file if `NULL`.
#' @param colors Named character vector of fill colors, keyed by
#'   [FF_STATE] names (`native`, `invader`, `fire`, `empty`,
#'   `empty_postfire`). Only the names you supply override the package
#'   default; you don't need to specify all five.
#' @param check_extinction Stop early once native vegetation goes extinct,
#'   matching [simulate_spatial()]'s own always-on behavior (see Details
#'   for what happens to the rendered video when this triggers). Set
#'   `FALSE` to keep simulating -- and rendering -- past extinction.
#' @param width,height Frame size in pixels.
#' @param verbose If `TRUE` (default), print progress as checkpoints are
#'   filled.
#' @return Invisibly, the path to the rendered video/GIF file (`file`,
#'   resolved to an absolute path), with attributes `n_engine_calls` (how
#'   many [simulate_spatial_from_grid()] calls it actually took to cover
#'   `[0, T]` -- compare to `n_frames` per Details), `stopped_early`
#'   (`TRUE` if `check_extinction` triggered before `T`) and `t_elapsed`
#'   (the actual simulated time reached, `== T` unless `stopped_early`).
#'
#' @examples
#' \dontrun{
#' # Needs the 'av' package: install.packages("av")
#' # T chosen large enough for these default rates to actually do
#' # something visible -- see Details.
#' path <- simulate_spatial_movie(T = 150, L = 80, density2 = 0.05, p = 1,
#'                                 n_frames = 150, fps = 15,
#'                                 file = "invasion.mp4", seed = 1)
#'
#' # Needs the 'gifski' package: install.packages("gifski")
#' simulate_spatial_movie(T = 150, L = 60, n_frames = 80,
#'                         file = "invasion.gif", seed = 1)
#'
#' # Fire spread, seeded directly (default L_30/Lig_23 make fire resolve
#' # in ~1e-4 years -- invisible to any video; see fig10/fig11's own
#' # reduced-dynamics recipe, and this function's Details). L_30 = 5500
#' # here is grounded in real fire-residence-time estimates for a 30x30 m
#' # cell (this package's own cell size), not a fitted package value --
#' # two independent estimates converge on it: (1) understory fire
#' # rate-of-spread in tropical forest, ~1-3 m/min, takes ~10-30 min to
#' # cross 30 m; (2) Saravia et al. 2025 (Oikos, doi:10.1111/oik.10764)
#' # fit a contact-spread fire model (same mechanism as this package's
#' # fire rule) on a 460 m/pixel Amazon grid with a fixed 1-day-per-site
#' # fire duration; rescaling that linearly by cell-size ratio
#' # (460/30 ~ 15x) gives ~1.6 h per 30 m cell. Both give a per-cell fire
#' # duration on the order of 10 min-2 h, i.e. L_30 ~ 5e3-5e4/year -- far
#' # from the package's L_30 = 1e6 default (effectively instantaneous)
#' # but also far from the earlier ad hoc guess of L_30 = 75 (~1 week per
#' # cell, too slow by ~2 orders of magnitude). T is rescaled down from
#' # years to ~15 days accordingly, so the fire front is still visibly
#' # spreading across the grid rather than resolving in a single frame or
#' # taking years. xi_inv still sets the spread-vs-burnout probability
#' # exactly as documented in xi2lambda() regardless of this rescaling.
#' ignite <- generate_landscape_layers(
#'   L = 80, fill_state = FF_STATE["invader"],
#'   layers = list(list(background = FF_STATE["invader"],
#'                       pattern = FF_STATE["fire"], density = 0.01, p = 1)),
#'   seed = 1
#' )
#' simulate_spatial_movie(T = 0.04, initial_grid = ignite,
#'                         xi_nat = 0, xi_inv = 0.8, eta_inv = 0,
#'                         L_01 = 0, L_02 = 0, L_12 = 0, L_21 = 0,
#'                         Lig_13 = 0, Lig_23 = 0, L_30 = 5500,
#'                         check_extinction = FALSE,
#'                         n_frames = 150, fps = 20,
#'                         file = "fire_spread.mp4", seed = 1)
#' }
#' @export
simulate_spatial_movie <- function(T, L = 100, density2 = 0.1, p = 1,
                                    xi_nat = 0.5, xi_inv = 0.6, eta_inv = 0.6,
                                    L_01 = 0.03, L_02 = 0.03, L_12 = 0.005, L_21 = 0.01,
                                    L_30 = 1e6, Lig_13 = 0, Lig_23 = 1e-4,
                                    periodic = TRUE, seed = NULL,
                                    initial_grid = NULL,
                                    n_frames = 150, fps = 15,
                                    file = NULL, colors = NULL,
                                    check_extinction = TRUE,
                                    width = 600, height = 600, verbose = TRUE) {
  if (n_frames < 1) stop("simulate_spatial_movie(): n_frames must be >= 1.")
  if (T <= 0) stop("simulate_spatial_movie(): T must be > 0.")

  file <- .ffr_movie_resolve_file(file)
  renderer <- .ffr_movie_pick_renderer(file)

  state_colors <- .ffr_movie_resolve_colors(colors)

  targets <- (seq_len(n_frames) / n_frames) * T  # checkpoint times, last one == T
  base_seed <- if (is.null(seed)) NULL else as.integer(seed)

  if (is.null(initial_grid)) {
    grid <- generate_landscape(L = L, density2 = density2, p = p, seed = base_seed)
  } else {
    if (nrow(initial_grid) != ncol(initial_grid)) {
      stop("simulate_spatial_movie(): initial_grid must be square.")
    }
    grid <- initial_grid
  }
  storage.mode(grid) <- "integer"

  frames <- vector("list", n_frames + 1L)
  frames[[1]] <- grid

  cur_time <- 0
  frame_idx <- 1L    # next unfilled checkpoint (targets[frame_idx] -> frames[[frame_idx + 1]])
  stopped_early <- FALSE
  call_idx <- 0L

  while (frame_idx <= n_frames) {
    chunk_T <- targets[frame_idx] - cur_time
    call_idx <- call_idx + 1L
    chunk_seed <- if (is.null(base_seed)) NULL else base_seed + call_idx

    r <- simulate_spatial_from_grid(
      T = chunk_T, initial_grid = grid,
      xi_nat = xi_nat, xi_inv = xi_inv, eta_inv = eta_inv,
      L_01 = L_01, L_02 = L_02, L_12 = L_12, L_21 = L_21,
      L_30 = L_30, Lig_13 = Lig_13, Lig_23 = Lig_23,
      periodic = periodic, seed = chunk_seed,
      record_dt = -1, record_grid = TRUE,
      check_extinction = check_extinction
    )

    new_time <- cur_time + r$time_sim
    grid <- r$final_grid
    storage.mode(grid) <- "integer"

    if (r$time_sim < chunk_T * (1 - 1e-9)) {
      # Absorbing state (or check_extinction) reached before the next
      # checkpoint -- freeze every remaining checkpoint at this grid.
      stopped_early <- TRUE
      cur_time <- new_time
      for (j in frame_idx:n_frames) frames[[j + 1L]] <- grid
      if (verbose) {
        message(sprintf(
          "  stopped early at t=%.3f/%.3f (extinction or absorbing state) after %d engine call(s) -- remaining %d checkpoint(s) frozen at the final state.",
          new_time, T, call_idx, n_frames - frame_idx + 1L))
      }
      frame_idx <- n_frames + 1L
      break
    }

    # This one call may have jumped past several pending checkpoints at
    # once (see Details) -- fill all of them with the resulting grid.
    filled_from <- frame_idx
    while (frame_idx <= n_frames && targets[frame_idx] <= new_time + 1e-9) {
      frames[[frame_idx + 1L]] <- grid
      frame_idx <- frame_idx + 1L
    }
    cur_time <- new_time

    if (verbose) {
      n_filled <- frame_idx - filled_from
      message(sprintf("  engine call %d: t=%.3f/%.3f -- filled checkpoint(s) %d-%d/%d%s",
                       call_idx, cur_time, T, filled_from, frame_idx - 1L, n_frames,
                       if (n_filled > 1) sprintf(" (one reaction spanned %d checkpoints)", n_filled) else ""))
    }
  }

  .ffr_movie_render_frames(frames, file = file, fps = fps, colors = state_colors,
                            width = width, height = height, renderer = renderer,
                            verbose = verbose)

  file <- normalizePath(file, mustWork = TRUE)
  attr(file, "n_engine_calls") <- call_idx
  attr(file, "stopped_early") <- stopped_early
  attr(file, "t_elapsed") <- cur_time
  invisible(file)
}

.ffr_movie_resolve_colors <- function(colors) {
  default_colors <- c(native = "darkgreen", invader = "mediumseagreen",
                       fire = "red", empty = "gray85", empty_postfire = "gray60")
  if (!is.null(colors)) {
    unknown <- setdiff(names(colors), names(default_colors))
    if (length(unknown) > 0) {
      stop("simulate_spatial_movie(): unknown color name(s): ", paste(unknown, collapse = ", "),
           " -- must be one of: ", paste(names(default_colors), collapse = ", "))
    }
    default_colors[names(colors)] <- colors
  }
  stats::setNames(unname(default_colors), as.character(FF_STATE[names(default_colors)]))
}

.ffr_movie_resolve_file <- function(file) {
  if (is.null(file)) file <- tempfile("ffr_movie_", fileext = ".mp4")
  file <- path.expand(file)
  dir.create(dirname(file), showWarnings = FALSE, recursive = TRUE)
  file
}

.ffr_movie_pick_renderer <- function(file) {
  ext <- tolower(tools::file_ext(file))
  if (identical(ext, "gif")) {
    if (!requireNamespace("gifski", quietly = TRUE)) {
      stop("simulate_spatial_movie(): rendering a .gif needs the 'gifski' package -- install.packages(\"gifski\").")
    }
    return("gifski")
  }
  if (!requireNamespace("av", quietly = TRUE)) {
    stop("simulate_spatial_movie(): rendering a video (.", ext, ") needs the 'av' package -- install.packages(\"av\").",
         " Use a '.gif' file (needs the 'gifski' package instead) if you'd rather avoid it.")
  }
  "av"
}

# Draws one grid to a PNG using explicit integer-state breaks (never
# image()'s default continuous-range binning, which silently reassigns
# colors when not every state 1:5 is present in a given grid -- the same
# bug fixed in this project's fig9/fig11 spatial-snapshot panels).
.ffr_movie_write_frame_png <- function(grid, path, colors, width, height) {
  states <- sort(as.integer(names(colors)))
  cols <- colors[as.character(states)]
  breaks <- c(states - 0.5, states[length(states)] + 0.5)

  grDevices::png(path, width = width, height = height)
  on.exit(grDevices::dev.off(), add = TRUE)
  op <- graphics::par(mar = c(0, 0, 0, 0))
  on.exit(graphics::par(op), add = TRUE)
  graphics::image(x = seq_len(nrow(grid)), y = seq_len(ncol(grid)), z = grid,
                   col = cols, breaks = breaks,
                   axes = FALSE, xlab = "", ylab = "", asp = 1, useRaster = TRUE)
}

.ffr_movie_render_frames <- function(frames, file, fps, colors, width, height, renderer, verbose) {
  tmpdir <- tempfile("ffr_movie_frames_")
  dir.create(tmpdir)
  on.exit(unlink(tmpdir, recursive = TRUE), add = TRUE)

  n <- length(frames)
  digits <- nchar(as.character(n))
  png_files <- file.path(tmpdir, sprintf(paste0("frame_%0", digits, "d.png"), seq_len(n)))

  if (verbose) message(sprintf("  rendering %d frame(s) to PNG...", n))
  for (i in seq_len(n)) {
    .ffr_movie_write_frame_png(frames[[i]], png_files[i], colors = colors,
                                width = width, height = height)
  }

  if (verbose) message(sprintf("  encoding via '%s' -> %s", renderer, file))
  if (identical(renderer, "gifski")) {
    gifski::gifski(png_files, gif_file = file, width = width, height = height, delay = 1 / fps)
  } else {
    av::av_encode_video(png_files, output = file, framerate = fps)
  }
  invisible(file)
}
