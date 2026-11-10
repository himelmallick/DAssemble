# =============================================================================
# Extra sticker motifs for differential-abundance / ensemble packages.
# Each returns a transparent ggplot to pass as `subplot =` to make_pkg_sticker().
# =============================================================================

library(ggplot2)

METHOD_COLS <- c("#2EC4B6", "#4EA8DE", "#9B8CFF", "#F4B942", "#FF7A7A")

# -----------------------------------------------------------------------------
# Bagging-style ensemble: parallel "learners" (one mini volcano per DA
#    method) whose arrows converge on a combined output, either a lightbulb
#    (output = "bulb") or a consensus volcano with gold hits (output = "volcano")
# -----------------------------------------------------------------------------
.mini_volcano <- function(seed, n = 55, xmin, xmax, ymin, ymax, shift = 0) {
  set.seed(seed)
  lfc <- c(rnorm(n, 0, 0.8), rnorm(5, 2.5 + shift, 0.4), rnorm(5, -2.5 - shift, 0.4))
  se  <- runif(length(lfc), 0.5, 1.2)
  y   <- pmin(-log10(2 * pnorm(-abs(lfc) / se)), 5.5)
  lfc <- pmax(pmin(lfc, 3.9), -3.9)
  data.frame(x = xmin + (lfc + 4) / 8 * (xmax - xmin),
             y = ymin + y / 5.8 * (ymax - ymin),
             sig = y > 2 & abs(lfc) > 1.3)
}

motif_bagging <- function(output = c("bulb", "volcano"),
                          cols = METHOD_COLS[c(1, 3, 2)],
                          gold = "#F4B942") {
  output <- match.arg(output)
  yc <- c(5.0, 3.0, 1.0); h <- 1.55; x0 <- 0.2; x1 <- 3.3
  cards <- data.frame(xmin = x0, xmax = x1, ymin = yc - h / 2, ymax = yc + h / 2,
                      m = factor(1:3))
  pts <- do.call(rbind, lapply(1:3, function(i) {
    d <- .mini_volcano(10 + i, xmin = x0 + 0.2, xmax = x1 - 0.2,
                       ymin = yc[i] - h / 2 + 0.15, ymax = yc[i] + h / 2 - 0.12,
                       shift = 0.2 * (i - 2))
    transform(d, m = factor(i))
  }))
  node <- c(5.9, 3.0)
  arr <- data.frame(x = x1 + 0.15, y = yc, xend = node[1] - 0.05, yend = node[2])

  p <- ggplot() +
    geom_rect(data = cards, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
                                colour = m), fill = "white", alpha = 0.07,
              linewidth = 0.9) +
    geom_point(data = pts[!pts$sig, ], aes(x, y), colour = "white",
               alpha = 0.45, size = 0.35) +
    geom_point(data = pts[pts$sig, ], aes(x, y, colour = m), size = 0.9) +
    geom_curve(data = arr[1, ], aes(x, y, xend = xend, yend = yend),
               curvature = -0.25, colour = "white", linewidth = 0.7) +
    geom_curve(data = arr[3, ], aes(x, y, xend = xend, yend = yend),
               curvature = 0.25, colour = "white", linewidth = 0.7) +
    geom_segment(data = arr[2, ], aes(x, y, xend = xend, yend = yend),
                 colour = "white", linewidth = 0.7) +
    geom_point(aes(node[1], node[2]), colour = "white", size = 2.4) +
    scale_colour_manual(values = cols, guide = "none")

  if (output == "bulb") {
    t  <- seq(0, 2 * pi, length.out = 200)
    bc <- c(8.25, 3.45); r <- 1.35
    glass <- data.frame(x = bc[1] + r * cos(t), y = bc[2] + r * sin(t))
    glass <- glass[!(glass$y < bc[2] - r * 0.72 & abs(glass$x - bc[1]) < r * 0.55), ]
    neck  <- data.frame(x = bc[1] + c(-0.62, 0.62, 0.5, -0.5),
                        y = bc[2] - c(1.12, 1.12, 1.85, 1.85))
    rays  <- data.frame(a = seq(20, 160, length.out = 5) * pi / 180)
    rays  <- transform(rays, x = bc[1] + 1.7 * cos(a), y = bc[2] + 1.7 * sin(a),
                       xend = bc[1] + 2.15 * cos(a), yend = bc[2] + 2.15 * sin(a))
    p <- p +
      geom_segment(aes(node[1] + 0.1, node[2], xend = bc[1] - r - 0.15, yend = node[2]),
                   colour = "white", linewidth = 0.7,
                   arrow = grid::arrow(length = grid::unit(0.03, "npc"), type = "closed")) +
      geom_polygon(data = glass, aes(x, y), fill = gold, colour = NA) +
      geom_polygon(data = neck, aes(x, y), fill = "white", alpha = 0.85) +
      geom_segment(aes(x = bc[1] - 0.5, xend = bc[1] + 0.5,
                       y = bc[2] - c(1.37, 1.6), yend = bc[2] - c(1.37, 1.6)),
                   colour = "#0E2A47", linewidth = 0.5) +
      geom_segment(data = rays, aes(x, y, xend = xend, yend = yend),
                   colour = gold, linewidth = 0.8, lineend = "round") +
      # tiny volcano "filament" inside the bulb: the combined result
      geom_path(aes(x = bc[1] + seq(-0.55, 0.55, length.out = 50),
                    y = bc[2] - 0.35 + 0.75 * abs(seq(-0.55, 0.55, length.out = 50) / 0.55)^1.6),
                colour = "#0E2A47", linewidth = 0.7) +
      geom_point(aes(x = bc[1] + c(-0.62, -0.5, 0.5, 0.65), y = bc[2] + c(0.55, 0.3, 0.35, 0.6)),
                 colour = "#0E2A47", size = 0.9)
  } else {
    ox0 <- 6.6; ox1 <- 9.9; oy0 <- 1.2; oy1 <- 4.8
    big <- .mini_volcano(33, n = 110, xmin = ox0 + 0.2, xmax = ox1 - 0.2,
                         ymin = oy0 + 0.15, ymax = oy1 - 0.15)
    p <- p +
      geom_segment(aes(node[1] + 0.1, node[2], xend = ox0 - 0.12, yend = node[2]),
                   colour = "white", linewidth = 0.7,
                   arrow = grid::arrow(length = grid::unit(0.03, "npc"), type = "closed")) +
      geom_rect(aes(xmin = ox0, xmax = ox1, ymin = oy0, ymax = oy1),
                fill = "white", alpha = 0.10, colour = gold, linewidth = 1.1) +
      geom_point(data = big[!big$sig, ], aes(x, y), colour = "white",
                 alpha = 0.5, size = 0.45) +
      geom_point(data = big[big$sig, ], aes(x, y), colour = gold, size = 1.25)
  }
  p + coord_cartesian(xlim = c(0, 10.1), ylim = c(-0.05, 6.05), expand = FALSE) +
    theme_void() + hexSticker::theme_transparent()
}

# -----------------------------------------------------------------------------
# More bagging-style layouts (same story: method volcanoes -> combined volcano)
# -----------------------------------------------------------------------------
.card_layers <- function(xmin, xmax, ymin, ymax, seed, col, fill = "#0E2A47",
                         combined = FALSE, n = 55, lw = 0.9, shift = 0) {
  d <- .mini_volcano(seed, n = n, xmin = xmin + 0.18, xmax = xmax - 0.18,
                     ymin = ymin + 0.14, ymax = ymax - 0.12, shift = shift)
  box <- data.frame(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)
  list(
    geom_rect(data = box, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = fill, colour = NA),
    geom_rect(data = box, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = "white", alpha = if (combined) 0.10 else 0.07,
              colour = col, linewidth = if (combined) lw * 1.25 else lw),
    geom_point(data = d[!d$sig, ], aes(x, y), colour = "white",
               alpha = 0.45, size = if (combined) 0.45 else 0.35),
    geom_point(data = d[d$sig, ], aes(x, y), colour = col,
               size = if (combined) 1.25 else 0.9)
  )
}
.arrow <- function(x, y, xend, yend, curv = 0, head = TRUE) {
  a <- if (head) grid::arrow(length = grid::unit(0.03, "npc"), type = "closed")
  if (curv == 0) {
    geom_segment(aes(x = x, y = y, xend = xend, yend = yend), colour = "white",
                 linewidth = 0.7, arrow = a)
  } else {
    geom_curve(aes(x = x, y = y, xend = xend, yend = yend), colour = "white",
               linewidth = 0.7, curvature = curv, arrow = a)
  }
}
.finish <- function(p, xlim, ylim) {
  p + coord_cartesian(xlim = xlim, ylim = ylim, expand = FALSE) +
    theme_void() + hexSticker::theme_transparent()
}

# Funnel: three method volcanoes side by side on top, merging downwards
#    into one combined volcano
motif_funnel <- function(cols = METHOD_COLS[c(1, 3, 2)], gold = "#F4B942") {
  xs <- c(0.1, 3.55, 7.0); w <- 2.9; y0 <- 5.4; y1 <- 7.4
  node <- c(5, 4.3)
  p <- ggplot()
  for (i in 1:3) p <- p + .card_layers(xs[i], xs[i] + w, y0, y1, 10 + i, cols[i],
                                      shift = 0.2 * (i - 2))
  p <- p +
    .arrow(xs[1] + w / 2, y0 - 0.1, node[1] - 0.05, node[2] + 0.05, curv = 0.3, head = FALSE) +
    .arrow(xs[3] + w / 2, y0 - 0.1, node[1] + 0.05, node[2] + 0.05, curv = -0.3, head = FALSE) +
    .arrow(5, y0 - 0.1, 5, node[2] + 0.1, head = FALSE) +
    geom_point(aes(node[1], node[2]), colour = "white", size = 2.4) +
    .arrow(5, node[2] - 0.1, 5, 3.62) +
    .card_layers(3.0, 7.0, 0.3, 3.45, 33, gold, combined = TRUE, n = 110)
  .finish(p, c(0, 10.1), c(0.1, 7.6))
}

# Hub: three method volcanoes around a central combined volcano,
#    each pointing inwards
motif_hub <- function(cols = METHOD_COLS[c(1, 3, 2)], gold = "#F4B942") {
  ctr <- c(3.55, 6.45, 3.3, 6.5)            # xmin, xmax, ymin, ymax of combined
  p <- ggplot() +
    .card_layers(0.05, 2.75, 4.2, 6.2, 11, cols[1]) +
    .card_layers(7.25, 9.95, 4.2, 6.2, 12, cols[2]) +
    .card_layers(3.65, 6.35, 0.2, 2.2, 13, cols[3]) +
    .arrow(2.85, 5.2, ctr[1] - 0.12, 5.2) +
    .arrow(7.15, 5.2, ctr[2] + 0.12, 5.2) +
    .arrow(5, 2.3, 5, ctr[3] - 0.12) +
    .card_layers(ctr[1], ctr[2], ctr[3], ctr[4], 33, gold, combined = TRUE, n = 110)
  .finish(p, c(0, 10), c(0.1, 6.6))
}

# Deck: three method volcanoes stacked like cards ("assembled"),
#    then an arrow to the combined volcano
motif_deck <- function(cols = METHOD_COLS[c(2, 3, 1)], gold = "#F4B942") {
  w <- 3.7; h <- 2.6; off <- c(0.95, 0.75)
  p <- ggplot()
  for (i in 3:1) {                          # back to front
    x0 <- 0.1 + (i - 1) * off[1]; y0 <- 0.4 + (i - 1) * off[2]
    p <- p + .card_layers(x0, x0 + w, y0, y0 + h, 10 + i, cols[i])
  }
  p <- p +
    .arrow(0.1 + 2 * off[1] + w + 0.15, 2.45, 7.05, 2.45) +
    .card_layers(7.2, 10.0, 0.9, 4.0, 33, gold, combined = TRUE, n = 90)
  .finish(p, c(0, 10.1), c(0.2, 4.75))
}
