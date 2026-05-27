## sim/intro_roadmap.R
##
## fig/intro_roadmap.pdf -- chapter-roadmap diagram for the
## Introduction.  Polished v2: consistent box sizing, cleaner
## arrow paths, more breathing room.

draw_box <- function(x, y, w, h, text,
                     fill = "white", border = "grey25",
                     text_cex = 0.95, text_font = 1,
                     header = NULL, header_cex = NULL,
                     header_font = 2, lwd = 1.3,
                     line_height = 0.27) {
  rect(x - w / 2, y - h / 2, x + w / 2, y + h / 2,
       col = fill, border = border, lwd = lwd)
  body_lines <- strsplit(text, "\n", fixed = TRUE)[[1L]]
  if (!is.null(header)) {
    text(x, y + h / 2 - 0.20, header,
         font = header_font,
         cex = if (is.null(header_cex)) text_cex * 1.05 else header_cex)
    y_top <- y + h / 2 - 0.20 - 0.30
  } else {
    n_lines <- length(body_lines)
    y_top   <- y + (n_lines - 1) / 2 * line_height
  }
  for (i in seq_along(body_lines)) {
    text(x, y_top - (i - 1) * line_height,
         body_lines[i], cex = text_cex, font = text_font)
  }
}

draw_arrow <- function(x0, y0, x1, y1, col = "grey30",
                       lwd = 1.4, length = 0.10) {
  arrows(x0, y0, x1, y1, length = length, angle = 22,
         col = col, lwd = lwd)
}

pdf("fig/intro_roadmap.pdf", width = 10.0, height = 6.8)
par(mar = c(0.4, 0.4, 0.4, 0.4))
plot.new()
plot.window(xlim = c(0, 10), ylim = c(0, 7.2), asp = NA)

## --- title ---
text(5.0, 7.00,
     "Empirical Likelihood under Censoring",
     font = 2, cex = 1.25)
text(5.0, 6.66,
     "Three censoring patterns  -  one Buckley-James thread  -  one R package",
     cex = 0.95, font = 3, col = "grey25")

## ============================================================
## TOP ROW -- three censoring patterns
## ============================================================
top_y <- 5.70
top_h <- 0.85
xpos  <- c(2.0, 5.0, 8.0)
top_labels <- list(
  "Right censoring\n(T_i, delta_i)",
  "Current-status /\nCase-1 interval-censoring",
  "Multinomial choice\n(rank order)"
)
for (i in seq_along(xpos)) {
  draw_box(xpos[i], top_y, 2.55, top_h, top_labels[[i]],
           fill = "grey94", text_font = 2, text_cex = 0.96)
}

## ============================================================
## MIDDLE ROW -- the methodological chapters
## ============================================================
mid_y <- 3.80
mid_h <- 1.55

draw_box(xpos[1], mid_y, 2.55, mid_h,
         text = "Recursive constrained\nKaplan-Meier estimator\n+ AFT regression EL",
         header = "Chapter 2 (KMC)",
         fill = "#dbeaf7", text_cex = 0.86, header_cex = 0.96)

draw_box(xpos[2], mid_y, 2.55, mid_h,
         text = "Wilks chi^2_{p-1}\nfor binary-choice EL\n+ centered (A7) variant",
         header = "Chapter 3 (BCM)",
         fill = "#fff3c0", text_cex = 0.86, header_cex = 0.96)

draw_box(xpos[3], mid_y, 2.55, mid_h,
         text = "Smoothed-rank EL\nchi^2_{p(J-1)-1}\nfor IIA + non-IIA",
         header = "Chapter 5 (multinomial)",
         fill = "#fcdada", text_cex = 0.86, header_cex = 0.96)

## arrows from censoring (top) to chapters (middle)
for (xx in xpos) {
  draw_arrow(xx, top_y - top_h / 2 - 0.05,
             xx, mid_y + mid_h / 2 + 0.05,
             col = "grey25")
}

## ============================================================
## BJ-RESIDUAL THREAD -- shaded band
## ============================================================
thread_y <- 2.50
rect(0.6, thread_y - 0.16, 9.4, thread_y + 0.16,
     col = "grey90", border = NA)
text(5.0, thread_y,
     "Shared thread: Buckley-James-style imputed residual",
     font = 3, cex = 0.92, col = "grey20")

## ============================================================
## CHAPTER 4 (AFT, high-dim) -- runs alongside the right-censored thread
## ============================================================
draw_box(xpos[1], 1.70, 2.55, 0.85,
         text = "Horowitz-style smoothing;\npost-2016 high-dim review",
         header = "Chapter 4 (high-dim AFT)",
         fill = "#e7ddf5", text_cex = 0.78, header_cex = 0.86)

## ============================================================
## BOTTOM -- the package
## ============================================================
draw_box(5.0, 0.55, 8.6, 0.70,
         text = "kmc.solve - kmc.solvelite - kmc.bjtest - kmc.bcm.test",
         header = "kmc R package v0.4-4 (CRAN)",
         fill = "#f4e2bf", text_cex = 0.86, header_cex = 0.96)

## arrows from each chapter to the package
draw_arrow(xpos[1], 1.70 - 0.85 / 2 - 0.05,
           4.0, 0.55 + 0.70 / 2 + 0.05, col = "grey30")
draw_arrow(xpos[2], thread_y - 0.16 - 0.05,
           5.0, 0.55 + 0.70 / 2 + 0.05, col = "grey30")
draw_arrow(xpos[3], thread_y - 0.16 - 0.05,
           6.0, 0.55 + 0.70 / 2 + 0.05, col = "grey30")

dev.off()
cat("Wrote fig/intro_roadmap.pdf (polished v2)\n")
