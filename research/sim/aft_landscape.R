## sim/aft_landscape.R
##
## fig/aft_landscape.pdf -- a categorical map of regularized-AFT
## methods organised by:
##   x-axis: era (pre-2016  |  post-2016)
##   y-axis: methodological family (loss type)
##
## Each method is a small box positioned in its cell; the
## point is to give the reader a one-glance view of how the
## field shifted between Yang's 2016 dissertation and the
## 2024-26 expansion.

draw_box <- function(x, y, w, h, text, fill = "white",
                     border = "black", text_cex = 0.85,
                     text_font = 1, lwd = 1.0) {
  rect(x - w/2, y - h/2, x + w/2, y + h/2,
       col = fill, border = border, lwd = lwd)
  ## handle multi-line text
  lines <- strsplit(text, "\n", fixed = TRUE)[[1L]]
  n_lines <- length(lines)
  y_top <- y + (n_lines - 1) / 2 * 0.10
  for (i in seq_along(lines)) {
    text(x, y_top - (i - 1) * 0.10,
         lines[i], cex = text_cex, font = text_font)
  }
}

pdf("fig/aft_landscape.pdf", width = 10.0, height = 6.6)
par(mar = c(0.4, 0.4, 0.4, 0.4))
plot.new()
plot.window(xlim = c(0, 10), ylim = c(0, 7.2))

## --- title ---
text(5.0, 7.0,
     "Regularized AFT: pre-2016 vs post-2016 landscape",
     font = 2, cex = 1.18)

## --- column headers ---
rect(0.5, 5.95, 4.7, 6.45, col = "grey92", border = "grey60", lwd = 0.8)
text(2.6, 6.20, "pre-2016 (dissertation era)", font = 4, cex = 1.00)
rect(5.3, 5.95, 9.5, 6.45, col = "grey92", border = "grey60", lwd = 0.8)
text(7.4, 6.20, "post-2016 (book-update era)", font = 4, cex = 1.00)

## --- row labels (left axis) ---
row_ys     <- c(5.3, 4.2, 3.1, 2.0, 0.9)
row_labels <- c("Penalty form",
                "Loss / score",
                "Variable-selection\nstrategy",
                "Inference\nmachinery",
                "Computation")
for (i in seq_along(row_ys)) {
  text(0.18, row_ys[i], row_labels[i],
       font = 2, cex = 0.84, adj = 0)
}

## --- horizontal separator lines ---
for (yy in c(5.90, 4.75, 3.65, 2.55, 1.45, 0.35)) {
  segments(0.5, yy, 9.5, yy, col = "grey80", lwd = 0.6, lty = 3)
}
## vertical column separator
segments(5.0, 0.35, 5.0, 5.90, col = "grey60", lwd = 0.6)

## --- entries: pre-2016 ---
draw_box(2.0, row_ys[1], 1.8, 0.46,
         "LASSO ('96)",
         fill = "lightyellow")
draw_box(3.7, row_ys[1], 1.8, 0.46,
         "SCAD ('01)\nAdaptive LASSO\n('06)", fill = "lightyellow", text_cex = 0.72)

draw_box(2.0, row_ys[2], 1.8, 0.46,
         "Buckley-James\n(KM-imputed LS)", fill = "lavender", text_cex = 0.75)
draw_box(3.7, row_ys[2], 1.8, 0.46,
         "log-rank / Gehan\n(Jin et al. '03)", fill = "lavender", text_cex = 0.72)

draw_box(2.0, row_ys[3], 1.8, 0.46,
         "stepwise + CV", fill = "lightblue", text_cex = 0.80)
draw_box(3.7, row_ys[3], 1.8, 0.46,
         "regularised BJ\n(Huang et al. '06)", fill = "lightblue", text_cex = 0.72)

draw_box(2.0, row_ys[4], 1.8, 0.46,
         "Wald + plug-in\nvariance", fill = "mistyrose", text_cex = 0.72)
draw_box(3.7, row_ys[4], 1.8, 0.46,
         "BJ-EL\n(Zhou-Li '08)", fill = "mistyrose", text_cex = 0.75)

draw_box(2.0, row_ys[5], 1.8, 0.46,
         "EM-EL\n(Zhou '05)", fill = "wheat", text_cex = 0.78)
draw_box(3.7, row_ys[5], 1.8, 0.46,
         "KMC recursion\n(Zhou-Yang '15)", fill = "wheat", text_cex = 0.72)

## --- entries: post-2016 ---
draw_box(6.6, row_ys[1], 1.8, 0.46,
         "MCP ('10)\nL0-AFT\n(Feng et al. '22)", fill = "lightyellow", text_cex = 0.72)
draw_box(8.4, row_ys[1], 1.8, 0.46,
         "Bayesian group\nlasso\n(Reeder '24)", fill = "lightyellow", text_cex = 0.72)

draw_box(6.6, row_ys[2], 1.8, 0.46,
         "DNN AFT\n(Norman '24)", fill = "lavender", text_cex = 0.75)
draw_box(8.4, row_ys[2], 1.8, 0.46,
         "smoothed-rank DNN\n(Kim '24)", fill = "lavender", text_cex = 0.72)

draw_box(6.6, row_ys[3], 1.8, 0.46,
         "SDAR / ASDAR\n(Feng et al. '22)", fill = "lightblue", text_cex = 0.72)
draw_box(8.4, row_ys[3], 1.8, 0.46,
         "broken adaptive\nridge (Lee '24)", fill = "lightblue", text_cex = 0.72)

draw_box(6.6, row_ys[4], 1.8, 0.46,
         "kmc.bjtest\n(v0.4-3, this book)", fill = "mistyrose", text_cex = 0.72)
draw_box(8.4, row_ys[4], 1.8, 0.46,
         "post-selection EL\n(open)", fill = "mistyrose", text_cex = 0.72, text_font = 3)

draw_box(6.6, row_ys[5], 1.8, 0.46,
         "kmc.solvelite\n(v0.4-3)", fill = "wheat", text_cex = 0.72)
draw_box(8.4, row_ys[5], 1.8, 0.46,
         "FFT smoother /\nGPU (open)", fill = "wheat", text_cex = 0.72, text_font = 3)

## --- bottom footnote ---
text(5.0, 0.10,
     paste0("Italic boxes mark open directions called out in Section 4.5 / Epilogue."),
     cex = 0.75, font = 3, col = "grey30")

dev.off()
cat("Wrote fig/aft_landscape.pdf\n")
