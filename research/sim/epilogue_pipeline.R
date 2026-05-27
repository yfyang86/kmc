## sim/epilogue_pipeline.R
##
## fig/epilogue_pipeline.pdf -- a 3-stage pipeline figure for
## the Epilogue's "integrated open direction":
##
##   1. AFT-SDAR / broken adaptive ridge -> active set hat S
##           (contribution of post-2016 high-dim AFT literature)
##   2. kmc.bjtest restricted to hat S -> chi^2_{|hat S|} test
##           (contribution of Chapter 2, KMC)
##   3. Centering correction for post-selection bias
##           (contribution of Chapter 3, BCM)

draw_box <- function(x, y, w, h, header, body,
                     fill = "white", border = "black",
                     header_cex = 1.0, body_cex = 0.83, lwd = 1.6) {
  rect(x - w/2, y - h/2, x + w/2, y + h/2,
       col = fill, border = border, lwd = lwd)
  text(x, y + h/2 - 0.18, header,
       font = 2, cex = header_cex)
  ## body: multi-line, justified
  lines <- strsplit(body, "\n", fixed = TRUE)[[1L]]
  for (i in seq_along(lines)) {
    text(x, y + h/2 - 0.43 - (i - 1) * 0.20,
         lines[i], cex = body_cex)
  }
}

big_arrow <- function(x0, y0, x1, y1, col = "grey25",
                      lwd = 2.0, length = 0.20) {
  arrows(x0, y0, x1, y1, length = length, angle = 22,
         col = col, lwd = lwd)
}

pdf("fig/epilogue_pipeline.pdf", width = 11.0, height = 5.0)
par(mar = c(0.3, 0.3, 0.3, 0.3))
plot.new()
plot.window(xlim = c(0, 11), ylim = c(0, 5.2))

## --- title ---
text(5.5, 4.95,
     "The integrated open direction",
     font = 2, cex = 1.22)
text(5.5, 4.60,
     "AFT-SDAR active-set selection  ->  kmc.bjtest EL inference  ->  BCM centering correction",
     font = 3, cex = 0.92, col = "grey30")

## --- the three boxes ---
box_y <- 2.4
box_h <- 2.10
box_w <- 3.10

draw_box(1.85, box_y, box_w, box_h,
         "(1) AFT-SDAR",
         "Estimate active set\nS-hat in {1,...,p}\nvia L0-regularised\nhigh-dim AFT.",
         fill = "lavender", header_cex = 1.05)

draw_box(5.50, box_y, box_w, box_h,
         "(2) kmc.bjtest",
         "Restrict to columns\nS-hat;  run Owen-EL\nwith BJ residuals.\n-2 log ELR ~ chi^2_{|S-hat|}.",
         fill = "lightblue", header_cex = 1.05)

draw_box(9.15, box_y, box_w, box_h,
         "(3) Centering correction",
         "Apply tilde g = (X - hat m) eps\nto discharge the\nno-bias condition\npost-selection.",
         fill = "lightyellow", header_cex = 1.05)

## --- arrows connecting them ---
big_arrow(1.85 + box_w/2 + 0.07, box_y,
          5.50 - box_w/2 - 0.07, box_y)
big_arrow(5.50 + box_w/2 + 0.07, box_y,
          9.15 - box_w/2 - 0.07, box_y)

## --- contribution attribution underneath each box ---
contrib_y <- 0.85
text(1.85, contrib_y, "(post-2016 AFT literature:",
     font = 3, cex = 0.80, col = "grey25")
text(1.85, contrib_y - 0.30,
     "Feng et al.\\ '22, Reeder '24)",
     font = 3, cex = 0.80, col = "grey25")
text(5.50, contrib_y, "(Chapter 2 of this book)",
     font = 3, cex = 0.85, col = "grey25")
text(9.15, contrib_y, "(Chapter 3 of this book,",
     font = 3, cex = 0.80, col = "grey25")
text(9.15, contrib_y - 0.30,
     "Theorem 3.10)",
     font = 3, cex = 0.80, col = "grey25")

## --- bottom summary line ---
text(5.5, 0.10,
     "Each step is now individually in place; the integration is open work.",
     font = 3, cex = 0.92, col = "grey20")

dev.off()
cat("Wrote fig/epilogue_pipeline.pdf\n")
