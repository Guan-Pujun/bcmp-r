source(".github/scripts/common.R")
record_environment()
pages <- list.files("man", pattern = "\\.Rd$", full.names = TRUE)
stopifnot(length(pages) > 0L)
for (page in pages) {
  rd <- tools::parse_Rd(page)
  issues <- tools::checkRd(rd)
  capture.output(print(issues), file = file.path(output, "Rd-diagnostics.txt"), append = TRUE)
  # Retain style diagnostics, but gate on parse issues, errors and warnings.
  problems <- capture.output(print(issues, minlevel = 0))
  if (length(problems)) stop("Rd issues in ", page, ": ", paste(problems, collapse = "\n"))
  tools::Rd2HTML(rd, out = file.path(output, paste0(basename(page), ".html")))
}
run_r(c("Rd2pdf", "--no-preview", paste0("--output=", file.path(output, "bcmp-manual.pdf")),
        root), "manual.log")
