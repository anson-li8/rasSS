test_that("read_sumstats standardizes a Savage-like schema and filters chr", {
  tf <- tempfile(fileext = ".tsv")
  writeLines(c(
    "SNP\tCHR\tPOS\tA1\tA2\tEAF_HRC\tZscore\tstdBeta\tSE\tP",
    "rs1\t21\t15000000\tA\tG\t0.30\t2.5\t0.05\t0.02\t0.012",
    "rs2\t21\t15000100\tG\tA\t0.70\t-2.5\t-0.05\t0.02\t0.012",
    "rs3\t22\t15000200\tA\tC\t0.10\t1.0\t0.02\t0.02\t0.30"
  ), tf)
  dt <- read_sumstats(tf, chr = 21)
  expect_equal(nrow(dt), 2)
  expect_equal(names(dt), c("snp","chr","bp","a1","a2","eaf","z","beta","se","p"))
  expect_equal(dt$z, c(2.5, -2.5))
})

test_that("read_sumstats derives beta and z from OR/SE (PGC-like schema)", {
  tf <- tempfile(fileext = ".tsv")
  writeLines(c(
    "SNP\tCHR\tBP\tA1\tA2\tOR\tSE\tP",
    "rs9\t6\t32000000\tA\tG\t1.2\t0.1\t0.05"
  ), tf)
  dt <- read_sumstats(tf)
  expect_equal(dt$beta, log(1.2), tolerance = 1e-12)
  expect_equal(dt$z, log(1.2) / 0.1, tolerance = 1e-12)
})

test_that("harmonize flips Z and w when the counted allele is the other allele", {
  G <- matrix(c(0,1,2,1, 2,1,0,1), nrow = 4, ncol = 2)
  colnames(G) <- c("rs1","rs2")
  ld <- list(G = G, R = stats::cor(G), ids = colnames(G),
             freq = colMeans(G)/2, counted = c("A","C"))
  bim <- data.frame(id = c("rs1","rs2"), bp = c(100L, 200L),
                    a1 = c("A","C"), a2 = c("G","T"), stringsAsFactors = FALSE)
  ss  <- data.frame(snp = c("rs1","rs2"), chr = 1L, bp = c(100L, 200L),
                    a1 = c("A","T"), a2 = c("G","C"),
                    z = c(2, 2), beta = c(0.5, 0.5), stringsAsFactors = FALSE)
  h <- harmonize(ss, bim, ld)
  expect_equal(h$n_kept, 2)
  expect_equal(h$Z, c(2, -2))
  expect_equal(h$w, c(0.5, -0.5))
})

test_that("harmonize drops palindromic SNPs", {
  G <- matrix(c(0,1,2,1), nrow = 4, ncol = 1); colnames(G) <- "rs1"
  ld <- list(G = G, R = matrix(1, 1, 1), ids = "rs1",
             freq = colMeans(G)/2, counted = "A")
  bim <- data.frame(id = "rs1", bp = 100L, a1 = "A", a2 = "T", stringsAsFactors = FALSE)
  ss  <- data.frame(snp = "rs1", chr = 1L, bp = 100L, a1 = "A", a2 = "T",
                    z = 3, beta = 1, stringsAsFactors = FALSE)
  h <- harmonize(ss, bim, ld)
  expect_equal(h$n_kept, 0)
})