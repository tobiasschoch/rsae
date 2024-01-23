library("rsae")
library("microbenchmark")

test_fun <- function(reps = 20, n = 4, g = 20, seed = 1024, ve.epsilon = 0,
                     vu.epsilon = 0, estimator = "ml", ...)
{
    set.seed(seed)
    for (i in 1:reps) {
        dat <- makedata(n = n, g = g, ve.epsilon = ve.epsilon,
                        vu.epsilon = vu.epsilon)
        fitsaemodel(estimator, dat, ...)
    }
}

# Small sample: mle and m-est
s1_ml <- microbenchmark(test_fun(), times = 20L)
s1_hub <- microbenchmark(test_fun(estimator = "huberm", k = 1.5),
               times = 20L)
# Medium sample: mle and m-est
m1_ml <- microbenchmark(test_fun(n = 20, g = 30), times = 10L)
m1_hub <- microbenchmark(test_fun(n = 20, g = 30, estimator = "huberm",
                                  k = 1.5), times = 10L)
# Large sample: mle and m-est
l1_ml <- microbenchmark(test_fun(reps = 10, n = 50, g = 50), times = 10L)
l1_hub <- microbenchmark(test_fun(reps = 10, n = 50, g = 50,
                                  estimator = "huberm", k = 1.5), times = 10L)
# Large sample: mle and m-est with contamination: ve
l2_ml <- microbenchmark(test_fun(reps = 10, n = 50, g = 50, ve.epsilon = 0.05),
               times = 10L)
l2_hub <- microbenchmark(test_fun(reps = 10, n = 50, g = 50, ve.epsilon = 0.05,
                        estimator = "huberm", k = 1.5), times = 10L)
# Large sample: mle and m-est with contamination: vu
l3_ml <- microbenchmark(test_fun(reps = 10, n = 50, g = 50, vu.epsilon = 0.05),
               times = 10L)
l3_hub <- microbenchmark(test_fun(reps = 10, n = 50, g = 50, vu.epsilon = 0.05,
                        estimator = "huberm", k = 1.5), times = 10L)

ml_avg <- NULL
for (v in ls(pattern = "_ml", sorted = TRUE)) {
    ml_avg <- c(ml_avg, mean(get(v)$time / 1e6))
}
ml_q75 <- NULL
for (v in ls(pattern = "_ml", sorted = TRUE)) {
    ml_q75 <- c(ml_q75, quantile(get(v)$time / 1e6, probs = 0.75))
}
hub_avg <- NULL
for (v in ls(pattern = "_hub", sorted = TRUE)) {
    hub_avg <- c(hub_avg, mean(get(v)$time / 1e6))
}
hub_q75 <- NULL
for (v in ls(pattern = "_hub", sorted = TRUE)) {
    hub_q75 <- c(hub_q75, quantile(get(v)$time / 1e6, probs = 0.75))
}

res <- cbind(ml_avg, ml_q75, hub_avg, hub_q75)
rownames(res) <- c("l1", "l2", "l3", "m1", "s1")
the_name <- paste0("benchmark_", packageVersion("rsae"))
assign(the_name, res)
save(list = the_name, file = paste0(the_name, ".RData"))

load("benchmark_0.2.RData")
cat("Changes in % (positive values: slower)\n")
print(round(100 * (benchmark_0.3 / benchmark_0.2 - 1), 2))
