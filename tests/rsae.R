library(rsae)

# functions that are required in the testing suite utility function
all_equal <- function(target, current, label,
                      tolerance = sqrt(.Machine$double.eps), scale = NULL,
                      check.attributes = FALSE)
{
    if (missing(label))
        stop("Argument 'label' is missing\n")
    res <- all.equal(target, current, tolerance, scale,
                     check.attributes = check.attributes)
    if (is.character(res))
        cat(paste0(label, ": ", res, "\n"))
}

data(landsat)
bhfmodel <- saemodel(formula = HACorn ~ PixelsCorn + PixelsSoybeans,
                     area = ~CountyName,
                     data = subset(landsat, subset = (outlier == FALSE)))

#-------------------------------------------------------------------------------
# maximum likelihood estimator
#-------------------------------------------------------------------------------
m <- fitsaemodel("ml", bhfmodel)

# check coefficients
ref_coef <- list(fixeff = structure(c(50.96755698661854, 0.32858050698728736,
    -0.13370990402151223), .Dim = c(1L, 3L), .Dimnames = list("fixeff",
    c("(Intercept)", "PixelsCorn", "PixelsSoybeans"))), raneff =
    structure(c(137.31355346825936, 121.06196596142701), .Dim = 1:2,
    .Dimnames = list("raneff", c("ResidualVar", "AreaVar"))))
all_equal(ref_coef, coef(m), label = "mle: coef")
# number of iterations
all_equal(1, m$converged, label = "mle: iterations")
# variance covariance matrix
ref_vcov <- structure(c(551.07980468822007, 20.130028525670017,
    13.462459120476971, -1.05455383869718, 0.0023024597019996413,
    3417.0328406402318, -1.1116340200125157, 0.0018440509470665526,
    0.0028156596100094009), .Dim = c(3L, 3L))
all_equal(ref_vcov, summary(m)$vcovbeta, label = "mle: vcov")
# prediction
d <- unique(landsat[-33, c("MeanPixelsCorn", "MeanPixelsSoybeans",
                           "CountyName")])
d <- cbind(rep(1,12), d); rownames(d) <- d$CountyName; d <- d[,1:3]
pr <- robpredict(m, areameans = d)
# prediction of fixed effects
ref_pred_fix <- structure(c(122.62932610201376, 123.37908865976928,
    118.67650271260092, 117.05345852448505, 130.37967556475763,
    102.42487738770292, 122.05168660183165, 120.35769616740677,
    104.07312831880023, 127.6710291896696, 121.73974193653029,
    134.40817795239562), .Dim = c(12L, 1L), .Dimnames = list(c("Cerro Gordo",
    "Hamilton", "Worth", "Humboldt", "Franklin", "Pocahontas", "Winnebago",
    "Wright", "Webster", "Hancock", "Kossuth", "Hardin"), NULL))
m_pred <- robpredict(m, d)
all_equal(ref_pred_fix, m_pred$fixeff, label = "mle: prediction: fixeff")
# prediction of random effects
ref_pred_raneff <- structure(c(-0.34795581561873956, 2.7306539845643276,
    -11.522100342616657, -8.3128316964299618, 13.641437929064862,
    9.5293713628629355, -9.043111288389424, 1.6482297982198533,
    11.082181455661823, -3.2293714931696194, -14.621108013778558,
    8.4446369759254374), .Dim = c(12L, 1L), .Dimnames = list(
    c("Cerro Gordo", "Hamilton", "Worth", "Humboldt", "Franklin",
    "Pocahontas", "Winnebago", "Wright", "Webster", "Hancock",
    "Kossuth", "Hardin"), NULL))
all_equal(ref_pred_raneff, m_pred$raneff, label = "mle: prediction: raneff")

#-------------------------------------------------------------------------------
# Huber M-estimator
#-------------------------------------------------------------------------------
m <- fitsaemodel("huberm", bhfmodel, k = 1.5)
# check coefficients
ref_coef <- list(fixeff = structure(c(50.284889482233901, 0.32845339396860779,
    -0.1391586558157463), .Dim = c(1L, 3L), .Dimnames = list("fixeff",
    c("(Intercept)", "PixelsCorn", "PixelsSoybeans"))), raneff =
    structure(c(152.8651163493895, 142.89195559015775), .Dim = 1:2,
    .Dimnames = list("raneff", c("ResidualVar", "AreaVar"))))
all_equal(ref_coef, coef(m), label = "huberm: coef")
# number of iterations
all_equal(1, m$converged, label = "huberm: iterations")
# variance covariance matrix
ref_vcov <- structure(c(617.34904107974603, 17.345477221076752,
    11.591770087408126, -1.1803136733091313, 0.0025778940009007264,
    2924.9265846687849, -1.2446871382399762, 0.0020639640773206873,
    0.0031559392715234426), .Dim = c(3L, 3L))
all_equal(ref_vcov, summary(m)$vcovbeta, label = "huberm: vcov")
# iterations
ref_optim <- structure(c(52.093820500254687, 55.259882759244917,
    50.733227002733955, 50.673399540505798, 50.3836345295719,
    50.320696700009591, 50.296102449042003, 50.288515631336672,
    50.285998734209528, 50.285149654076797, 50.284889482233901,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0.33358789961705571, 0.31834691276578125,
    0.32787363920213675, 0.32757047675141948, 0.32824617644909881,
    0.32837542034815403, 0.32842922081416226, 0.32844553684656014,
    0.32845098990884902, 0.32845283842511241, 0.32845339396860779,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, -0.15048524843533653, -0.15670431074104482,
    -0.14131898923730463, -0.14036253446157937, -0.13949260826911247,
    -0.13927474232364837, -0.13919537313766475, -0.13917044569139958,
    -0.1391622499352278, -0.13915950386737561, -0.1391586558157463,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 111.68227068477097, 144.14531496856631,
    148.78831340912174, 151.69589657133869, 152.46714498183184,
    152.73850847496223, 152.82404328674116, 152.85198329617418,
    152.86129812449576, 152.86424319274357, 152.8651163493895,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 1.2664842622335379,
    1.0643518882295631, 0.97117251506981217, 0.94665700291617505,
    0.93855257400764402, 0.9359886215259523, 0.93515434519498242,
    0.93488494813235012, 0.93479575384092883, 0.93476704272567468,
    0.93475842626883543, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0), .Dim = c(40L, 5L))
all_equal(ref_optim, attr(m, "optim")$tau, label = "huberm: optim")
# prediction
d <- unique(landsat[-33, c("MeanPixelsCorn", "MeanPixelsSoybeans",
                           "CountyName")])
d <- cbind(rep(1,12), d); rownames(d) <- d$CountyName; d <- d[,1:3]
pr <- robpredict(m, areameans = d)
# prediction of fixed effects
ref_pred_fix <- structure(c(120.87549517897703, 121.58673936423716,
    116.83850350968632, 115.1339100609233, 128.63186716427532,
    100.36297019739538, 120.32189621188968, 118.43059889784308,
    102.01080351347107, 125.86596357433179, 119.90424302449875,
    132.71937136988248), .Dim = c(12L, 1L), .Dimnames = list(c("Cerro Gordo",
    "Hamilton", "Worth", "Humboldt", "Franklin", "Pocahontas", "Winnebago",
    "Wright", "Webster", "Hancock", "Kossuth", "Hardin"), NULL))
m_pred <- robpredict(m, d)
all_equal(ref_pred_fix, m_pred$fixeff, label = "huberm: prediction: fixeff")
# prediction of random effects
ref_pred_raneff <- structure(c(0.17828694239225723, 4.7943357066320891,
    -13.97284439057621, -5.1623984071383386, 14.857085562188589,
    14.051391580853528, -10.295293125720921, 3.7889490911737225,
    16.550921841260163, -2.0976595103240858, -17.100553853635315,
    13.08070103738811), .Dim = c(12L, 1L), .Dimnames = list(c("Cerro Gordo",
    "Hamilton", "Worth", "Humboldt", "Franklin", "Pocahontas", "Winnebago",
    "Wright", "Webster", "Hancock", "Kossuth", "Hardin"), NULL))
all_equal(ref_pred_raneff, m_pred$raneff, label = "huberm: prediction: raneff")

#-------------------------------------------------------------------------------
# Huber M-estimator with initialization by least trimmed squares (lts)
# estimator in pkg robustbase
#-------------------------------------------------------------------------------
if (requireNamespace("robustbase", quietly = TRUE)) {
    m <- fitsaemodel("huberm", bhfmodel, k = 1.5, init = "lts")
    # check initialization
    ref_init <- c(45.978179340302198, 0.25190758887769144,
                  -0.034000084155002705, 231.41425594147111, 1)
    all_equal(ref_init, attr(m, "init"), "lts initialization of huberm")
    # check coefficients
    ref_init_coef <- list(fixeff = structure(c(50.284274527447955,
        0.32845491966527657, -0.13915705270413778), dim = c(1L, 3L),
        dimnames = list("fixeff", c("(Intercept)", "PixelsCorn",
        "PixelsSoybeans"))), raneff = structure(c(152.86658297422545,
        142.89116135291644), dim = 1:2, dimnames = list("raneff",
        c("ResidualVar", "AreaVar"))))
    all_equal(ref_init_coef, coef(m), "lts init: huberm: coef")
}
