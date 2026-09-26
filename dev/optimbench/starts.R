# Stored starting points, by name. A cell's `start` field is `default:<seed>`
# (inits = NULL, set.seed(seed) immediately before ctFit: the user's path),
# `zeros`, or `stored:<name>` for one of these. Each entry is a function of
# (npar, rawnames) so a start can be defined by parameter name as well as by
# value; a stored vector whose length is not npar is refused, never recycled.
#
# Supplying inits switches the prior warm-up off (ctJuliaBackend.R), so a
# stored start is never the default path.
#
# Provenance of the vectors (all raw scale, the raw layout of juliaFit at
# 2026-09-24, npar 11 for families A-C, 9 for D, N and S):
#   anomS1_spurious  AnomAuth S1 'gated_zeros' fit; anomS2_spurious the S2
#                    'gated_seed2' fit -- the stored spurious Laplace maxima
#                    of review/LAPLACE-gated-gaps-2026-09-24.md section 3b
#                    (lapdef/M5_S1.rds, M5_S2.rds in session 00fd0b41).
#   hist_<config>    the best-known estimate for that config: the estimate
#                    with the highest penalised exact log likelihood that any
#                    floor reached from any start in the gaps sweeps M1 to M5
#                    (g2/res/final/out, softcut-3.5 reference). Its value is
#                    in references.csv.

BENCH_STARTS <- list(
  flatdrift8 = function(npar, rawnames) {
    # test-julia-convergence.R: drift inside its own transform's flat region.
    v <- rep(0, npar)
    hit <- which(rawnames == "drift")
    if (length(hit) != 1L) stop("flatdrift8 needs exactly one raw parameter named drift")
    v[hit] <- 8
    v
  },
  anomS1_spurious = function(npar, rawnames) .bench_fixed(c(
    c(0.2594856522638761, 5.3257409918419176, -6.0265586198841987, -1.1534740666333883, 
    0.0043894005909815466, -0.99595729559154345, 1.0201002138790745, -2.1416595716397064, 
    0.506005506505948)), npar),
  anomS2_spurious = function(npar, rawnames) .bench_fixed(c(
    c(0.26138160257507836, 5.6369503837343053, -4.8291200586487042, -1.1666413730068887, 
    0.0031223268411726206, -0.99073669620635796, 1.1091302032058568, -2.1557926624457751, 
    0.65897461428693715)), npar),
  hist_gA1 = function(npar, rawnames) .bench_fixed(c(
    c(0.11211357925405711, 1.1898202615963815, -6.4225211002935261, -1.0347363183679084, 
    0.24108548016395775, -0.47048921963525803, 0.84583910580834032, -0.34633266340415592, 
    0.78209870388749092, 0.97759830395515002, -0.051654882910698301)), npar),
  hist_gA10 = function(npar, rawnames) .bench_fixed(c(
    c(0.11049892358772781, 0.66121943451236831, -2.6676946800937524, -1.0563068070995463, 
    0.24985920048607935, -0.40261408400024978, 1.672741517625542, -0.19732460936350232, 
    0.10756440228513785, 0.60803188726199864, 0.11353953725529531)), npar),
  hist_gA14 = function(npar, rawnames) .bench_fixed(c(
    c(0.11273573878776515, 2.3580020302626266, -6.2899909740991902, -1.0408762719903606, 
    0.23690007139605662, -0.40776285981629756, 1.5238584841173688, -0.18699847654921764, 
    0.31981946682399748, 0.56119930386355577, 0.15382335595421368)), npar),
  hist_gB1 = function(npar, rawnames) .bench_fixed(c(
    c(0.20092881726473197, -0.76900023506252002, 0.93359208550670925, -4.2034272359772702, 
    0.2524385812096458, -0.62719269619845752, 0.16693001660611556, -0.56931670343952312, 
    0.72652269929563862, 0.40401823926194158, -0.52933204295252034)), npar),
  hist_gB2 = function(npar, rawnames) .bench_fixed(c(
    c(0.20389018455834682, 0.31610550728072179, 0.14849292977695891, -1.2628220705944142, 
    0.12321677744799002, -0.53031761915745212, 1.0181067107630748, -0.74434921451920311, 
    -0.15315102598306937, -0.16156265609640399, -0.17547597386158056)), npar),
  hist_gB8 = function(npar, rawnames) .bench_fixed(c(
    c(0.20557297297690003, 0.46746039952191787, -1.1443209943006718, -1.2445848027899815, 
    0.11490585224878765, -0.56825320314586658, 1.3005494339719654, -0.75053564477938661, 
    0.061670701411447482, -0.13243598163802134, 0.014207268078530389)), npar),
  hist_gC1 = function(npar, rawnames) .bench_fixed(c(
    c(0.20594548502855414, -0.30942045927055611, -1.14154522860593, -0.012533808577554702, 
    0.18092870063159344, -0.73257213510680419, 0.47729913926703393, -0.77429507409385623, 
    0.74291381723500682, 0.30431761924383643, -0.61591076957014623)), npar),
  hist_gC2 = function(npar, rawnames) .bench_fixed(c(
    c(0.19410425054642491, 0.10510879887419423, -5.9866542641041747, 0.31184364986377588, 
    0.14919305200537197, -0.44308784591881273, 0.47264722470263604, -0.62911310530446951, 
    -0.21864359757380827, -0.15103551871992466, -0.33081233796123816)), npar),
  hist_gC8 = function(npar, rawnames) .bench_fixed(c(
    c(0.20384085302070734, 0.27895427632453201, -1.4074921644370737, -0.67037533429960061, 
    0.12961683919096909, -0.60972888042909557, 1.0256422360220598, -0.70363457178177291, 
    0.01510086133959457, -0.099707609228174537, -0.074771305679053726)), npar),
  hist_gD1 = function(npar, rawnames) .bench_fixed(c(
    c(0.03066698976194723, 1.4430236323139731, -1.0983004726276682, 0.023393159172322594, 
    -0.063160010359700841, -0.0033635939395050764, -4.1543634957057334, -0.78650708044610629, 
    -4.1428419699078826e-06)), npar),
  hist_gD3 = function(npar, rawnames) .bench_fixed(c(
    c(0.0086006842520587636, 1.9646495552686818, -1.1235621307267956, 0.049952819920549348, 
    -0.049124895452136663, 0.012548034281602193, 1.0712705121988109, -0.85260496687947995, 
    -0.25564134755647872)), npar),
  hist_gN1 = function(npar, rawnames) .bench_fixed(c(
    c(0.0027101708406166948, 1.0444046947168477, -1.4174100098355282, -1.0650846683102619, 
    0.095313030797729778, 0.66020154741091264, -0.75089683950360608, 0.098183915307533773, 
    -0.94556931674601441)), npar),
  hist_gN2 = function(npar, rawnames) .bench_fixed(c(
    c(0.019128203266238316, 1.0684132639300874, -1.5531756554504557, -1.0286111062364121, 
    0.22228457814177663, 0.98813806434605522, -0.74298592046876055, 0.25885483672415527, 
    -0.86117819635117954)), npar),
  hist_gN3 = function(npar, rawnames) .bench_fixed(c(
    c(0.0097153515570728918, 0.52790307222832467, -1.4343436193645438, -1.1127294644570047, 
    0.1826321191014787, 0.93541604232839504, -0.52454924307136508, 0.47013480002368918, 
    -0.73719314807595293)), npar),
  hist_gN4 = function(npar, rawnames) .bench_fixed(c(
    c(-0.00019522346110405003, 1.265836328627898, -1.4449780754499084, -1.072623533610545, 
    0.20668858864472064, 0.72543396096273927, -0.57715176395926127, 0.56965760969365886, 
    -0.84169807974452915)), npar),
  hist_anomS1 = function(npar, rawnames) .bench_fixed(c(
    c(0.25841485579787438, 3.9344502667419685, -2.0669563608701198, -1.1697546347293633, 
    0.0070872082449264292, -0.98701235022892519, -4.4929536447377387, -6.6331247233414095, 
    1.3231526865941474e-06)), npar),
  hist_anomS2 = function(npar, rawnames) .bench_fixed(c(
    c(0.26051551402977685, 4.631709070752807, -2.1001569771203314, -1.1807685800451622, 
    0.0036940185321996453, -0.98858955673283333, -4.4734710616139584, -6.9161636212329292, 
    6.4754681921140952e-06)), npar)
)

.bench_fixed <- function(v, npar) {
  if (length(v) != npar) stop("stored start has ", length(v), " values; the model has ",
    npar, " raw parameters")
  v
}

bench_start <- function(start, npar, rawnames) {
  if (identical(start, "zeros")) return(list(inits = rep(0, npar), seed = 1L))
  if (grepl("^default:[0-9]+$", start))
    return(list(inits = NULL, seed = as.integer(sub("^default:", "", start))))
  if (grepl("^stored:", start)) {
    nm <- sub("^stored:", "", start)
    f <- BENCH_STARTS[[nm]]
    if (is.null(f)) stop("no stored start named ", nm)
    return(list(inits = f(npar, rawnames), seed = 1L))
  }
  stop("start must be default:<seed>, zeros or stored:<name>, not ", start)
}
