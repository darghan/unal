library(shiny)
library(DT)
library(ggplot2)

# ══════════════════════════════════════════════════════════════
# Helper functions
# ══════════════════════════════════════════════════════════════

# Build the design: t = 3 (C, T1, T2) on a fixed 9x9 field grid (81 cells,
# 27 replicates per treatment). The grid only defines positions and the
# neighbourhood structure; treatments are assigned completely at random
# (completely randomised design, see below).
build_design <- function(neighbor_type = c("rook", "queen"), seed_treat = 1) {
  neighbor_type <- match.arg(neighbor_type)
  grid_size <- 9
  n <- grid_size * grid_size  # 81
  
  cell_row <- rep(1:grid_size, times = grid_size)
  cell_col <- rep(1:grid_size, each = grid_size)
  
  # True CRD: treatments are randomly assigned across all 81 cells (27 per
  # level), NOT arranged in a systematic spatial block/pattern -- systematic
  # placement would confound treatment with spatial position and violate
  # complete randomisation. Any spatial clustering the reader perceives must
  # come from the response Y itself (through the true active-treatment lag
  # theta_act), not from the treatment layout.
  set.seed(seed_treat)
  levels_vec <- sample(rep(c("C", "T1", "T2"), each = n / 3))
  
  Xtau <- cbind(
    as.numeric(levels_vec == "C"),
    as.numeric(levels_vec == "T1"),
    as.numeric(levels_vec == "T2")
  )
  Xact <- Xtau[, 2:3, drop = FALSE]
  
  J <- matrix(1, n, 1)
  X <- cbind(J, Xtau)  # n x 4, rank-DEFICIENT: rho(X) = 3, since Xtau %*% c(1,1,1) = J
  
  # --- SVD-based reparametrisation X -> Z (full column rank), as in the
  #     original companion note. Vr collects the right singular vectors
  #     associated with the non-zero singular values of X. ---
  sv <- svd(X)
  tol <- max(dim(X)) * max(sv$d) * .Machine$double.eps
  r_rank <- sum(sv$d > tol)              # should be 3
  Vr <- sv$v[, 1:r_rank, drop = FALSE]   # 4 x 3 transformation matrix
  Z <- X %*% Vr                          # n x 3, full column rank, C(Z) = C(X)
  
  MZ <- Z %*% solve(t(Z) %*% Z) %*% t(Z)  # ppm onto C(Z) = C(X)
  
  # --- Spatial weights matrix (rook or queen contiguity), row-standardised ---
  row_diff <- outer(cell_row, cell_row, "-")
  col_diff <- outer(cell_col, cell_col, "-")
  if (neighbor_type == "rook") {
    Adj <- (abs(row_diff) + abs(col_diff)) == 1
  } else {
    Adj <- (abs(row_diff) <= 1 & abs(col_diff) <= 1) & !(row_diff == 0 & col_diff == 0)
  }
  Adj <- matrix(as.numeric(Adj), n, n)
  W <- sweep(Adj, 1, rowSums(Adj), FUN = "/")  # row-standardised: W %*% 1_n = 1_n
  W_rank <- qr(W)$rank  # informational only: S3 (W non-singular) is a SUFFICIENT,
  # not necessary, route to Proposition 1; a contiguity W on a
  # bipartite lattice graph can be rank-deficient without that
  # breaking estimability of the active-treatment lag itself.
  
  XWact <- W %*% Xact
  
  I_n <- diag(n)
  A <- (I_n - MZ) %*% XWact          # the actual object Condition S4 / Proposition 2 concern
  AtA <- t(A) %*% A                  # 2x2 -- the coefficients that verify Condition S4
  A_rank <- qr(A)$rank
  det_AtA <- det(AtA)
  
  if (A_rank < 2) {
    warning("Rank(A) < 2: Condition S4 fails for this W / treatment layout -- ",
            "the active-treatment lag is not fully estimable with this design.")
  }
  
  Mlag <- A %*% solve(AtA) %*% t(A)
  ML <- MZ + Mlag
  
  # Weak (rank) sufficiency, already verified above via rho(A) and det(A'A),
  # guarantees the correct degrees of freedom. It does NOT guarantee that
  # X_{W,act} is exactly orthogonal to C(Z); that stronger property is what
  # keeps the Treatment-effect expected mean square free of any contribution
  # from theta_act. This single ratio reports how far the design is from
  # that stronger, idealised condition: 0 = fully orthogonal (ideal),
  # closer to 1 = X_{W,act} is largely explained by C(Z) already.
  leak_Z <- norm(MZ %*% XWact, type = "F") / norm(XWact, type = "F")
  
  list(
    grid_size = grid_size, n = n, t_levels = 3, r = n / 3,
    cell_row = cell_row, cell_col = cell_col, levels_vec = levels_vec,
    X = X, sv = sv, r_rank = r_rank, Vr = Vr, Z = Z,
    Xtau = Xtau, Xact = Xact, XWact = XWact, W = W, W_rank = W_rank,
    neighbor_type = neighbor_type,
    J = J, MJ = (J %*% t(J)) / n, MZ = MZ,
    A = A, AtA = AtA, A_rank = A_rank, det_AtA = det_AtA,
    Mlag = Mlag, ML = ML, I_n = I_n, leak_Z = leak_Z,
    df_treat = 3 - 1, df_lag = 2, df_err = n - 5, df_WL = n - 3
  )
}

# Generate one realised response vector Y under given true parameters
generate_Y <- function(design, tau, theta_act, sigma, seed_y) {
  set.seed(seed_y)
  mu <- as.vector(design$Xtau %*% tau + design$XWact %*% theta_act)
  Y <- mu + rnorm(design$n, 0, sigma)
  Y
}

# The usual (lag-ignoring) ANOVA table
anova_usual <- function(design, Y) {
  SSt  <- as.numeric(t(Y) %*% (design$MZ - design$MJ) %*% Y)
  SSwl <- as.numeric(t(Y) %*% (design$I_n - design$MZ) %*% Y)
  MSt  <- SSt / design$df_treat
  MSwl <- SSwl / design$df_WL
  Fwl  <- MSt / MSwl
  pwl  <- pf(Fwl, design$df_treat, design$df_WL, lower.tail = FALSE)
  SStot <- SSt + SSwl
  
  data.frame(
    Source = c("Treatment", "Error (unlagged)", "Total (corrected)"),
    df = c(design$df_treat, design$df_WL, design$df_treat + design$df_WL),
    SS = c(SSt, SSwl, SStot),
    MS = c(MSt, MSwl, NA),
    F = c(Fwl, NA, NA),
    p_value = c(pwl, NA, NA)
  )
}

# The full spatial (SAV) ANOVA table, active-treatment-only lag
anova_spatial <- function(design, Y) {
  SSmu <- as.numeric(t(Y) %*% design$MJ %*% Y)
  SSt  <- as.numeric(t(Y) %*% (design$MZ - design$MJ) %*% Y)
  SSl  <- as.numeric(t(Y) %*% design$Mlag %*% Y)
  SSe  <- as.numeric(t(Y) %*% (design$I_n - design$ML) %*% Y)
  SStot <- SSt + SSl + SSe
  
  MSt <- SSt / design$df_treat
  MSl <- SSl / design$df_lag
  MSe <- SSe / design$df_err
  
  Ft <- MSt / MSe
  Fl <- MSl / MSe
  pt <- pf(Ft, design$df_treat, design$df_err, lower.tail = FALSE)
  pl <- pf(Fl, design$df_lag, design$df_err, lower.tail = FALSE)
  
  data.frame(
    Source = c("Mean", "Treatment", "Spatial lag (active)", "Error", "Total (corrected)"),
    Halmos_term = c("P0", "P1 - P0", "P2 - P1", "I - P2", "I - P0"),
    Subspace = c("V0 = C(J)", "V1 \u2229 V0\u22a5", "V2 \u2229 V1\u22a5", "V2\u22a5", "V0\u22a5"),
    df = c(1, design$df_treat, design$df_lag, design$df_err,
           design$df_treat + design$df_lag + design$df_err),
    SS = c(SSmu, SSt, SSl, SSe, SStot),
    MS = c(NA, MSt, MSl, MSe, NA),
    F = c(NA, Ft, Fl, NA, NA),
    p_value = c(NA, pt, pl, NA, NA)
  )
}

# Monte Carlo validation grid (Section 3.1.2 of the paper):
# for each lag magnitude c, generates B replications under
# theta_act = c(1, 1)' and the given tau, with sigma^2 = 1, and returns
# empirical rejection rates plus the bias check of the lag-ignoring residual
# mean square, E(MS_WL) = sigma^2 + phi_lag/(n - t) (Section 2.4.1).
# The row c = 0 is the complete null hypothesis (empirical size).
run_grid <- function(design, tau, cs, B, seed = 123) {
  set.seed(seed)
  n <- design$n
  Xtau <- design$Xtau; XWact <- design$XWact; A <- design$A
  MJ <- design$MJ; MZ <- design$MZ; Mlag <- design$Mlag; ML <- design$ML
  I_n <- design$I_n
  df_treat <- design$df_treat; df_lag <- design$df_lag
  df_err <- design$df_err; df_WL <- design$df_WL
  
  crit_treat    <- qf(0.95, df_treat, df_err)
  crit_lag      <- qf(0.95, df_lag, df_err)
  crit_treat_WL <- qf(0.95, df_treat, df_WL)
  
  qform <- function(Ymat, P) rowSums((Ymat %*% P) * Ymat)
  mu_treat <- as.vector(Xtau %*% tau)
  
  rows <- lapply(cs, function(c) {
    theta_act <- c(c, c)
    mu <- mu_treat + as.vector(XWact %*% theta_act)
    eps <- matrix(rnorm(B * n, 0, 1), nrow = B, ncol = n)
    Ymat <- sweep(eps, 2, mu, "+")
    
    SSt  <- qform(Ymat, MZ - MJ)
    SSl  <- qform(Ymat, Mlag)
    SSe  <- qform(Ymat, I_n - ML)
    SSwl <- qform(Ymat, I_n - MZ)
    
    F_treat    <- (SSt / df_treat) / (SSe / df_err)
    F_lag      <- (SSl / df_lag) / (SSe / df_err)
    F_treat_WL <- (SSt / df_treat) / (SSwl / df_WL)
    
    phiL <- as.numeric(t(theta_act) %*% (t(A) %*% A) %*% theta_act)
    theory_wl <- 1 + phiL / (n - 3)
    
    data.frame(
      c = c,
      rej_treat_correct = mean(F_treat > crit_treat),
      rej_treat_naive = mean(F_treat_WL > crit_treat_WL),
      rej_lag = mean(F_lag > crit_lag),
      phiL_theory = phiL,
      EMS_E_emp = mean(SSe / df_err),
      EMS_WL_emp = mean(SSwl / df_WL),
      EMS_WL_theory = theory_wl
    )
  })
  do.call(rbind, rows)
}

# Moran's I for a residual vector e, given spatial weights W.
# I = (n/S0) * (e'We)/(e'e); since W is row-standardised here, S0 = n, so I = e'We/e'e.
moran_I <- function(e, W) {
  n_e <- length(e)
  S0 <- sum(W)
  as.numeric((n_e / S0) * (t(e) %*% W %*% e) / (t(e) %*% e))
}

# Permutation (Monte Carlo) significance test for Moran's I: shuffles residuals
# across locations (keeping W fixed) to build the null distribution.
# One-sided, H1: positive spatial autocorrelation.
moran_test <- function(e, W, nsim = 999) {
  I_obs <- moran_I(e, W)
  I_perm <- replicate(nsim, moran_I(sample(e), W))
  p_val <- (sum(I_perm >= I_obs) + 1) / (nsim + 1)
  list(I = I_obs, p.value = p_val)
}

# Render the grid as a ready-to-paste LaTeX table (booktabs style)
grid_to_latex <- function(df, B, seed_mc, seed_treat = 1) {
  header <- paste0(
    "\\begin{table}[h]\n\\centering\n",
    "\\caption{Monte Carlo validation grid for the active-treatment-only ",
    "spatial ANOVA (B = ", B, " replications per row; rook contiguity; ",
    "treatment-randomization seed = ", seed_treat, "; simulation seed = ", seed_mc,
    "; $\\tau = \\mathbf{0}$; $\\sigma^2 = 1$).}\n",
    "\\begin{tabular}{@{}rrrrrrrr@{}}\n\\toprule\n",
    "$c$ & Rej.\\ Treat.\\ (correct) & Rej.\\ Treat.\\ (naive) & Rej.\\ Lag & ",
    "$\\phi_L$ & $E[MS_E]$ emp. & $E[MS_{WL}]$ emp. & $E[MS_{WL}]$ theory \\\\\n",
    "\\midrule\n"
  )
  body <- paste0(
    sprintf(
      "%.1f & %.4f & %.4f & %.4f & %.3f & %.4f & %.4f & %.4f \\\\\n",
      df$c, df$rej_treat_correct, df$rej_treat_naive, df$rej_lag,
      df$phiL_theory, df$EMS_E_emp, df$EMS_WL_emp, df$EMS_WL_theory
    ),
    collapse = ""
  )
  footer <- "\\bottomrule\n\\end{tabular}\n\\end{table}\n"
  paste0(header, body, footer)
}

# ══════════════════════════════════════════════════════════════
# UI
# ══════════════════════════════════════════════════════════════

ui <- fluidPage(
  tags$head(tags$style(HTML("
    body { font-family: 'Segoe UI', sans-serif; font-size: 14px; color: #2c3e50; background-color: #f8f9fa; }
    .well { background-color: #ffffff; border: 1px solid #e1e8ed; border-radius: 12px;
            box-shadow: 0 2px 8px rgba(0,0,0,0.05); padding: 20px; }
    h2, h3, h4 { font-weight: 600; color: #1a1a2e; }
    .nav-tabs>li>a { font-weight: 600; color: #7f8c8d; }
    .nav-tabs>li.active>a { color: #3498db; }
    .tab-content { background-color: #ffffff; padding: 22px; border-radius: 12px;
                   margin-top: 15px; box-shadow: 0 2px 12px rgba(0,0,0,0.08); }
    pre.latex-box { background:#f4f4f4; padding:14px; border-radius:8px; font-size:12px;
                    white-space: pre-wrap; word-break: break-word; }
    .btn-primary { background-color:#2980b9; border-color:#2980b9; }
  "))),
  
  div(style = "max-width:1300px; margin:0 auto; padding:20px;",
      h2("Active-Treatment-Only Spatial ANOVA \u2014 Design, ANOVA Tables, and Monte Carlo Grid"),
      p(style = "color:#7f8c8d;",
        "Companion Shiny app to the one-way SAV note with a non-lagged control (t = 3: Control, T1, T2)."),
      
      sidebarLayout(
        sidebarPanel(
          width = 3,
          h4("Design"),
          p(style = "color:#7f8c8d; font-size:12px;",
            "Fixed 9\u00d79 field grid (81 cells), used only to define positions and neighbourhood. Treatments are assigned at random (completely randomized design, 27 replicates per level)."),
          numericInput("seed_treat", "Treatment randomization seed", value = 1, step = 1),
          selectInput("neighbor_type", "Neighbourhood type (W)",
                      choices = c("Rook (4 neighbours)" = "rook", "Queen (8 neighbours)" = "queen"),
                      selected = "queen"),
          numericInput("sigma", "Error standard deviation (\u03c3)", value = 0.31, min = 0.01, step = 0.05),
          hr(),
          h4("True parameters (to generate Y)"),
          numericInput("tauC", "\u03c4 Control (baseline mean)", value = 3.0, step = 0.1),
          numericInput("tauT1", "\u03c4 T1", value = 3.2, step = 0.1),
          numericInput("tauT2", "\u03c4 T2", value = 2.8, step = 0.1),
          numericInput("theta1", "\u03b8 T1 (active lag)", value = 0.3, step = 0.05),
          numericInput("theta2", "\u03b8 T2 (active lag)", value = 0.3, step = 0.05),
          numericInput("seed_y", "Seed (realization of Y)", value = 10039, step = 1),
          actionButton("gen", "Generate design and data", class = "btn-primary", icon = icon("play"),
                       style = "width:100%; margin-top:10px;"),
          br(), br(),
          actionButton("findDemo", "Find illustrative sample",
                       class = "btn-primary", icon = icon("magnifying-glass"),
                       style = "width:100%;"),
          helpText("Searches, among error resamples, for one in which Moran (usual) gives p < 0.05 and Moran (SAV) gives p \u2265 0.05 \u2014 that is, a case in which incorporating the lag corrects the spatial dependence of the residuals."),
          textOutput("demoStatusText")
        ),
        
        mainPanel(
          width = 9,
          tabsetPanel(
            type = "tabs",
            
            tabPanel("Design and data",
                     br(),
                     uiOutput("design_summary"),
                     plotOutput("field_map", height = "420px"),
                     br(),
                     plotOutput("response_map", height = "420px"),
                     br(),
                     DTOutput("data_table")
            ),
            
            tabPanel("Matrices and verification",
                     br(),
                     h4("1. Original design matrix X = [J, X\u03c4] (rank-deficient)"),
                     verbatimTextOutput("mat_X_info"),
                     br(),
                     h4("2. Singular values of X and transformation matrix V_r"),
                     p(style = "color:#7f8c8d;",
                       "X = UDV', with rank \u03c1(X) = 3 < 4 (since X\u03c4 \u00b7 1\u2083 = J). V_r holds the right singular vectors associated with the non-zero singular values."),
                     verbatimTextOutput("mat_svd_info"),
                     br(),
                     h4("3. Full-rank design matrix Z = X V_r"),
                     verbatimTextOutput("mat_Z_info"),
                     br(),
                     h4("4. Spatial weights matrix W"),
                     p(style = "color:#7f8c8d;",
                       "Note: Condition S3 (W non-singular) is sufficient but not necessary; what the model actually requires is that the matrix A below has rank 2, which is verified directly in item 6, even if W itself is rank-deficient (common for contiguity matrices on regular grids)."),
                     verbatimTextOutput("mat_W_info"),
                     br(),
                     h4("5. Active-treatment lag matrix X_{W,act} = W X_act"),
                     verbatimTextOutput("mat_XWact_info"),
                     br(),
                     h4("6. Coefficients that verify Condition S4: A = (I - M_Z) X_{W,act} and A'A"),
                     p(style = "color:#7f8c8d;",
                       "If \u03c1(A) = 2 (equivalently, det(A'A) \u2260 0), Condition S4 holds numerically for this design and this W: the active lag is estimable with its full 2 degrees of freedom."),
                     verbatimTextOutput("mat_A_info")
            ),
            
            tabPanel("Usual ANOVA (no lag)",
                     br(),
                     p(style = "color:#7f8c8d;",
                       "Model that ignores the active spatial lag: Y = Z\u03b4 + \u03b5. The residual mean square of this model is a biased estimator of \u03c3\u00b2 when the true lag is non-zero (Section 2.4.1)."),
                     DTOutput("table_usual"),
                     br(),
                     div(class = "well",
                         h5(style = "margin-top:0;", "Spatial dependence of the residuals (Moran's I)"),
                         verbatimTextOutput("moran_usual")
                     )
            ),
            
            tabPanel("Spatial ANOVA (SAV)",
                     br(),
                     p(style = "color:#7f8c8d;",
                       "Full decomposition via the Halmos theorem: Y = Z\u03b4 + X_{W,act}\u03b8_act + \u03b5."),
                     DTOutput("table_spatial"),
                     br(),
                     div(class = "well",
                         h5(style = "margin-top:0;", "Spatial dependence of the residuals (Moran's I)"),
                         verbatimTextOutput("moran_spatial")
                     )
            ),
            
            tabPanel("Simulation (Monte Carlo grid)",
                     br(),
                     div(style = "background:#eef6fb; padding:12px; border-radius:8px; border-left: 4px solid #2980b9;",
                         p(style="margin:0; color:#2c3e50; font-size:13px;",
                           strong("Note: "), "this tab uses its own fixed design (Rook, treatment seed = 1, \u03c4 = 0, \u03c3\u00b2 = 1), independent of the values in the left panel, to reproduce exactly Table 2 / Figure 2 of the paper. The sidebar fields (Queen, \u03c4\u22600, \u03c3=0.31) apply only to the \"Design and data\", \"Usual ANOVA\" and \"Spatial ANOVA (SAV)\" tabs, which reproduce Tables 3/4.")
                     ),
                     br(),
                     div(style = "background:#f8f9fa; padding:16px; border-radius:10px;",
                         fluidRow(
                           column(3, textInput("cs_grid", "Grid of c (comma-separated)",
                                               value = "0, 1, 2, 3, 4, 6, 8")),
                           column(3, numericInput("B_mc", "Monte Carlo replications (B)",
                                                  value = 20000, min = 500, max = 100000, step = 500)),
                           column(3, numericInput("seed_mc", "Seed (simulation)", value = 123, step = 1)),
                           column(3, br(), actionButton("run_mc", "Run simulation",
                                                        class = "btn-primary", icon = icon("play"),
                                                        style = "width:100%;"))
                         )
                     ),
                     br(),
                     DTOutput("table_grid"),
                     br(),
                     h4("Power curve"),
                     p(style = "color:#7f8c8d;",
                       "Empirical rejection rate as a function of the true active-lag magnitude (c), for the correctly specified Treatment test (via SAV), the naive Treatment test (which ignores the lag), and the Active-lag test. The dashed line marks the nominal 5% level."),
                     plotOutput("power_plot", height = "420px"),
                     br(),
                     downloadButton("download_plot", "Download plot (PNG)"),
                     br(), br(),
                     h4("LaTeX code (for the appendix)"),
                     downloadButton("download_tex", "Download .tex table"),
                     br(), br(),
                     verbatimTextOutput("latex_preview")
            )
          )
        )
      )
  )
)

# ══════════════════════════════════════════════════════════════
# Server
# ══════════════════════════════════════════════════════════════

server <- function(input, output, session) {
  
  rv <- reactiveValues(design = NULL, Y = NULL, grid = NULL, grid_B = NULL, grid_seed = NULL)
  
  observeEvent(input$gen, {
    d <- build_design(neighbor_type = input$neighbor_type, seed_treat = input$seed_treat)
    tau <- c(input$tauC, input$tauT1, input$tauT2)
    theta_act <- c(input$theta1, input$theta2)
    Y <- generate_Y(d, tau, theta_act, input$sigma, input$seed_y)
    rv$design <- d
    rv$Y <- Y
  }, ignoreNULL = FALSE)
  
  # observeEvent above already runs once on load, since ignoreNULL = FALSE
  
  output$design_summary <- renderUI({
    req(rv$design)
    d <- rv$design
    div(class = "well",
        p(HTML(paste0(
          "<b>Grid</b> = ", d$grid_size, "\u00d7", d$grid_size, " (n = ", d$n, ") &nbsp; | &nbsp; ",
          "<b>t</b> = ", d$t_levels, " (Control, T1, T2) &nbsp; | &nbsp; ",
          "<b>r</b> = ", d$r, " replicates/level &nbsp; | &nbsp; ",
          "<b>Neighbourhood</b> = ", ifelse(d$neighbor_type == "rook", "Rook (4)", "Queen (8)"), " &nbsp; | &nbsp; ",
          "df Treatment = ", d$df_treat, " &nbsp; | &nbsp; ",
          "df Spatial lag = ", d$df_lag, " &nbsp; | &nbsp; ",
          "df Error (SAV) = ", d$df_err, " &nbsp; | &nbsp; ",
          "df Error (usual) = ", d$df_WL
        )))
    )
  })
  
  output$field_map <- renderPlot({
    req(rv$design)
    d <- rv$design
    df <- data.frame(row = d$cell_row, col = d$cell_col, level = d$levels_vec)
    ggplot(df, aes(x = col, y = row, fill = level)) +
      geom_tile(color = "white", linewidth = 0.8) +
      geom_text(aes(label = level), color = "white", fontface = "bold", size = 3) +
      scale_fill_manual(values = c("C" = "#7f8c8d", "T1" = "#2980b9", "T2" = "#e74c3c")) +
      scale_y_reverse(breaks = 1:9) +
      scale_x_continuous(breaks = 1:9) +
      coord_fixed() +
      labs(title = "Treatment allocation (completely randomized design, 81 cells, 27/level)",
           subtitle = "Scattered on purpose: the randomization of a CRD should show no spatial pattern.",
           x = "Column", y = "Row", fill = "Level") +
      theme_minimal(base_size = 13)
  })
  
  output$response_map <- renderPlot({
    req(rv$design, rv$Y)
    d <- rv$design
    df <- data.frame(row = d$cell_row, col = d$cell_col, level = d$levels_vec, Y = rv$Y)
    ggplot(df, aes(x = col, y = row, fill = Y)) +
      geom_tile(color = "white", linewidth = 0.8) +
      geom_text(aes(label = level), color = "white", fontface = "bold", size = 2.8) +
      scale_fill_viridis_c(option = "magma") +
      scale_y_reverse(breaks = 1:9) +
      scale_x_continuous(breaks = 1:9) +
      coord_fixed() +
      labs(title = "Map of the response Y (this is where spatial clustering should be visible)",
           subtitle = "Colour shows the value of Y; treatment labels are kept for reference only.",
           x = "Column", y = "Row", fill = "Y") +
      theme_minimal(base_size = 13)
  })
  
  output$data_table <- renderDT({
    req(rv$design, rv$Y)
    d <- rv$design
    df <- data.frame(
      Cell = seq_len(d$n),
      Row = d$cell_row,
      Column = d$cell_col,
      Level = d$levels_vec,
      Y_response = round(rv$Y, 4)
    )
    datatable(df, options = list(pageLength = 10), rownames = FALSE)
  })
  
  output$mat_X_info <- renderPrint({
    req(rv$design)
    d <- rv$design
    cat("dim(X) =", dim(d$X)[1], "x", dim(d$X)[2], "\n")
    cat("rank(X) =", d$r_rank, "(deficient: 4 columns, rank 3)\n\n")
    cat("First 8 rows of X = [J, X_C, X_1, X_2]:\n")
    print(head(d$X, 8))
  })
  
  output$mat_svd_info <- renderPrint({
    req(rv$design)
    d <- rv$design
    cat("Singular values of X (d):\n")
    print(round(d$sv$d, 6))
    cat("\nTransformation matrix V_r (4 x", d$r_rank, "):\n")
    print(round(d$Vr, 6))
  })
  
  output$mat_Z_info <- renderPrint({
    req(rv$design)
    d <- rv$design
    cat("dim(Z) =", dim(d$Z)[1], "x", dim(d$Z)[2], "\n")
    cat("rank(Z) =", qr(d$Z)$rank, "(full)\n\n")
    cat("First 8 rows of Z:\n")
    print(round(head(d$Z, 8), 6))
  })
  
  output$mat_W_info <- renderPrint({
    req(rv$design)
    d <- rv$design
    cat("dim(W) =", dim(d$W)[1], "x", dim(d$W)[2], "\n")
    cat("rank(W) =", d$W_rank, "of", d$n, "\n")
    cat("Neighbourhood type:", d$neighbor_type, "\n\n")
    cat("Submatrix W[1:8, 1:8] (upper-left corner):\n")
    print(round(d$W[1:8, 1:8], 3))
  })
  
  output$mat_XWact_info <- renderPrint({
    req(rv$design)
    d <- rv$design
    cat("dim(X_{W,act}) =", dim(d$XWact)[1], "x", dim(d$XWact)[2], "\n\n")
    cat("First 10 rows of X_{W,act} = [W X_1, W X_2]:\n")
    print(round(head(d$XWact, 10), 4))
  })
  
  output$mat_A_info <- renderPrint({
    req(rv$design)
    d <- rv$design
    cat("rank(A) =", d$A_rank, "(2 is required for S4 to hold)\n\n")
    cat("A'A (2 x 2) -- the coefficients that verify Condition S4:\n")
    print(round(d$AtA, 6))
    cat("\ndet(A'A) =", d$det_AtA, "\n")
    if (d$A_rank == 2 && abs(d$det_AtA) > 1e-8) {
      cat("\n=> Condition S4 (weak, rank-based) verified: theta_1 and theta_2 are\n")
      cat("   individually and jointly estimable under this design and this W.\n")
    } else {
      cat("\n=> WARNING: Condition S4 does NOT hold for this W / this treatment\n")
      cat("   layout. Try the other neighbourhood type.\n")
    }
    cat("\n---\n")
    cat("Orthogonality leakage ||M_Z X_{W,act}|| / ||X_{W,act}|| =", round(d$leak_Z, 4), "\n")
    cat("(0 = X_{W,act} fully orthogonal to Z, the ideal condition for\n")
    cat(" E(MS_treatment) to be free of theta_act; larger values indicate that the\n")
    cat(" Treatment mean square may be influenced by the active lag.)\n")
  })
  
  output$moran_usual <- renderPrint({
    req(rv$design, rv$Y)
    d <- rv$design
    e_usual <- as.numeric((d$I_n - d$MZ) %*% rv$Y)  # residuals WITHOUT the lag term
    mt <- moran_test(e_usual, d$W, nsim = 999)
    cat("Moran's I =", round(mt$I, 4), "\n")
    cat("p-value (permutation, 999 resamples) =", round(mt$p.value, 4), "\n")
    if (mt$p.value < 0.05) {
      cat("=> SIGNIFICANT spatial autocorrelation in the residuals of the usual model.\n")
    } else {
      cat("=> No evidence of spatial autocorrelation in the residuals of the usual model.\n")
    }
  })
  
  output$moran_spatial <- renderPrint({
    req(rv$design, rv$Y)
    d <- rv$design
    e_sav <- as.numeric((d$I_n - d$ML) %*% rv$Y)  # residuals WITH the active-treatment lag term
    mt <- moran_test(e_sav, d$W, nsim = 999)
    cat("Moran's I =", round(mt$I, 4), "\n")
    cat("p-value (permutation, 999 resamples) =", round(mt$p.value, 4), "\n")
    if (mt$p.value < 0.05) {
      cat("=> SIGNIFICANT spatial autocorrelation persists even after incorporating the lag.\n")
    } else {
      cat("=> After incorporating the active lag, the spatial autocorrelation of the residuals is no longer significant.\n")
    }
  })
  
  observeEvent(input$findDemo, {
    req(rv$design)
    d <- rv$design
    tau <- c(input$tauC, input$tauT1, input$tauT2)
    theta_act <- c(input$theta1, input$theta2)
    
    max_tries <- 300
    found <- FALSE
    result_seed <- NA
    result_attempt <- NA
    result_Y <- NULL
    
    withProgress(message = "Searching for an illustrative sample...", value = 0, {
      for (i in seq_len(max_tries)) {
        incProgress(1 / max_tries, detail = paste("Attempt", i, "of", max_tries))
        s <- 10000 + i
        Y_try <- generate_Y(d, tau, theta_act, input$sigma, s)
        e_usual <- as.numeric((d$I_n - d$MZ) %*% Y_try)
        e_sav <- as.numeric((d$I_n - d$ML) %*% Y_try)
        p_usual <- moran_test(e_usual, d$W, nsim = 199)$p.value
        p_sav <- moran_test(e_sav, d$W, nsim = 199)$p.value
        if (p_usual < 0.05 && p_sav >= 0.05) {
          found <- TRUE
          result_seed <- s
          result_attempt <- i
          result_Y <- Y_try
          break
        }
      }
    })
    
    if (found) {
      rv$Y <- result_Y
      updateNumericInput(session, "seed_y", value = result_seed)
      output$demoStatusText <- renderText(
        paste0("\u2713 Found at attempt ", result_attempt,
               " (seed ", result_seed, "). Check the ANOVA tabs.")
      )
    } else {
      output$demoStatusText <- renderText(
        paste0("No sample meeting both conditions was found in ", max_tries,
               " attempts. Try increasing \u03b8 (lag) or \u03c4 (treatment effects).")
      )
    }
  })
  
  output$table_usual <- renderDT({
    req(rv$design, rv$Y)
    tab <- anova_usual(rv$design, rv$Y)
    tab[, c("SS", "MS", "F", "p_value")] <- round(tab[, c("SS", "MS", "F", "p_value")], 4)
    datatable(tab, options = list(dom = "t", pageLength = 10), rownames = FALSE)
  })
  
  output$table_spatial <- renderDT({
    req(rv$design, rv$Y)
    tab <- anova_spatial(rv$design, rv$Y)
    tab[, c("SS", "MS", "F", "p_value")] <- round(tab[, c("SS", "MS", "F", "p_value")], 4)
    datatable(tab, options = list(dom = "t", pageLength = 10), rownames = FALSE)
  })
  
  observeEvent(input$run_mc, {
    cs <- as.numeric(trimws(strsplit(input$cs_grid, ",")[[1]]))
    # Self-contained design for the Monte Carlo grid: always rook contiguity,
    # treatment-randomization seed = 1, tau = 0 -- matching Sections 3.1.1/3.1.2
    # of the paper (Table 2 / Figure 2), regardless of whatever the sidebar
    # is currently set to for the separate Tables 3/4 illustrative example.
    d_mc <- build_design(neighbor_type = "rook", seed_treat = 1)
    tau_mc <- c(0, 0, 0)
    withProgress(message = "Running Monte Carlo simulation...", value = 0.3, {
      grid <- run_grid(d_mc, tau_mc, cs, B = input$B_mc, seed = input$seed_mc)
      incProgress(1)
    })
    rv$grid <- grid
    rv$grid_B <- input$B_mc
    rv$grid_seed <- input$seed_mc
  })
  
  output$table_grid <- renderDT({
    req(rv$grid)
    tab <- rv$grid
    num_cols <- setdiff(names(tab), "c")
    tab[, num_cols] <- round(tab[, num_cols], 4)
    datatable(tab, options = list(dom = "t", pageLength = 15), rownames = FALSE)
  })
  
  # ── Power curve plot (reactive object, reused by both renderPlot and the
  #    PNG download handler so the on-screen plot and the downloaded file
  #    are guaranteed to be identical) ──────────────────────────────────
  power_plot_obj <- reactive({
    req(rv$grid)
    df_long <- data.frame(
      c = rep(rv$grid$c, 3),
      rejection = c(rv$grid$rej_treat_correct, rv$grid$rej_treat_naive, rv$grid$rej_lag),
      Test = rep(c("Treatment (correct, SAV)", "Treatment (naive, unlagged)", "Active lag"),
                 each = nrow(rv$grid))
    )
    ggplot(df_long, aes(x = c, y = rejection, color = Test, shape = Test)) +
      geom_line(linewidth = 0.9) +
      geom_point(size = 2.8) +
      geom_hline(yintercept = 0.05, linetype = "dashed", color = "black") +
      scale_color_manual(values = c("Treatment (correct, SAV)" = "#2980b9",
                                    "Treatment (naive, unlagged)" = "#e74c3c",
                                    "Active lag" = "#27ae60")) +
      labs(title = "Empirical rejection rate vs. active-treatment lag magnitude (c)",
           subtitle = "Dashed line: nominal 5% level",
           x = "Lag magnitude (c)", y = "Empirical rejection rate") +
      theme_minimal(base_size = 13) +
      theme(legend.position = "bottom", legend.title = element_blank())
  })
  
  output$power_plot <- renderPlot({ power_plot_obj() })
  
  output$download_plot <- downloadHandler(
    filename = function() paste0("sav_power_curve_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png"),
    content = function(file) {
      ggsave(file, plot = power_plot_obj(), width = 8, height = 5.5, dpi = 300)
    }
  )
  
  output$latex_preview <- renderText({
    req(rv$grid)
    grid_to_latex(rv$grid, B = rv$grid_B, seed_mc = rv$grid_seed)
  })
  
  output$download_tex <- downloadHandler(
    filename = function() paste0("sav_montecarlo_grid_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".tex"),
    content = function(file) {
      req(rv$grid)
      writeLines(grid_to_latex(rv$grid, B = rv$grid_B, seed_mc = rv$grid_seed), file)
    }
  )
}

shinyApp(ui = ui, server = server)