library(shiny)
library(ggplot2)
library(dplyr)

# ── CONFUSION MATRIX (potato crop MLP — fixed) ───────────────────────────────
TP_fix <- 9; FP_fix <- 10; FN_fix <- 18; TN_fix <- 175
N_fix  <- TP_fix + FP_fix + FN_fix + TN_fix   # 212

# ── METRIC FUNCTIONS ──────────────────────────────────────────────────────────

CA_fn <- function(Cse, Csp, P) P*Cse + (1-P)*Csp

BG_fn <- function(Cse, Csp, P) {
  A <- P*Cse + (1-P)*(1-Csp)
  B <- P*(1-Cse) + (1-P)*Csp
  if (is.na(A)||is.na(B)||A<1e-9||B<1e-9) return(NA_real_)
  sqrt((P*Cse/A) * ((1-P)*Csp/B))
}

CJ_fn <- function(Cse, Csp) Cse + Csp - 1

BJ_fn <- function(Cse, Csp, P) {
  A <- P*Cse + (1-P)*(1-Csp)
  B <- P*(1-Cse) + (1-P)*Csp
  if (is.na(A)||is.na(B)||A<1e-9||B<1e-9) return(NA_real_)
  P*Cse/A + (1-P)*Csp/B - 1
}

Pstars_fn <- function(Cse, Csp) {
  f_P <- function(P) {
    CA <- P*Cse + (1-P)*Csp
    A  <- P*Cse + (1-P)*(1-Csp)
    B  <- P*(1-Cse) + (1-P)*Csp
    CA^2 * A * B - P*(1-P)*Cse*Csp
  }
  grid  <- seq(0.001, 0.999, by = 0.0002)
  fvals <- vapply(grid, f_P, numeric(1))
  chg   <- which(diff(sign(fvals)) != 0)
  roots <- vapply(chg, function(i) {
    tryCatch(
      uniroot(f_P, lower=grid[i], upper=grid[i+1], tol=1e-12)$root,
      error=function(e) NA_real_)
  }, numeric(1))
  roots <- sort(roots[!is.na(roots) & roots > 1e-6 & roots < 1-1e-6])
  list(p1=if(length(roots)>=1) roots[1] else NA_real_,
       p2=if(length(roots)>=2) roots[2] else NA_real_)
}

# ── COMPARISON METRICS (from fixed confusion matrix) ─────────────────────────

metrics_table <- function(Cse, Csp, P) {
  TP <- TP_fix; FP <- FP_fix; FN <- FN_fix; TN <- TN_fix; N <- N_fix
  
  # classical accuracy
  CA  <- (TP + TN) / N
  
  # F1
  prec <- TP / (TP + FP)
  rec  <- TP / (TP + FN)
  F1   <- if ((prec + rec) > 0) 2 * prec * rec / (prec + rec) else 0
  
  # MCC
  denom_mcc <- sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
  MCC <- if (denom_mcc > 0) (TP*TN - FP*FN) / denom_mcc else 0
  
  # G-mean (intrinsic space: sqrt(Se * Sp))
  Gmean <- sqrt(Cse * Csp)
  
  # BG at observed prevalence
  BG <- BG_fn(Cse, Csp, P)
  
  # PPV and NPV at observed prevalence (Bayesian)
  A   <- P*Cse + (1-P)*(1-Csp)
  B   <- P*(1-Cse) + (1-P)*Csp
  PPV <- if (!is.na(A) && A > 1e-9) P*Cse/A else NA
  NPV <- if (!is.na(B) && B > 1e-9) (1-P)*Csp/B else NA
  
  list(CA=CA, F1=F1, MCC=MCC, Gmean=Gmean, BG=BG, PPV=PPV, NPV=NPV,
       prec=prec, rec=rec)
}

# ── COLOURS ───────────────────────────────────────────────────────────────────
CC   <- "#1b4f72"
CB   <- "#922b21"
CP1  <- "#1e8449"
CP2  <- "#7d6608"
CO   <- "#117a65"
COJ  <- "#1a5276"
CBJ  <- "#6e2f1a"
CPU  <- "#6c3483"
CGR  <- "#aaaaaa"
CTX  <- "#444444"
CSB  <- "#555555"
CCP  <- "#777777"
CAX  <- "#333333"

# ── THEME ─────────────────────────────────────────────────────────────────────
tp <- function() {
  theme_minimal(base_size = 12) +
    theme(
      plot.title    = element_text(face="bold", hjust=.5, size=13, colour="#1a1a2e"),
      plot.subtitle = element_text(hjust=.5, size=10, colour=CSB),
      plot.caption  = element_text(size=8, colour=CCP, hjust=0),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour="#ebebeb", linewidth=.35),
      panel.border     = element_rect(colour="#aaaaaa", fill=NA, linewidth=.7),
      axis.title       = element_text(face="bold", size=11, colour="#1a1a2e"),
      axis.text        = element_text(size=9, colour=CAX),
      legend.position  = "bottom",
      legend.title     = element_blank(),
      legend.text      = element_text(size=9),
      legend.key.width = unit(1.4, "cm"),
      plot.margin      = margin(10, 14, 8, 10)
    )
}

# ── PLOT 1 ────────────────────────────────────────────────────────────────────
make_plot1 <- function(Cse, Csp, Pobs) {
  ps    <- Pstars_fn(Cse, Csp)
  prevs <- seq(0.005, 0.995, by=0.005)
  p1    <- if(!is.na(ps$p1)) ps$p1 else 0.35
  p2    <- if(!is.na(ps$p2)) ps$p2 else 0.75
  
  df <- data.frame(
    P   = c(prevs, prevs),
    val = c(sapply(prevs, function(p) CA_fn(Cse,Csp,p)),
            sapply(prevs, function(p) BG_fn(Cse,Csp,p))),
    met = rep(c("Ca","Bg"), each=length(prevs))
  )
  
  if (Cse >= Csp) { fL <- "#fde8e8"; fM <- "#e8f8ee"; fR <- "#fde8e8"
  } else          { fL <- "#e8f8ee"; fM <- "#fde8e8"; fR <- "#e8f8ee" }
  
  lbl_L <- if(Cse>=Csp) "Bg > Ca" else "Ca > Bg"
  lbl_M <- if(Cse>=Csp) "Ca > Bg" else "Bg > Ca"
  
  g <- ggplot(df, aes(x=P, y=val, colour=met)) +
    annotate("rect", xmin=0,  xmax=p1, ymin=0, ymax=1, fill=fL, alpha=.45) +
    annotate("rect", xmin=p1, xmax=p2, ymin=0, ymax=1, fill=fM, alpha=.45) +
    annotate("rect", xmin=p2, xmax=1,  ymin=0, ymax=1, fill=fL, alpha=.45) +
    annotate("text", x=(0+p1)/2,  y=.95, label=lbl_L,
             size=2.8, colour=CTX, hjust=.5) +
    annotate("text", x=(p1+p2)/2, y=.95, label=lbl_M,
             size=2.8, colour=CTX, hjust=.5) +
    annotate("text", x=(p2+1)/2,  y=.95, label=lbl_L,
             size=2.8, colour=CTX, hjust=.5) +
    geom_vline(xintercept=.5, linetype="dotted", colour=CGR, linewidth=.5) +
    annotate("text", x=.51, y=.03, label="P=0.5",
             size=2.6, colour=CGR, hjust=0) +
    geom_line(linewidth=1.7, na.rm=TRUE) +
    scale_colour_manual(
      values = c(Ca=CC, Bg=CB),
      labels = c(
        Ca = "Ca = P*Cse + (1-P)*Csp   (classical accuracy)",
        Bg = "Bg = sqrt(PPV*NPV)   (Bayesian geometric accuracy)"
      )
    ) +
    scale_x_continuous(breaks=seq(0,1,.1), limits=c(0,1), expand=c(.005,.005)) +
    scale_y_continuous(breaks=seq(0,1,.1), limits=c(0,1)) +
    labs(
      title    = "Ca vs Bg as a function of Prevalence",
      subtitle = paste0("Cse=",Cse,"   Csp=",Csp,
                        "   CJ=",round(CJ_fn(Cse,Csp),3)),
      x = "Prevalence  P",
      y = "Metric value",
      caption = "Green region: Bg >= Ca  |  Red region: Ca > Bg (overestimation)"
    ) + tp()
  
  if (!is.na(ps$p1)) {
    y1ca <- CA_fn(Cse,Csp,ps$p1)
    y1bg <- BG_fn(Cse,Csp,ps$p1)
    nx   <- ps$p1 + ifelse(ps$p1 < .55,  .02, -.02)
    hj   <- ifelse(ps$p1 < .55, 0, 1)
    g <- g +
      geom_vline(xintercept=ps$p1, linetype="dashed", colour=CP1, linewidth=.9) +
      annotate("point", x=ps$p1, y=y1ca, colour=CP1, size=3.5, shape=18) +
      annotate("point", x=ps$p1, y=y1bg, colour=CP1, size=3.5, shape=18) +
      annotate("text", x=nx, y=min(y1ca+.08, .90),
               label=paste0("P*A-1 = ",round(ps$p1,3)),
               size=3.0, hjust=hj, colour=CP1, fontface="bold") +
      annotate("text", x=nx, y=min(y1bg+.07, .82),
               label=paste0("Popt = ",round(ps$p1,3)),
               size=2.7, hjust=hj, colour=CPU)
  }
  
  if (!is.na(ps$p2)) {
    y2ca <- CA_fn(Cse,Csp,ps$p2)
    y2bg <- BG_fn(Cse,Csp,ps$p2)
    nx   <- ps$p2 + ifelse(ps$p2 < .75,  .02, -.02)
    hj   <- ifelse(ps$p2 < .75, 0, 1)
    g <- g +
      geom_vline(xintercept=ps$p2, linetype="dotdash", colour=CP2, linewidth=.9) +
      annotate("point", x=ps$p2, y=y2ca, colour=CP2, size=3.5, shape=18) +
      annotate("point", x=ps$p2, y=y2bg, colour=CP2, size=3.5, shape=18) +
      annotate("text", x=nx, y=max(y2ca-.09, .07),
               label=paste0("P*A-2 = ",round(ps$p2,3)),
               size=3.0, hjust=hj, colour=CP2, fontface="bold")
  }
  
  if (!is.na(Pobs) && Pobs > 0 && Pobs < 1) {
    cao <- CA_fn(Cse,Csp,Pobs)
    bgo <- BG_fn(Cse,Csp,Pobs)
    nx  <- Pobs + ifelse(Pobs < .65, .02, -.02)
    hj  <- ifelse(Pobs < .65, 0, 1)
    g <- g +
      geom_vline(xintercept=Pobs, linetype="longdash", colour=CO, linewidth=.9) +
      annotate("point", x=Pobs, y=cao, colour=CO, size=3.5, shape=16) +
      annotate("point", x=Pobs, y=bgo, colour=CO, size=3.5, shape=16) +
      annotate("text", x=nx, y=(cao+bgo)/2+.03,
               label=paste0("Pobs = ",round(Pobs,3)),
               size=3.0, hjust=hj, colour=CO, fontface="bold") +
      annotate("text", x=nx, y=cao+.05,
               label=paste0("Ca = ",round(cao,3)),
               size=2.7, hjust=hj, colour=CC) +
      annotate("text", x=nx, y=bgo-.06,
               label=paste0("Bg = ",round(bgo,3)),
               size=2.7, hjust=hj, colour=CB)
  }
  g
}

# ── DOWNLOAD ──────────────────────────────────────────────────────────────────
save_plot <- function(p, file, w=10, h=6, dpi=150) {
  if (grepl("\\.pdf$", file)) {
    pdf(file, width=w, height=h); print(p); dev.off()
  } else {
    png(file, width=round(w*dpi), height=round(h*dpi), res=dpi)
    print(p); dev.off()
  }
}

# ── UI ────────────────────────────────────────────────────────────────────────
ui <- fluidPage(
  tags$head(tags$style(HTML("
    body  { background:#f4f6f8; font-family:'Georgia',serif; font-size:13px; }
    .well { background:#ffffff; border:1px solid #d5d8dc; border-radius:7px; padding:12px; }
    h3    { color:#1a2535; font-weight:bold; font-size:16px; margin:0 0 3px 0; }
    h5    { color:#1a2535; font-weight:bold; font-size:12px; margin:8px 0 2px 0; }
    .mbox { background:#eaf0fb; border-left:3px solid #2471a3; padding:5px 9px;
            margin:2px 0; border-radius:3px; font-size:11px; font-family:monospace; }
    .pbox { background:#fef9e7; border-left:3px solid #d4ac0d; padding:5px 9px;
            margin:2px 0; border-radius:3px; font-size:11px; font-family:monospace; }
    .gbox { background:#eafaf1; border-left:3px solid #1e8449; padding:5px 9px;
            margin:2px 0; border-radius:3px; font-size:11px; font-family:monospace; }
    .rbox { background:#fdedec; border-left:3px solid #922b21; padding:5px 9px;
            margin:2px 0; border-radius:3px; font-size:11px; font-family:monospace; }
    .nav-tabs>li>a { font-size:12px; font-weight:bold; color:#1a2535; padding:6px 12px; }
    .dl-row { display:flex; gap:6px; padding:6px 0 4px 0; }
    .btn-xs { font-size:11px !important; padding:3px 10px !important; }

    /* ── comparison table ── */
    .tbl-wrap  { overflow-x:auto; padding:6px 0; }
    table.ctbl { border-collapse:collapse; width:100%; font-size:12px; }
    table.ctbl th {
      background:#1b4f72; color:#ffffff; padding:8px 10px;
      text-align:center; font-weight:bold; border:1px solid #aaaaaa; }
    table.ctbl td {
      padding:7px 10px; border:1px solid #cccccc;
      text-align:center; }
    table.ctbl tr:nth-child(even) td { background:#f2f6fc; }
    table.ctbl tr:last-child   td { background:#eafaf1; font-weight:bold; }
    .yes  { color:#1e8449; font-weight:bold; }
    .no   { color:#922b21; }
    .part { color:#7d6608; font-style:italic; }
    .note-box {
      background:#fef9e7; border-left:4px solid #d4ac0d;
      padding:10px 14px; margin:10px 0; border-radius:4px;
      font-size:11.5px; line-height:1.6; }
    .ref-cm {
      background:#eaf0fb; border-left:4px solid #2471a3;
      padding:8px 12px; margin:6px 0 10px 0; border-radius:4px;
      font-size:11px; font-family:monospace; line-height:1.8; }
  "))),
  
  div(style="padding:10px 0 4px 0;",
      h3("Bg = sqrt(PPV*NPV) — Bayesian Geometric Accuracy Explorer"),
      p(style="color:#666666;font-size:11px;margin:0;",
        "Darghan et al. (2026) — Potato crop MLP classifier")
  ),
  
  sidebarLayout(
    sidebarPanel(width=3,
                 wellPanel(
                   h5("Classifier parameters"),
                   sliderInput("Cse",  HTML("C<sub>se</sub> (sensitivity)"),
                               0.10, 0.99, 0.333, step=0.001),
                   sliderInput("Csp",  HTML("C<sub>sp</sub> (specificity)"),
                               0.10, 0.99, 0.946, step=0.001),
                   sliderInput("Pobs", HTML("P<sub>obs</sub> (observed prevalence)"),
                               0.01, 0.99, 0.127, step=0.001)
                 ),
                 wellPanel(h5("Intrinsic metrics"),  uiOutput("sb_met")),
                 wellPanel(h5("Equilibrium points"), uiOutput("sb_ps")),
                 wellPanel(h5(HTML("At P<sub>obs</sub>")), uiOutput("sb_obs"))
    ),
    
    mainPanel(width=9,
              tabsetPanel(type="tabs",
                          
                          # ── TAB 1: Plot ──────────────────────────────────────────
                          tabPanel(HTML("Plot — C<sub>A</sub> vs B<sub>G</sub>"),
                                   br(),
                                   div(class="dl-row",
                                       downloadButton("dl1p","PNG",class="btn-xs btn-default"),
                                       downloadButton("dl1f","PDF",class="btn-xs btn-default")),
                                   plotOutput("p1", height="450px")
                          ),
                          
                          # ── TAB 2: Metric comparison table ───────────────────────
                          tabPanel("Metric Comparison",
                                   br(),
                                   div(class="ref-cm",
                                       strong("Fixed confusion matrix — Potato crop MLP (n = 212)"), br(),
                                       "TP = 9  |  FP = 10  |  FN = 18  |  TN = 175", br(),
                                       HTML("C<sub>se</sub> = 9/27 &asymp; 0.333  &nbsp;|&nbsp;
                  C<sub>sp</sub> = 175/185 &asymp; 0.946  &nbsp;|&nbsp;
                  P<sub>obs</sub> = 27/212 &asymp; 0.127")
                                   ),
                                   
                                   h5(style="color:#1a2535;margin:14px 0 6px 0;",
                                      "Numerical values at current slider settings"),
                                   div(class="tbl-wrap", tableOutput("num_table")),
                                   
                                   br(),
                                   h5(style="color:#1a2535;margin:10px 0 6px 0;",
                                      "Formal properties — structural comparison"),
                                   div(class="tbl-wrap", tableOutput("prop_table")),
                                   
                                   div(class="note-box",
                                       HTML("<b>Reading the table.</b>
            F1, MCC and G-mean are <i>scalar summaries</i> computed at a
            fixed operating point: they do not express prevalence as a formal
            argument and cannot show how performance changes as the deployment
            distribution shifts.
            G-mean is the closest structural analog to B<sub>G</sub> — both are
            geometric means — but it operates in the intrinsic space of
            (C<sub>se</sub>, C<sub>sp</sub>), answering whether the classifier
            <i>discriminates</i>, rather than in the predictive space of
            (PPV, NPV), answering whether <i>predictions can be trusted</i>.
            B<sub>G</sub> complements rather than replaces these metrics: it adds
            a prevalence-sensitive, prediction-space perspective that scalar
            summaries structurally cannot provide, and uniquely generates an
            algebraically grounded operating interval identifying which
            prevalence regimes render classical accuracy misleading.")
                                   )
                          )
              )
    )
  )
)

# ── SERVER ────────────────────────────────────────────────────────────────────
server <- function(input, output, session) {
  
  v <- reactive({
    req(input$Cse, input$Csp, input$Pobs)
    ps <- Pstars_fn(input$Cse, input$Csp)
    list(Cse=input$Cse, Csp=input$Csp, Pobs=input$Pobs,
         ps=ps, CJ=CJ_fn(input$Cse,input$Csp))
  })
  
  # ── sidebar panels ──────────────────────────────────────────────
  output$sb_met <- renderUI({
    x <- v()
    tagList(
      div(class="mbox", sprintf("CJ = %.4f", x$CJ)),
      div(class="mbox", sprintf("u = Cse(1-Cse) = %.4f", x$Cse*(1-x$Cse))),
      div(class="mbox", sprintf("v = Csp(1-Csp) = %.4f", x$Csp*(1-x$Csp))),
      div(class="mbox", sprintf("sqrt(uv) = %.4f",
                                sqrt(x$Cse*(1-x$Cse)*x$Csp*(1-x$Csp))))
    )
  })
  
  output$sb_ps <- renderUI({
    x <- v()
    w <- if(!is.na(x$ps$p1)&&!is.na(x$ps$p2)) x$ps$p2-x$ps$p1 else NA
    tagList(
      div(class="pbox", sprintf("P*A-1 = %.4f", ifelse(is.na(x$ps$p1),NA,x$ps$p1))),
      div(class="pbox", sprintf("P*A-2 = %.4f", ifelse(is.na(x$ps$p2),NA,x$ps$p2))),
      div(class="pbox", sprintf("Popt  = %.4f", ifelse(is.na(x$ps$p1),NA,x$ps$p1))),
      div(class="pbox", sprintf("Width = %.4f", ifelse(is.na(w),NA,w)))
    )
  })
  
  output$sb_obs <- renderUI({
    x    <- v()
    cao  <- CA_fn(x$Cse, x$Csp, x$Pobs)
    bgo  <- BG_fn(x$Cse, x$Csp, x$Pobs)
    bgop <- if(!is.na(x$ps$p1)) BG_fn(x$Cse,x$Csp,x$ps$p1) else NA
    gap  <- cao - bgo
    p1   <- if(!is.na(x$ps$p1)) x$ps$p1 else 0
    p2   <- if(!is.na(x$ps$p2)) x$ps$p2 else 1
    over <- if(x$Cse>=x$Csp) x$Pobs>p1 && x$Pobs<p2 else x$Pobs<p1 || x$Pobs>p2
    tagList(
      div(class="mbox", sprintf("Ca(Pobs) = %.4f", cao)),
      div(class="mbox", sprintf("Bg(Pobs) = %.4f", bgo)),
      div(class="mbox", sprintf("Bg(Popt) = %.4f  [max]", ifelse(is.na(bgop),NA,bgop))),
      div(class=if(gap>0)"rbox" else "gbox",
          sprintf("Ca - Bg = %+.4f  (%s)", gap,
                  if(gap>0) "Ca overestimates" else "Bg >= Ca")),
      div(class=if(over)"rbox" else "gbox",
          if(over) "\u26a0 Pobs in overestimation region"
          else    "\u2713 Pobs outside overestimation region")
    )
  })
  
  # ── plot 1 ─────────────────────────────────────────────────────
  output$p1 <- renderPlot({
    make_plot1(v()$Cse, v()$Csp, v()$Pobs)
  }, res=96)
  
  output$dl1p <- downloadHandler("plot1_Ca_vs_Bg.png",
                                 content = function(f) save_plot(make_plot1(v()$Cse,v()$Csp,v()$Pobs), f))
  output$dl1f <- downloadHandler("plot1_Ca_vs_Bg.pdf",
                                 content = function(f) save_plot(make_plot1(v()$Cse,v()$Csp,v()$Pobs), f))
  
  # ── numeric table ──────────────────────────────────────────────
  output$num_table <- renderTable({
    x  <- v()
    m  <- metrics_table(x$Cse, x$Csp, x$Pobs)
    ps <- x$ps
    
    data.frame(
      Metric = c(
        "Classical Accuracy (Ca)",
        "Precision (PPV at Pobs)",
        "Recall = Sensitivity (Cse)",
        "F1-score",
        "MCC",
        "G-mean  [sqrt(Cse * Csp)]",
        "PPV  [Bayesian — at Pobs]",
        "NPV  [Bayesian — at Pobs]",
        "Bg = sqrt(PPV * NPV)  [at Pobs]",
        "Bg at Popt (= P*A-1)",
        "Ca - Bg gap  [at Pobs]"
      ),
      Value = c(
        round(m$CA,    4),
        round(m$prec,  4),
        round(x$Cse,   4),
        round(m$F1,    4),
        round(m$MCC,   4),
        round(m$Gmean, 4),
        round(m$PPV,   4),
        round(m$NPV,   4),
        round(m$BG,    4),
        round(ifelse(!is.na(ps$p1), BG_fn(x$Cse,x$Csp,ps$p1), NA), 4),
        round(m$CA - m$BG, 4)
      ),
      stringsAsFactors = FALSE
    )
  }, striped=TRUE, hover=TRUE, bordered=TRUE, width="100%")
  
  # ── properties table ───────────────────────────────────────────
  output$prop_table <- renderUI({
    yes  <- function(t="Yes") sprintf('<span class="yes">%s</span>', t)
    no   <- function(t="No")  sprintf('<span class="no">%s</span>',  t)
    part <- function(t)       sprintf('<span class="part">%s</span>', t)
    
    rows <- list(
      list("Prevalence as formal argument",
           no(), no(), no(), yes()),
      list("Vanishes at P \u2192 0 and P \u2192 1",
           no(), no(), part("Partially"), yes()),
      list("Anchored in prediction space (PPV / NPV)",
           no(), no(), no(), yes()),
      list("Symmetric under class exchange",
           no("No \u2014 asymmetric"), yes(), yes(), yes()),
      list("Generates algebraic operating interval",
           no(), no(), no(), yes()),
      list("Closed-form prevalence curve",
           no(), no(), no(), yes()),
      list("Geometric mean structure",
           no(), no(), yes("Yes \u2014 intrinsic"), yes("Yes \u2014 predictive")),
      list("Penalizes asymmetry between components",
           part("Partial"), yes(), yes(), yes("Yes \u2014 AM-GM"))
    )
    
    header <- paste0(
      "<table class='ctbl'>",
      "<thead><tr>",
      "<th style='text-align:left;width:34%'>Property</th>",
      "<th>F1</th><th>MCC</th><th>G-mean</th>",
      "<th style='background:#117a65;'>B<sub>G</sub></th>",
      "</tr></thead><tbody>"
    )
    
    body <- paste(sapply(rows, function(r) {
      cells <- paste(sprintf("<td>%s</td>", r), collapse="")
      # first cell left-aligned
      cells <- sub("<td>", "<td style='text-align:left;'>", cells)
      paste0("<tr>", cells, "</tr>")
    }), collapse="")
    
    HTML(paste0(header, body, "</tbody></table>"))
  })
}

shinyApp(ui, server)