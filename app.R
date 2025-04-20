# --- Load Necessary Libraries ---
# shiny: Web application framework for R
# AER: Provides functions for applied econometrics, including ivreg() for IV estimation.
# sandwich: Provides functions for robust covariance matrix estimators.
# lmtest: Provides functions for hypothesis testing, including coeftest() for displaying coefficients with robust SEs.
# DT: Provides an interface to the DataTables library for interactive tables.
# rms: Provides functions for regression modeling strategies, including rcs() for restricted cubic splines.

# Specific References Provided by User:
# 1. Brookhart MA, Rassen JA, Schneeweiss S. Instrumental Variables in Randomized Trials. N Engl J Med. 2010;363(25):e39.
#    Link: https://www.dropbox.com/scl/fi/lsnmm2ge2q1ircugjgqsf/nejm_iv.pdf?rlkey=8xbouk44mnryh2o98fdjbvlfx&raw=1
# 2. Maron DJ, Hochman JS, Reynolds HR, et al. Initial Invasive or Conservative Strategy for Stable Coronary Disease. N Engl J Med. 2020;382(15):1395-1407. (ISCHEMIA Trial)
#    Link: https://www.nejm.org/doi/full/10.1056/NEJMoa1915922

# Install packages if you don't have them (uncomment the lines below)
# install.packages("shiny")
# install.packages("AER")
# install.packages("sandwich")
# install.packages("lmtest")
# install.packages("DT")
# install.packages("rms")

library(shiny)
library(AER)
library(sandwich)
library(lmtest)
library(DT)
library(rms) # Load rms for rcs()

# --- UI Definition (ä½¿ç¨èä»é¢) ---
ui <- fluidPage(
  tags$head(
    tags$meta(charset = "UTF-8"),
    tags$html(lang = "zh-Hant-TW"),
    tags$script(src = "https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.7/MathJax.js?config=TeX-MML-AM_CHTML")
  ),
  titlePanel("æè¨±å¤ééµå¾èä¹é¨æ©åéå°ç§è¨åºè©¦é© RCT çå·¥å·è®æ¸åæ"),
  sidebarLayout(
    sidebarPanel(
      h4("æ§å¶é¢æ¿"),
      h5("æ¨¡æ¬åæ¸è¨­å®:"),
      numericInput("n_sim", "æ¨£æ¬æ¸ (Sample Size, n):", value = 5000, min = 100, max = 20000, step = 100),
      numericInput("beta_T_sim", "çå¯¦æ²»çææ (LATE/CACE):", value = 5.5, min = -10, max = 20, step = 0.5),
      sliderInput("compliance_inv", "éµå¾ç - ä»å¥çµ (P(T=1|Z=1)):", min = 0, max = 1, value = 0.80, step = 0.01),
      sliderInput("compliance_con", "ééµå¾ç - å°ç§çµ (P(T=1|Z=0)):", min = 0, max = 1, value = 0.12, step = 0.01),
      hr(),
      h5("æ¨¡åè¨­å®:"),
      numericInput("rcs_knots", "åºç· SAQ ç RCS çµæ¸ (Number of Knots for RCS):", value = 4, min = 3, max = 7, step = 1),
      helpText("è¨­çº 0 æ 1 åä½¿ç¨ç·æ§é ãå»ºè­° 3-5 åçµã"),
      hr(),
      actionButton("run_analysis", "éæ°æ¨¡æ¬æ¸æä¸¦å·è¡åæ", icon = icon("play"), class = "btn-primary"),
      hr(),
      downloadButton("download_data", "ä¸è¼æ¨¡æ¬æ¸æ (.csv)", class = "btn-success"),
      hr(),
      h5("èªªæ:"),
      p("æ­¤æç¨ç¨å¼æ¼ç¤º RCT ä¸­å­å¨ä¸éµå¾æ§æï¼å¦ä½ä½¿ç¨ IV æ¹æ³ä¼°è¨æ²»çææã"),
      p("æ¨¡æ¬æ¸ææ¨å¨åæ é¡ä¼¼ ISCHEMIA è©¦é©çææ³ã"),
      p("å¯é¸æä½¿ç¨åéä¸æ¬¡æ¨£æ¢ (Restricted Cubic Splines, RCS) ä¾æ¨¡æ¬åºç· SAQ çéç·æ§ææã"),
      p("IV æ¨¡åçµæç¾å¨åå«è¨ºæ·æª¢å®ã")
    ),
    mainPanel(
      tabsetPanel(
        id = "results_tabs",
        tabPanel("åæèªªæèæ¨¡æ¬è¨­å®",
                 h3("å·¥å·è®æ¸ (IV) èé¨æ©è©¦é©"),
                 p("ç¶å­å¨ä¸éµå¾æ§æï¼ç´æ¥æ¯è¼å¯¦éæ¥åæ²»ç (T=1) åæªæ¥åæ²»ç (T=0) çæ£èï¼å³ãä¾æ²»çåæã, OLSï¼å¯è½æå çºãé¸æåèª¤ãèç¢çèª¤å°ï¼å çºæ±ºå®æ¥åæ²»ççå ç´ å¯è½ä¹èçµæç¸éã"),
                 p("å·¥å·è®æ¸ (IV) åæå©ç¨æåçãé¨æ©åæ´¾ã(Z) ä½çºãå·¥å·ãï¼ä¾ä¼°è¨ãå¯¦éæ¥åæ²»çã(T) å°ãçµæã(Y) çå æææãå¨å­å¨ä¸éµå¾æ§ç RCT ä¸­ï¼ééå¸¸ä¼°è¨çæ¯ãéµå¾èçå¹³åå æææã(Complier Average Causal Effect, CACE)ï¼ä¹ç¨±çºãå±é¨å¹³åæ²»çææã(Local Average Treatment Effect, LATE)ã"),
                 hr(),
                 h4("å·¥å·è®æ¸åæçæ ¸å¿åè¨­:"),
                 tags$ul(
                   tags$li("1. éè¯æ§: å·¥å·è®æ¸ Z èå§çè®æ¸ï¼è·èª¤å·®é ç¸éçå³å´è®é ï¼T ç¸é (å¨æ§å¶å±è®æ¸å¾)ãéæ¯å¯ä»¥æª¢é©ç (è¦ãéµå¾æ§æª¢æ¥ãåé )ã"),
                   tags$li("2. ç¨ç«æ§/é¨æ©æ§: å·¥å·è®æ¸ Z èçµæ Y çèª¤å·®é ç¡éï¼äº¦å³ Z æ¯å¤çæ§è®æ¸ãå¨ RCT ä¸­ï¼é¨æ©åæ´¾éå¸¸è½æ»¿è¶³æ­¤åè¨­ã"),
                   tags$li("3. æé¤éå¶ï¼Z âXâYï¼: å·¥å·è®æ¸ Z åªè½ééå§çè®æ¸ T å½±é¿çµæ Yï¼æ²æç´æ¥ç¶ç±æå¶ä»éæ¥è·¯å¾ (å¨æ§å¶å±è®æ¸å¾) å½±é¿ Yãéååè¨­ç¡æ³ç´æ¥æª¢é©ï¼ä½å¨ RCT ä¸­éå¸¸è¢«èªçºæ¯åçç (äº¦å³åæ´¾æ¬èº«ä¸å½±é¿çµæ)ã"),
                   tags$li("4. å®èª¿æ§ï¼æ²æãåæèãï¼å³æ²æäººæå çºè¢«åæ´¾å°ä»å¥çµåèä¸æ¥åæ²»çï¼åæå¦æè¢«åæ´¾å°å°ç§çµåèæ¥åæ²»çãéæå³è Z å° T çå½±é¿æ¹åå°ææäººé½æ¯ä¸è´ç (æèçºé¶)ã")
                 ),
                 hr(),
                 h4("å¨æ­¤æ¨¡æ¬ä¸­ä½¿ç¨çè®æ¸ï¼"),
                 tags$ul(
                   tags$li("Y: çµæè®æ¸ (ä¸å¹´å¾ç SAQ çæ´»åè³ªåæ¸, 0-100, é«å=å¥½)"),
                   tags$li("T: å§çè®æ¸ (å¯¦éæ¥åçè¡ç®¡éå»ºæ²»çï¼1=æ¯, 0=å¦)"),
                   tags$li("Z: å·¥å·è®æ¸ (é¨æ©åæ´¾çæ²»çç­ç¥ï¼1=ä»å¥æ§, 0=ä¿å®æ§)"),
                   tags$li("X: æ§å¶è®æ¸/å±è®æ¸ (æ¨¡æ¬åºç· SAQ åæ¸ `baseline_SAQ`ãå°å `region`)")
                 ),
                 h4("åéä¸æ¬¡æ¨£æ¢ (Restricted Cubic Splines, RCS)"),
                 p("æ¨å¯ä»¥é¸æä½¿ç¨ RCS ä¾å°åºç· SAQ (`baseline_SAQ`) èçµæ Y ææ²»ç T ä¹éçéä¿é²è¡å»ºæ¨¡ãRCS æ¯ä¸ç¨®éæ´»çæ¹å¼ï¼å¯ä»¥å¨ä¸éåº¦æ¬åæ¸æçææ³ä¸æææ½å¨çéç·æ§éä¿ãçµ (knots) çæ¸éæ±ºå®äºæ²ç·çéæ´»æ§ï¼éå¸¸é¸æ 3 å° 5 åçµï¼ã"),
                 hr(),
                 h4("æ¸æçææ¨¡å (Ground Truth Formula for Outcome Y)"),
                 uiOutput("ground_truth_formula_display"),
                 h5("ç®åæ¨¡æ¬ä½¿ç¨çåæ¸å¼:"),
                 verbatimTextOutput("ground_truth_params_display"),
                 hr(),
                 h4("ä¸»è¦åææ¹æ³èå¬å¼ï¼"),
                 tags$ol(
                   tags$li("ç¬¬ä¸éæ®µ (First Stage): T ~ Z + Covariates", uiOutput("formula_first_stage")),
                   tags$li("ç°¡åå¼ (Reduced Form) / ITT: Y ~ Z + Covariates", uiOutput("formula_reduced_form")),
                   tags$li("å©éæ®µæå°å¹³æ¹æ³ (2SLS / IV): Y ~ T + Covariates | Z + Covariates", uiOutput("formula_2sls")),
                   tags$li("æ®éæå°å¹³æ¹æ³ (OLS) / As-Treated: Y ~ T + Covariates", uiOutput("formula_ols"))
                 )
        ),
        tabPanel("æ¸æé è¦½èéµå¾æ§æª¢æ¥",
                 h4("æ¨¡æ¬æ¸æé è¦½ (å 6 ç­)"),
                 tableOutput("data_head"),
                 hr(),
                 h4("ç¬¬ä¸éæ®µè¿´æ­¸ (éµå¾æ§æª¢æ¥ï¼Z å° T çå½±é¿)"),
                 h5("å·¥å·è®æ¸å¼·åº¦èå¼±å·¥å·è®æ¸åé¡"),
                 p("å·¥å·è®æ¸çãéè¯æ§ãåè¨­è¦æ± Z å¿é è T ç¸éãæåéå¸¸ä½¿ç¨ç¬¬ä¸éæ®µè¿´æ­¸ä¸­å·¥å·è®æ¸ Z ç F çµ±è¨éä¾è©ä¼°å·¥å·è®æ¸çå¼·åº¦ã"),
                 tags$ul(
                    tags$li("è¨ç®æ¹å¼ï¼æª¢é©ç¬¬ä¸éæ®µè¿´æ­¸ä¸­ææå·¥å·è®æ¸ï¼å¨æ­¤ä¾ä¸­åªæ Zï¼ä¿æ¸æ¯å¦è¯åçºé¶ç F æª¢å®çµ±è¨éã"),
                    tags$li("ç¶é©æ³åï¼å¸¸ç¨çå¤æ·æ¨æºæ¯ F çµ±è¨é > 10ãéè¡¨ç¤ºå·¥å·è®æ¸èå§çè®æ¸æè¶³å¤ å¼·çéè¯ã"),
                    tags$li("å¼±å·¥å·è®æ¸ï¼è¥ F < 10ï¼åèªçºå·¥å·è®æ¸ Z è T çéè¯æ§è¼å¼±ï¼ç¨±çºãå¼±å·¥å·è®æ¸ãã"),
                    tags$li("å¼±å·¥å·è®æ¸çå¾æï¼"),
                    tags$ul(
                        tags$li("2SLS ä¼°è¨éæç¢çåèª¤ (Bias)ï¼ä¸åèª¤æ¹åææå OLS ä¼°è¨éã"),
                        tags$li("2SLS ä¼°è¨éçæ¨æºèª¤æè®å¾ä¸å¯é ï¼å°è´ä¿¡è³´åéååèªªæª¢å® (å¦ t æª¢å®) ççµæä¸æºç¢ºã"),
                        tags$li("å³ä½¿æ¨£æ¬æ¸å¾å¤§ï¼å¼±å·¥å·è®æ¸åé¡ä»ç¶å¯è½å­å¨ã")
                    ),
                    tags$li("å¦ä½èçå¼±å·¥å·è®æ¸ï¼"),
                    tags$ul(
                        tags$li("å°æ¾æ´å¼·çå·¥å·è®æ¸ (å¦æå¯è½)ã"),
                        tags$li("ä½¿ç¨å°å¼±å·¥å·è®æ¸è¼ä¸ææçä¼°è¨æ¹æ³ï¼ä¾å¦ãæéè³è¨æå¤§æ¦ä¼¼æ³ãLIML)ï¼ä½ AER å¥ä»¶ä¸­ç `ivreg` ä¸ç´æ¥æä¾ LIMLã"),
                        tags$li("ä½¿ç¨å°å¼±å·¥å·è®æ¸ç©©å¥çæ¨è«æ¹æ³ï¼ä¾å¦ Anderson-Rubin æª¢å®ã"),
                        tags$li("æ¿èªå¼±å·¥å·è®æ¸çå­å¨ï¼ä¸¦å¨è§£éçµæææ´å è¬¹æã")
                    ),
                    tags$li("éè¦æ§ï¼æª¢æ¥ç¬¬ä¸éæ®µ F çµ±è¨éæ¯ IV åæä¸­éå¸¸éè¦çä¸æ­¥ã`ivreg` ç `summary` å½æ¸å ä¸ `diagnostics = TRUE` æèªåå ±åæ­¤ F çµ±è¨éã")
                 ),
                 h5("ç·æ§åºç· SAQ æ¨¡åçµæ:"),
                 verbatimTextOutput("first_stage_summary_linear"),
                 h5("RCS åºç· SAQ æ¨¡åçµæ:"),
                 verbatimTextOutput("first_stage_summary_rcs")
        ),
        tabPanel("æ¨¡åçµææ¯è¼",
                 h4("ç°¡åå¼è¿´æ­¸ (ITT ææï¼Z å° Y çå½±é¿)"),
                 h5("ç·æ§åºç· SAQ æ¨¡åçµæ:"),
                 verbatimTextOutput("reduced_form_summary_linear"),
                 h5("RCS åºç· SAQ æ¨¡åçµæ:"),
                 verbatimTextOutput("reduced_form_summary_rcs"),
                 hr(),
                 h4("å·¥å·è®æ¸è¿´æ­¸ (2SLS çµæï¼T å° Y ç LATE/CACE)"),
                 p("æ³¨æï¼ä»¥ä¸çµæåå«ä½¿ç¨ `summary(..., diagnostics = TRUE)` ç¢ççè¨ºæ·æª¢å®ã"),
                 tags$ul(
                    tags$li(strong("å¼±å·¥å·è®æ¸æª¢å®:"), "å ±åç¬¬ä¸éæ®µç F çµ±è¨é (åãéµå¾æ§æª¢æ¥ãåé )ãF > 10 éå¸¸è¡¨ç¤ºå·¥å·è®æ¸å¼·åº¦è¶³å¤ ã"),
                    tags$li(strong("å§çæ§æª¢å® (Wu-Hausman Test):"), "æª¢å®å§çè®æ¸ (T) æ¯å¦ççå·æå§çæ§ (å³ OLS æ¯å¦æåèª¤)ãèç¡åè¨­ (H0) æ¯ T çºå¤çãè¥ p å¼é¡¯è (ä¾å¦ < 0.05)ï¼åæçµ H0ï¼æ¯æ T å·æå§çæ§ï¼ä½¿ç¨ IV æ¯é©ç¶çã"),
                    tags$li(strong("éåº¦è­å¥æª¢å® (Sargan Test of Overidentifying Restrictions):"), "åå¨å·¥å·è®æ¸æ¸é > å§çè®æ¸æ¸éæé©ç¨ (å¨æ­¤æ¨¡æ¬ä¸­ä¸é©ç¨ï¼å çºåªæä¸åå·¥å· Z åä¸åå§çè®æ¸ T)ãæª¢å®ãé¡å¤çãå·¥å·è®æ¸æ¯å¦ææ (å³æ¯å¦æ»¿è¶³æé¤éå¶ä¸èèª¤å·®é ç¡é)ãèç¡åè¨­ (H0) æ¯ææå·¥å·è®æ¸é½ææãè¥ p å¼ä¸é¡¯èï¼åç¡æ³æçµ H0ï¼æ¯æå·¥å·è®æ¸çæææ§ã")
                 ),
                 h5("ç·æ§åºç· SAQ æ¨¡åçµæ (å«è¨ºæ·):"),
                 verbatimTextOutput("iv_summary_linear"),
                 h5("RCS åºç· SAQ æ¨¡åçµæ (å«è¨ºæ·):"),
                 verbatimTextOutput("iv_summary_rcs"),
                 hr(),
                 h4("å·¥å·è®æ¸è¿´æ­¸ (å«ç©©å¥æ¨æºèª¤)"),
                 p("æ³¨æï¼ä½¿ç¨ `coeftest` æ­é `vcovHC` è¨ç®ç©©å¥æ¨æºèª¤æï¼ä¸æé¡¯ç¤º `summary` çè¨ºæ·æª¢å®ãè¨ºæ·æª¢å®æåèä¸é¢çæ¨æº IV çµæã"),
                 h5("ç·æ§åºç· SAQ æ¨¡åçµæ (ç©©å¥ SE):"),
                 verbatimTextOutput("iv_robust_summary_linear"),
                 h5("RCS åºç· SAQ æ¨¡åçµæ (ç©©å¥ SE):"),
                 verbatimTextOutput("iv_robust_summary_rcs"),
                 hr(),
                 h4("æ®éæå°å¹³æ¹æ³è¿´æ­¸ (OLS / As-Treated)"),
                 h5("ç·æ§åºç· SAQ æ¨¡åçµæ:"),
                 verbatimTextOutput("ols_summary_linear"),
                 h5("RCS åºç· SAQ æ¨¡åçµæ:"),
                 verbatimTextOutput("ols_summary_rcs"),
                 hr(),
                 h4("æ¨¡åä¼°è¨å¼æ¯è¼ (åºæ¼ç·æ§ SAQ æ¨¡å)"),
                 uiOutput("comparison_estimates")
        ),
        tabPanel("å®æ´æ¨¡æ¬æ¸æ",
                 h4("å®æ´çæ¨¡æ¬æ¸æé"),
                 DT::dataTableOutput("full_data_table")
        ),
        # --- NEW: References Tab ---
        tabPanel("åèæç»",
                 h4("ä¸»è¦åèæç»"),
                 tags$ul(
                   tags$li(
                     tags$a(href = "https://www.dropbox.com/scl/fi/lsnmm2ge2q1ircugjgqsf/nejm_iv.pdf?rlkey=8xbouk44mnryh2o98fdjbvlfx&raw=1",
                            target = "_blank", # Open in new tab
                            "Brookhart MA, Rassen JA, Schneeweiss S. Instrumental Variables in Randomized Trials. N Engl J Med. 2010;363(25):e39.")
                   ),
                   tags$li(
                     tags$a(href = "https://www.nejm.org/doi/full/10.1056/NEJMoa1915922",
                            target = "_blank", # Open in new tab
                            "Maron DJ, Hochman JS, Reynolds HR, et al. Initial Invasive or Conservative Strategy for Stable Coronary Disease. N Engl J Med. 2020;382(15):1395-1407. (ISCHEMIA Trial)")
                   ) 
                 )
        )
        # --- End of New Tab ---
      )
    )
  )
)


# --- Server Logic (ä¼ºæå¨éè¼¯) ---
server <- function(input, output, session) {

  analysis_results <- eventReactive(input$run_analysis, {

    showNotification("æ­£å¨æ¨¡æ¬æ¸æä¸¦å·è¡åæ...", type = "message", duration = 3)

    # --- 1. Simulate Data ---
    set.seed(as.integer(Sys.time()))
    n <- input$n_sim
    beta_T_true_late <- input$beta_T_sim
    prob_revasc_if_Z1 <- input$compliance_inv
    prob_revasc_if_Z0 <- input$compliance_con
    nk <- input$rcs_knots # Get number of knots for RCS

    # Z: Random assignment
    assignment_Z <- rbinom(n, 1, 0.5)

    # X: Covariates
    baseline_SAQ <- rnorm(n, mean = 60, sd = 15)
    region <- factor(sample(1:3, n, replace = TRUE), labels = c("å°åA", "å°åB", "å°åC"))

    # T: Treatment Received (Endogenous)
    # Introduce some correlation between baseline_SAQ and compliance to make OLS biased
    compliance_propensity <- pnorm(-1.5 + 0.5*assignment_Z + 0.02*(baseline_SAQ - 60)) # Higher SAQ slightly more likely to comply if Z=1
    prob_revasc_adj <- ifelse(assignment_Z == 1,
                              pmin(pmax(prob_revasc_if_Z1 + (compliance_propensity - mean(compliance_propensity[assignment_Z==1])), 0.01), 0.99),
                              pmin(pmax(prob_revasc_if_Z0 + (compliance_propensity - mean(compliance_propensity[assignment_Z==0])), 0.01), 0.99))

    treatment_revasc_T <- rbinom(n, 1, prob_revasc_adj)


    # Y: Outcome (Final SAQ Score) - Ground Truth (linear baseline effect)
    intercept <- 76
    beta_baseline <- 0.1 # True effect is linear for baseline SAQ
    beta_regionB <- 2
    beta_regionC <- -1
    error_sd <- 20
    # Add confounding: Higher baseline SAQ also independently leads to better outcome
    confounding_effect <- 0.05 * (baseline_SAQ - 60)

    outcome_Y <- intercept +
                 beta_T_true_late * treatment_revasc_T +
                 beta_baseline * baseline_SAQ + # Effect of baseline SAQ
                 ifelse(region == "å°åB", beta_regionB, 0) +
                 ifelse(region == "å°åC", beta_regionC, 0) +
                 # confounding_effect + # Add confounding related to baseline_SAQ
                 rnorm(n, mean = 0, sd = error_sd)

    # Create data frame
    ischemia_sim <- data.frame(
      Y = outcome_Y,
      T = treatment_revasc_T,
      Z = assignment_Z,
      baseline_SAQ = baseline_SAQ,
      region = region
    )

    # --- 2. Perform Analyses ---

    # Define covariate terms strings
    formula_linear_cov_terms <- "baseline_SAQ + region"
    formula_rcs_cov_terms <- if (nk >= 3) paste("rcs(baseline_SAQ, nk=", nk, ") + region") else formula_linear_cov_terms

    # Define formulas (Linear and RCS)
    formula_fs_linear <- as.formula(paste("T ~ Z +", formula_linear_cov_terms))
    formula_rf_linear <- as.formula(paste("Y ~ Z +", formula_linear_cov_terms))
    formula_ols_linear <- as.formula(paste("Y ~ T +", formula_linear_cov_terms))
    formula_iv_linear <- as.formula(paste("Y ~ T +", formula_linear_cov_terms, "| Z +", formula_linear_cov_terms))

    formula_fs_rcs <- if (nk >= 3) as.formula(paste("T ~ Z +", formula_rcs_cov_terms)) else NULL
    formula_rf_rcs <- if (nk >= 3) as.formula(paste("Y ~ Z +", formula_rcs_cov_terms)) else NULL
    formula_ols_rcs <- if (nk >= 3) as.formula(paste("Y ~ T +", formula_rcs_cov_terms)) else NULL
    formula_iv_rcs <- if (nk >= 3) as.formula(paste("Y ~ T +", formula_rcs_cov_terms, "| Z +", formula_rcs_cov_terms)) else NULL


    # --- Run Linear Models ---
    first_stage_linear_model <- lm(formula_fs_linear, data = ischemia_sim)
    first_stage_linear_summary <- summary(first_stage_linear_model) # Summary for F-stat later

    reduced_form_linear_model <- lm(formula_rf_linear, data = ischemia_sim)
    reduced_form_linear_summary <- summary(reduced_form_linear_model)

    ols_linear_model <- lm(formula_ols_linear, data = ischemia_sim)
    ols_linear_summary <- summary(ols_linear_model)

    iv_linear_model <- tryCatch({
        ivreg(formula_iv_linear, data = ischemia_sim)
      }, error = function(e) { message("IV Linear Model Error: ", e$message); NULL })

    # **MODIFIED**: Use diagnostics = TRUE in summary() for ivreg object
    iv_linear_summary <- if (!is.null(iv_linear_model) && inherits(iv_linear_model, "ivreg")) {
        tryCatch({
            summary(iv_linear_model, diagnostics = TRUE)
        }, error = function(e) { message("IV Linear Summary Error: ", e$message); "IV (Linear) æ¨¡åæè¦å¤±æ (å«è¨ºæ·)"})
    } else {
        "IV (Linear) æ¨¡åä¼°è¨å¤±æ"
    }


    iv_linear_robust_summary <- NULL
    if (!is.null(iv_linear_model) && inherits(iv_linear_model, "ivreg")) {
      iv_linear_robust_summary <- tryCatch({
          coeftest(iv_linear_model, vcov. = vcovHC(iv_linear_model, type = "HC1"))
        }, error = function(e) { message("Robust SE (Linear) Error: ", e$message); NULL })
      if(is.null(iv_linear_robust_summary)) iv_linear_robust_summary <- "ç¡æ³è¨ç®ç©©å¥æ¨æºèª¤ (Linear)ã"
    } else {
        iv_linear_robust_summary <- "ç¡æ³è¨ç®ç©©å¥æ¨æºèª¤ (Linear)ã"
    }

    # --- Run RCS Models (if nk >= 3) ---
    first_stage_rcs_model <- NULL
    first_stage_rcs_summary <- "æªå·è¡ RCS æ¨¡å (çµæ¸ < 3)ã"
    reduced_form_rcs_model <- NULL
    reduced_form_rcs_summary <- "æªå·è¡ RCS æ¨¡å (çµæ¸ < 3)ã"
    ols_rcs_model <- NULL
    ols_rcs_summary <- "æªå·è¡ RCS æ¨¡å (çµæ¸ < 3)ã"
    iv_rcs_model <- NULL
    iv_rcs_summary <- "æªå·è¡ RCS æ¨¡å (çµæ¸ < 3)ã"
    iv_rcs_robust_summary <- "æªå·è¡ RCS æ¨¡å (çµæ¸ < 3)ã"

    if (!is.null(formula_fs_rcs)) {
        first_stage_rcs_model <- tryCatch({ lm(formula_fs_rcs, data = ischemia_sim) }, error = function(e){ message("FS RCS Error: ", e$message); NULL})
        first_stage_rcs_summary <- if (!is.null(first_stage_rcs_model)) summary(first_stage_rcs_model) else "ç¬¬ä¸éæ®µ (RCS) æ¨¡åä¼°è¨å¤±æ" # Summary for F-stat later

        reduced_form_rcs_model <- tryCatch({ lm(formula_rf_rcs, data = ischemia_sim) }, error = function(e){ message("RF RCS Error: ", e$message); NULL})
        reduced_form_rcs_summary <- if (!is.null(reduced_form_rcs_model)) summary(reduced_form_rcs_model) else "ç°¡åå¼ (RCS) æ¨¡åä¼°è¨å¤±æ"

        ols_rcs_model <- tryCatch({ lm(formula_ols_rcs, data = ischemia_sim) }, error = function(e){ message("OLS RCS Error: ", e$message); NULL})
        ols_rcs_summary <- if (!is.null(ols_rcs_model)) summary(ols_rcs_model) else "OLS (RCS) æ¨¡åä¼°è¨å¤±æ"

        iv_rcs_model <- tryCatch({
            ivreg(formula_iv_rcs, data = ischemia_sim)
        }, error = function(e) { message("IV RCS Model Error: ", e$message); NULL })

        # **MODIFIED**: Use diagnostics = TRUE in summary() for ivreg object
        iv_rcs_summary <- if (!is.null(iv_rcs_model) && inherits(iv_rcs_model, "ivreg")) {
             tryCatch({
                summary(iv_rcs_model, diagnostics = TRUE)
             }, error = function(e) { message("IV RCS Summary Error: ", e$message); "IV (RCS) æ¨¡åæè¦å¤±æ (å«è¨ºæ·)"})
        } else {
             "IV (RCS) æ¨¡åä¼°è¨å¤±æ"
        }


        if (!is.null(iv_rcs_model) && inherits(iv_rcs_model, "ivreg")) {
          iv_rcs_robust_summary <- tryCatch({
              coeftest(iv_rcs_model, vcov. = vcovHC(iv_rcs_model, type = "HC1"))
            }, error = function(e) { message("Robust SE (RCS) Error: ", e$message); NULL })
           if(is.null(iv_rcs_robust_summary)) iv_rcs_robust_summary <- "ç¡æ³è¨ç®ç©©å¥æ¨æºèª¤ (RCS)ã"
        } else {
            iv_rcs_robust_summary <- "ç¡æ³è¨ç®ç©©å¥æ¨æºèª¤ (RCS)ã"
        }
    }


    # Store ground truth parameters
    ground_truth_params <- list(
        n = n,
        intercept = intercept,
        beta_T_true_late = beta_T_true_late,
        beta_baseline = beta_baseline,
        beta_regionB = beta_regionB,
        beta_regionC = beta_regionC,
        error_sd = error_sd,
        prob_revasc_if_Z1 = prob_revasc_if_Z1,
        prob_revasc_if_Z0 = prob_revasc_if_Z0
    )

    # Return all results
    list(
      data_full = ischemia_sim,
      # Linear models
      first_stage_linear_model = first_stage_linear_model, # Pass model object
      first_stage_linear_summary = first_stage_linear_summary, # Pass summary object
      reduced_form_linear_model = reduced_form_linear_model,
      reduced_form_linear_summary = reduced_form_linear_summary,
      ols_linear_model = ols_linear_model,
      ols_linear_summary = ols_linear_summary,
      iv_linear_model = iv_linear_model, # Pass model object
      iv_linear_summary = iv_linear_summary, # Pass summary object (now includes diagnostics)
      iv_linear_robust_summary = iv_linear_robust_summary,
      # RCS models
      first_stage_rcs_model = first_stage_rcs_model,
      first_stage_rcs_summary = first_stage_rcs_summary,
      reduced_form_rcs_model = reduced_form_rcs_model,
      reduced_form_rcs_summary = reduced_form_rcs_summary,
      ols_rcs_model = ols_rcs_model,
      ols_rcs_summary = ols_rcs_summary,
      iv_rcs_model = iv_rcs_model,
      iv_rcs_summary = iv_rcs_summary, # Pass summary object (now includes diagnostics)
      iv_rcs_robust_summary = iv_rcs_robust_summary,
      # Params
      ground_truth_params = ground_truth_params,
      nk = nk
    )
  }, ignoreNULL = FALSE)

  # --- 3. Render Outputs ---

  output$data_head <- renderTable({
    results <- analysis_results()
    req(results$data_full)
    head(results$data_full)
  }, rownames = TRUE)

  output$full_data_table <- DT::renderDataTable({
    results <- analysis_results()
    req(results$data_full)
    DT::datatable(results$data_full,
                  options = list(pageLength = 10, scrollX = TRUE, searching = TRUE),
                  rownames = FALSE,
                  filter = 'top')
  })

  # Helper function to print summary and F-stat for FIRST STAGE (summary.lm)
  print_first_stage_summary <- function(summary_obj, model_type = "Linear") {
      if (is.character(summary_obj)) { # Handle error messages
          cat(paste("ç¬¬ä¸éæ®µ (", model_type, "): ", summary_obj, "\n", sep=""))
          return()
      }
      # Check if summary object is valid before printing
      req(summary_obj, inherits(summary_obj, "summary.lm"))
      cat(paste("--- ç¬¬ä¸éæ®µæ¨¡å (", model_type, ") ---\n", sep=""))
      print(summary_obj)

      # Extract and display F-statistic for the instrument Z
      # We need the F-stat for the exclusion of the instrument(s)
      # For a single instrument Z, this is equivalent to the square of the t-statistic for Z
      if ("Z" %in% rownames(summary_obj$coefficients)) {
          t_stat_Z <- summary_obj$coefficients["Z", "t value"]
          f_stat_Z <- t_stat_Z^2
          cat("\n---\n")
          cat(paste("å·¥å·è®æ¸ Z ç F çµ±è¨é (è¿ä¼¼å¼ t^2):", round(f_stat_Z, 2), "\n"))
          if (f_stat_Z < 10) {
              cat("è­¦åï¼ç¬¬ä¸éæ®µ F çµ±è¨é < 10ãå¯è½å­å¨å¼±å·¥å·è®æ¸åé¡ã\n")
          } else {
              cat("ç¬¬ä¸éæ®µ F çµ±è¨é >= 10ï¼å·¥å·è®æ¸å¼·åº¦å°å¯æ¥åã\n")
          }
      } else {
          cat("\n---\nç¡æ³æåå·¥å·è®æ¸ Z ç F çµ±è¨éã\n")
      }
  }

  # First Stage Output
  output$first_stage_summary_linear <- renderPrint({
      results <- analysis_results()
      req(results$first_stage_linear_summary)
      print_first_stage_summary(results$first_stage_linear_summary, "Linear")
  })
  output$first_stage_summary_rcs <- renderPrint({
      results <- analysis_results()
      req(results$first_stage_rcs_summary)
      print_first_stage_summary(results$first_stage_rcs_summary, "RCS")
  })

  # Reduced Form Output
  output$reduced_form_summary_linear <- renderPrint({
      results <- analysis_results()
      req(results$reduced_form_linear_summary)
      cat("--- ç°¡åå¼æ¨¡å (Linear) ---\n")
      if(is.character(results$reduced_form_linear_summary)) {cat(results$reduced_form_linear_summary)} else {print(results$reduced_form_linear_summary)}
  })
  output$reduced_form_summary_rcs <- renderPrint({
      results <- analysis_results()
      req(results$reduced_form_rcs_summary)
      cat("--- ç°¡åå¼æ¨¡å (RCS) ---\n")
      if(is.character(results$reduced_form_rcs_summary)) {cat(results$reduced_form_rcs_summary)} else {print(results$reduced_form_rcs_summary)}
  })

  # IV / 2SLS Output (Now includes diagnostics)
  output$iv_summary_linear <- renderPrint({
    results <- analysis_results()
    req(results$iv_linear_summary)
    cat("--- IV/2SLS æ¨¡å (Linear, å«è¨ºæ·) ---\n")
    # The summary object itself now contains diagnostics when printed
    if (is.character(results$iv_linear_summary)) { cat(results$iv_linear_summary) } else { print(results$iv_linear_summary) }
  })
  output$iv_summary_rcs <- renderPrint({
    results <- analysis_results()
    req(results$iv_rcs_summary)
    cat("--- IV/2SLS æ¨¡å (RCS, å«è¨ºæ·) ---\n")
     if (is.character(results$iv_rcs_summary)) { cat(results$iv_rcs_summary) } else { print(results$iv_rcs_summary) }
  })

  # IV / 2SLS Robust Output (No diagnostics here)
  output$iv_robust_summary_linear <- renderPrint({
     results <- analysis_results()
     req(results$iv_linear_robust_summary)
     cat("--- IV/2SLS æ¨¡å (Linear, ç©©å¥ SE) ---\n")
     if (is.character(results$iv_linear_robust_summary)) { cat(results$iv_linear_robust_summary) } else { print(results$iv_linear_robust_summary) }
  })
   output$iv_robust_summary_rcs <- renderPrint({
     results <- analysis_results()
     req(results$iv_rcs_robust_summary)
     cat("--- IV/2SLS æ¨¡å (RCS, ç©©å¥ SE) ---\n")
     if (is.character(results$iv_rcs_robust_summary)) { cat(results$iv_rcs_robust_summary) } else { print(results$iv_rcs_robust_summary) }
  })

  # OLS Output
  output$ols_summary_linear <- renderPrint({
     results <- analysis_results()
     req(results$ols_linear_summary)
     cat("--- OLS æ¨¡å (Linear) ---\n")
     if(is.character(results$ols_linear_summary)) {cat(results$ols_linear_summary)} else {print(results$ols_linear_summary)}
  })
  output$ols_summary_rcs <- renderPrint({
     results <- analysis_results()
     req(results$ols_rcs_summary)
     cat("--- OLS æ¨¡å (RCS) ---\n")
     if(is.character(results$ols_rcs_summary)) {cat(results$ols_rcs_summary)} else {print(results$ols_rcs_summary)}
  })

  # --- Render Formulas ---
  output$ground_truth_formula_display <- renderUI({
      withMathJax(helpText(
        "$$ Y_i = \\beta_0 + \\beta_{T, LATE} T_i + \\mathbf{X}_i'\\boldsymbol{\\beta}_X + \\epsilon_i $$"
      ))
  })
  output$formula_first_stage <- renderUI({
      withMathJax(helpText(
        "$$ T_i = \\gamma_0 + \\gamma_Z Z_i + \\mathbf{X}_i'\\boldsymbol{\\gamma}_X + \\nu_i $$"
      ))
  })
  output$formula_reduced_form <- renderUI({
      withMathJax(helpText(
        "$$ Y_i = \\pi_0 + \\pi_Z Z_i + \\mathbf{X}_i'\\boldsymbol{\\pi}_X + \\omega_i $$"
      ))
  })
  output$formula_2sls <- renderUI({
      # Representing the two stages implicitly
      withMathJax(helpText(
        "$$ Y_i = \\beta_0 + \\beta_{T, 2SLS} \\hat{T}_i + \\mathbf{X}_i'\\boldsymbol{\\beta}_X + e_i \\quad \\text{where } \\hat{T}_i \\text{ is predicted from the first stage.} $$"
      ))
  })
  output$formula_ols <- renderUI({
      withMathJax(helpText(
        "$$ Y_i = \\delta_0 + \\delta_T T_i + \\mathbf{X}_i'\\boldsymbol{\\delta}_X + \\mu_i $$"
      ))
  })

  # Render ground truth parameters
  output$ground_truth_params_display <- renderPrint({
      results <- analysis_results()
      req(results$ground_truth_params)
      params <- results$ground_truth_params
      prop_compliers = params$prob_revasc_if_Z1 - params$prob_revasc_if_Z0
      cat(paste("æ¨£æ¬æ¸ (n):", params$n, "\n"))
      cat(paste("çå¯¦ LATE/CACE (Î²_T):", params$beta_T_true_late, "\n"))
      cat(paste("çå¯¦åºç· SAQ ææ (Î²_baseline, ç·æ§):", params$beta_baseline, "\n"))
      cat(paste("å°å B ææ (Î²_regionB):", params$beta_regionB, "\n"))
      cat(paste("å°å C ææ (Î²_regionC):", params$beta_regionC, "\n"))
      cat(paste("èª¤å·®æ¨æºå·® (Ï_Îµ):", params$error_sd, "\n"))
      cat(paste("P(T=1 | Z=1) (å¹³å):", params$prob_revasc_if_Z1, "\n")) # Note: This is the input slider value, actual average might differ slightly due to adjustment
      cat(paste("P(T=1 | Z=0) (å¹³å):", params$prob_revasc_if_Z0, "\n")) # Note: This is the input slider value, actual average might differ slightly due to adjustment
      cat(paste("éµå¾èæ¯ä¾ä¼°è¨ (åºæ¼è¼¸å¥):", round(prop_compliers, 3), "\n"))
  })

  # Comparison uses linear models
  output$comparison_estimates <- renderUI({
      results <- analysis_results()
      # Ensure linear models ran and coefficients can be extracted
      req(results$reduced_form_linear_model, inherits(results$reduced_form_linear_model, "lm"),
          results$iv_linear_model, inherits(results$iv_linear_model, "ivreg"),
          results$ols_linear_model, inherits(results$ols_linear_model, "lm"),
          results$ground_truth_params)

      # Safely extract coefficients
      itt_coef <- tryCatch(coef(results$reduced_form_linear_model)["Z"], error = function(e) NA_real_)
      iv_coef <- tryCatch(coef(results$iv_linear_model)["T"], error = function(e) NA_real_)
      ols_coef <- tryCatch(coef(results$ols_linear_model)["T"], error = function(e) NA_real_)
      true_late <- results$ground_truth_params$beta_T_true_late

      # Check if extraction failed
      if (anyNA(c(itt_coef, iv_coef, ols_coef))) {
          return(p("é¯èª¤ï¼ç¡æ³å¾ä¸åæå¤åç·æ§æ¨¡åä¸­æåä¿æ¸é²è¡æ¯è¼ãè«æª¢æ¥æ¨¡åçµæã"))
      }

      # Create comparison text
      tagList(
          p(paste0("å¨æ­¤æ¨¡æ¬ä¸­ï¼çå¯¦çéµå¾èå¹³åå æææ (LATE/CACE) è¨­å®çº: ", round(true_late, 3))),
          hr(),
          h5("ä¸ååææ¹æ³çä¼°è¨å¼ (åºæ¼ç·æ§ SAQ æ¨¡å)ï¼"),
          tags$ul(
              tags$li(paste0("ITT (æåæ²»ç) ææ (Z å° Y): ", round(itt_coef, 3), ". é¨æ©ãåæ´¾ãå°çµæ Y çå¹³åå½±é¿ã")),
              tags$li(paste0("OLS (ä¾æ²»çåæ) éè¯æ§ (T å° Y): ", round(ols_coef, 3), ". æ¯è¼å¯¦éãæ¥åãæ²»çèççµæå·®ç°ï¼å¯è½æåèª¤ã")),
              tags$li(paste0("2SLS (å·¥å·è®æ¸) ææ (T å° Y ç LATE/CACE): ", round(iv_coef, 3), ". å©ç¨é¨æ©åæ´¾ Z ä¼°è¨å¯¦éãæ¥åãæ²»ç T å°ãéµå¾èãçå æææã"))
          ),
          hr(),
          h5("æ¯è¼èè§£éï¼"),
          tags$ul(
              tags$li(paste0("2SLS ä¼°è¨å¼ (", round(iv_coef, 3), ") éå¸¸ææ¥è¿çå¯¦ LATE (", round(true_late, 3), ")ï¼å å¶è©¦åæ ¡æ­£ç± T çå§çæ§ï¼é¸æåèª¤ï¼å¼èµ·çåèª¤ã")),
              tags$li(paste0("OLS ä¼°è¨å¼ (", round(ols_coef, 3), ") èçå¯¦ LATE çå·®ç°åæ äºãéµå¾è¡çºãèçµæä¹éçæ··æ·ãå¨æ­¤æ¨¡æ¬ä¸­ï¼åºç· SAQ å¯è½åæå½±é¿æ²»çé¸æåçµæï¼å°è´ OLS åèª¤ã")),
              tags$li(paste0("ITT ä¼°è¨å¼ (", round(itt_coef, 3), ") ä»£è¡¨æ²»çãç­ç¥ãçå¹³åææï¼å¶å¤§å°åéµå¾èæ¯ä¾å½±é¿ãçè«ä¸ï¼ITT â LATE Ã (P(T=1|Z=1) - P(T=1|Z=0))ã"))
          ),
          hr(),
          p("æ³¨æï¼ä»¥ä¸æ¯è¼åºæ¼åè¨­åºç· SAQ ææçºç·æ§çæ¨¡åãæ¨å¯ä»¥æ¥çåå« RCS çæ¨¡åçµæï¼ä»¥è©ä¼°éç·æ§åè¨­æ¯å¦å½±é¿ä¼°è¨ãåæï¼è«æª¢æ¥ IV æ¨¡åçè¨ºæ·æª¢å® (å¦å¼±å·¥å·è®æ¸æª¢å®)ã")
      )
  })

  # --- Download Handler ---
  output$download_data <- downloadHandler(
    filename = function() {
      results <- analysis_results()
      req(results$nk)
      paste0("simulated_ischemia_data_", Sys.Date(), "_nk", results$nk, ".csv")
    },
    content = function(file) {
      results <- analysis_results()
      if (!is.null(results$data_full)) {
        write.csv(results$data_full, file, row.names = FALSE, fileEncoding = "UTF-8")
      } else {
        # Provide a minimal CSV if data generation failed
        write.csv(data.frame(Error = "No simulation data generated yet."), file, row.names = FALSE)
      }
    }
  )

}

# --- Run the application ---
shinyApp(ui = ui, server = server)
