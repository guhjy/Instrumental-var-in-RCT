# --- Load Necessary Libraries ---
# shiny: Web application framework for R
# AER: Provides functions for applied econometrics, including ivreg() for IV estimation.
# sandwich: Provides functions for robust covariance matrix estimators.
# lmtest: Provides functions for hypothesis testing, including coeftest() for displaying coefficients with robust SEs.
# DT: Provides an interface to the DataTables library for interactive tables.
# rms: Provides functions for regression modeling strategies, including rcs() for restricted cubic splines.

# --- References ---
# General IV Concepts in RCTs:
# Hernán MA, Robins JM. Instruments for Causal Inference: An Epidemiologist’s Dream? Epidemiology. 2006;17(4):360-372.
# Baiocchi M, Cheng J, Small DS. Instrumental Variable Methods for Causal Inference. Stat Med. 2014;33(13):2297-2340.
# Weak Instrument / F-statistic Rule of Thumb:
# Staiger D, Stock JH. Instrumental Variables Regression with Weak Instruments. Econometrica. 1997;65(3):557-586.
# Stock JH, Yogo M. Testing for Weak Instruments in Linear IV Regression. In: Andrews DWK, Stock JH, eds. Identification and Inference for Econometric Models: Essays in Honor of Thomas Rothenberg. Cambridge University Press; 2005:80-108.

# Specific References Provided by User:
# 1. Brookhart MA, Rassen JA, Schneeweiss S. Instrumental Variables in Randomized Trials. N Engl J Med. 2010;363(25):e39. (Equivalent to Ref 25 from PDF link)
#    Link: https://www.dropbox.com/scl/fi/lsnmm2ge2q1ircugjgqsf/nejm_iv.pdf?rlkey=8xbouk44mnryh2o98fdjbvlfx&raw=1
# 2. Maron DJ, Hochman JS, Reynolds HR, et al. Initial Invasive or Conservative Strategy for Stable Coronary Disease. N Engl J Med. 2020;382(15):1395-1407. (ISCHEMIA Trial - Ref 20)
#    Link: https://www.nejm.org/doi/full/10.1056/NEJMoa1915922

# Install packages if you don't have them (uncomment the lines below)
# install.packages("shiny")
# install.packages("AER")
# install.packages("sandwich")
# install.packages("lmtest")
# install.packages("DT")
# install.packages("rms") # Added rms

library(shiny)
library(AER)
library(sandwich)
library(lmtest)
library(DT)
library(rms) # Load rms for rcs()

# --- UI Definition (使用者介面) ---
ui <- fluidPage(
  tags$head(
    tags$meta(charset = "UTF-8"),
    tags$html(lang = "zh-Hant-TW"),
    tags$script(src = "https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.7/MathJax.js?config=TeX-MML-AM_CHTML")
  ),
  titlePanel("工具變數分析 Shiny 應用程式 (含 RCS 選項)"),
  sidebarLayout(
    sidebarPanel(
      h4("控制面板"),
      h5("模擬參數設定:"),
      numericInput("n_sim", "樣本數 (Sample Size, n):", value = 5000, min = 100, max = 20000, step = 100),
      numericInput("beta_T_sim", "真實治療效應 (LATE/CACE):", value = 5.5, min = -10, max = 20, step = 0.5),
      sliderInput("compliance_inv", "遵從率 - 介入組 (P(T=1|Z=1)):", min = 0, max = 1, value = 0.80, step = 0.01),
      sliderInput("compliance_con", "非遵從率 - 對照組 (P(T=1|Z=0)):", min = 0, max = 1, value = 0.12, step = 0.01),
      hr(),
      # NEW: Input for RCS knots
      h5("模型設定:"),
      numericInput("rcs_knots", "基線 SAQ 的 RCS 結數 (Number of Knots for RCS):", value = 4, min = 3, max = 7, step = 1),
      helpText("設為 0 或 1 則使用線性項。建議 3-5 個結。"),
      hr(),
      actionButton("run_analysis", "重新模擬數據並執行分析", icon = icon("play"), class = "btn-primary"),
      hr(),
      downloadButton("download_data", "下載模擬數據 (.csv)", class = "btn-success"),
      hr(),
      h5("說明:"),
      p("此應用程式演示 RCT 中存在不遵從性時，如何使用 IV 方法估計治療效應。"),
      p("模擬數據旨在反映類似 ISCHEMIA 試驗的情況。"),
      p("可選擇使用受限三次樣條 (Restricted Cubic Splines, RCS) 來模擬基線 SAQ 的非線性效應。")
    ),
    mainPanel(
      tabsetPanel(
        id = "results_tabs",
        tabPanel("分析說明與模擬設定",
                 h3("工具變數 (IV) 與隨機試驗"),
                 # ... (existing explanation text remains largely the same) ...
                 p("當存在不遵從性時，直接比較實際接受治療 (T=1) 和未接受治療 (T=0) 的患者（即『依治療分析』, as-treated analysis using OLS）可能會因為『選擇偏誤』而產生誤導..."),
                 p("工具變數 (IV) 分析利用最初的「隨機分派」(Z) 作為「工具」，來估計「實際接受治療」(T) 對「結果」(Y) 的因果效應 (LATE/CACE)..."),
                 hr(),
                 h4("工具變數分析的核心假設:"),
                 # ... (assumptions list) ...
                 hr(),
                 h4("在此模擬中使用的變數："),
                 tags$ul(
                   tags$li("Y: 結果變數 (模擬一年後的 SAQ 生活品質分數, 0-100, 高分=好)"),
                   tags$li("T: 內生變數 (實際接受的血管重建治療，1=是, 0=否)"),
                   tags$li("Z: 工具變數 (隨機分派的治療策略，1=介入性, 0=保守性)"),
                   tags$li("X: 控制變數/共變數 (模擬基線 SAQ 分數 `baseline_SAQ`、地區 `region`)")
                 ),
                 # NEW: Explanation of RCS
                 h4("受限三次樣條 (Restricted Cubic Splines, RCS)"),
                 p("您可以選擇使用 RCS 來對基線 SAQ (`baseline_SAQ`) 與結果 Y 或治療 T 之間的關係進行建模。RCS 是一種靈活的方式，可以在不過度擬合數據的情況下捕捉潛在的非線性關係。它假設關係在第一個結之前和最後一個結之後是線性的。結 (knots) 的數量決定了曲線的靈活性（通常選擇 3 到 5 個結）。"),
                 hr(),
                 h4("數據生成模型 (Ground Truth Formula for Outcome Y)"),
                 uiOutput("ground_truth_formula_display"),
                 h5("目前模擬使用的參數值:"),
                 verbatimTextOutput("ground_truth_params_display"),
                 hr(),
                 h4("主要分析方法與公式："),
                 # ... (formulas remain the same, but interpretation changes slightly if RCS is used) ...
                 tags$ol(
                   tags$li("第一階段 (First Stage): T ~ Z + X (或 rcs(X))", uiOutput("formula_first_stage")),
                   tags$li("簡化式 (Reduced Form) / ITT: Y ~ Z + X (或 rcs(X))", uiOutput("formula_reduced_form")),
                   tags$li("兩階段最小平方法 (2SLS / IV): Y ~ T + X (或 rcs(X)) | Z + X (或 rcs(X))", uiOutput("formula_2sls")),
                   tags$li("普通最小平方法 (OLS) / As-Treated: Y ~ T + X (或 rcs(X))", uiOutput("formula_ols"))
                 )
        ),
        tabPanel("數據預覽與遵從性檢查",
                 h4("模擬數據預覽 (前 6 筆)"),
                 tableOutput("data_head"),
                 hr(),
                 h4("第一階段迴歸 (遵從性檢查：Z 對 T 的影響)"),
                 p("此模型檢視隨機分派 (Z) 是否能有效預測實際接受的治療 (T)。"),
                 # UPDATED: Explanation of strong IV criteria
                 h5("工具變數強度 (Instrument Strength)"),
                 p("工具變數的『關聯性』假設要求 Z 必須與 T 相關。我們通常使用第一階段迴歸的 F 統計量來評估工具變數的強度。"),
                 tags$ul(
                    tags$li("一個常用的經驗法則是，如果排除內生迴歸變數後，檢驗工具變數聯合顯著性的 F 統計量大於 10 (Staiger & Stock, 1997; Stock & Yogo, 2005)，則認為工具變數足夠強。"),
                    tags$li("若 F 統計量小於 10，則可能存在『弱工具變數』(Weak Instrument) 問題。"),
                    tags$li("弱工具變數會導致：(1) 2SLS 估計量產生較大的有限樣本偏誤 (finite-sample bias)，可能偏向 OLS 估計量；(2) 常規的標準誤和假設檢定變得不可靠。"),
                    tags$li("因此，檢查第一階段 F 統計量是 IV 分析的重要步驟。")
                 ),
                 h5("線性基線 SAQ 模型結果:"),
                 verbatimTextOutput("first_stage_summary_linear"),
                 h5("RCS 基線 SAQ 模型結果:"),
                 verbatimTextOutput("first_stage_summary_rcs")
        ),
        tabPanel("模型結果比較",
                 h4("簡化式迴歸 (ITT 效應：Z 對 Y 的影響)"),
                 p("此模型估計了『意向治療』(ITT) 效應，即隨機『分派』對結果 Y 的平均影響。"),
                 h5("線性基線 SAQ 模型結果:"),
                 verbatimTextOutput("reduced_form_summary_linear"),
                 h5("RCS 基線 SAQ 模型結果:"),
                 verbatimTextOutput("reduced_form_summary_rcs"),
                 hr(),
                 h4("工具變數迴歸 (2SLS 結果：T 對 Y 的 LATE/CACE)"),
                 p("此模型使用 2SLS 估計了『局部平均處理效應』(LATE)。"),
                 h5("線性基線 SAQ 模型結果:"),
                 verbatimTextOutput("iv_summary_linear"),
                 h5("RCS 基線 SAQ 模型結果:"),
                 verbatimTextOutput("iv_summary_rcs"),
                 hr(),
                 h4("工具變數迴歸 (含穩健標準誤)"),
                 p("使用穩健標準誤可以使推論在異方差存在時更為可靠。"),
                 h5("線性基線 SAQ 模型結果 (穩健 SE):"),
                 verbatimTextOutput("iv_robust_summary_linear"),
                 h5("RCS 基線 SAQ 模型結果 (穩健 SE):"),
                 verbatimTextOutput("iv_robust_summary_rcs"),
                 hr(),
                 h4("普通最小平方法迴歸 (OLS / As-Treated)"),
                 p("此模型直接比較實際接受治療 (T=1) 與未接受治療 (T=0) 患者的結果 Y，可能存在偏誤。"),
                 h5("線性基線 SAQ 模型結果:"),
                 verbatimTextOutput("ols_summary_linear"),
                 h5("RCS 基線 SAQ 模型結果:"),
                 verbatimTextOutput("ols_summary_rcs"),
                 hr(),
                 # UPDATED: Comparison focuses on linear model for simplicity, but acknowledges RCS
                 h4("模型估計值比較 (基於線性 SAQ 模型)"),
                 uiOutput("comparison_estimates")
        ),
        tabPanel("完整模擬數據",
                 h4("完整的模擬數據集"),
                 DT::dataTableOutput("full_data_table")
        )
      )
    )
  )
)

# --- Server Logic (伺服器邏輯) ---
server <- function(input, output, session) {

  analysis_results <- eventReactive(input$run_analysis, {

    showNotification("正在模擬數據並執行分析...", type = "message", duration = 3)

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
    region <- factor(sample(1:3, n, replace = TRUE), labels = c("地區A", "地區B", "地區C"))

    # T: Treatment Received (Endogenous)
    prob_revasc <- ifelse(assignment_Z == 1, prob_revasc_if_Z1, prob_revasc_if_Z0)
    treatment_revasc_T <- rbinom(n, 1, prob_revasc)

    # Y: Outcome (Final SAQ Score) - Ground Truth (linear baseline effect)
    intercept <- 76
    beta_baseline <- 0.1 # True effect is linear for baseline SAQ
    beta_regionB <- 2
    beta_regionC <- -1
    error_sd <- 20

    outcome_Y <- intercept +
                 beta_T_true_late * treatment_revasc_T +
                 beta_baseline * baseline_SAQ + # Linear effect in generation
                 ifelse(region == "地區B", beta_regionB, 0) +
                 ifelse(region == "地區C", beta_regionC, 0) +
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

    # Create formulas
    formula_linear_cov <- "~ baseline_SAQ + region"
    formula_rcs_cov <- if (nk >= 3) paste("~ rcs(baseline_SAQ, nk=", nk, ") + region") else formula_linear_cov

    # Define formulas for each model type (Linear and RCS)
    # Linear Models
    formula_fs_linear <- as.formula(paste("T ~ Z +", formula_linear_cov))
    formula_rf_linear <- as.formula(paste("Y ~ Z +", formula_linear_cov))
    formula_ols_linear <- as.formula(paste("Y ~ T +", formula_linear_cov))
    formula_iv_linear <- as.formula(paste("Y ~ T +", formula_linear_cov, "| Z +", formula_linear_cov))

    # RCS Models (only if nk >= 3)
    formula_fs_rcs <- if (nk >= 3) as.formula(paste("T ~ Z +", formula_rcs_cov)) else NULL
    formula_rf_rcs <- if (nk >= 3) as.formula(paste("Y ~ Z +", formula_rcs_cov)) else NULL
    formula_ols_rcs <- if (nk >= 3) as.formula(paste("Y ~ T +", formula_rcs_cov)) else NULL
    formula_iv_rcs <- if (nk >= 3) as.formula(paste("Y ~ T +", formula_rcs_cov, "| Z +", formula_rcs_cov)) else NULL


    # --- Run Linear Models ---
    first_stage_linear_model <- lm(formula_fs_linear, data = ischemia_sim)
    first_stage_linear_summary <- summary(first_stage_linear_model)

    reduced_form_linear_model <- lm(formula_rf_linear, data = ischemia_sim)
    reduced_form_linear_summary <- summary(reduced_form_linear_model)

    ols_linear_model <- lm(formula_ols_linear, data = ischemia_sim)
    ols_linear_summary <- summary(ols_linear_model)

    iv_linear_model <- tryCatch({
        ivreg(formula_iv_linear, data = ischemia_sim)
      }, error = function(e) { message("IV Linear Model Error: ", e$message); NULL })
    iv_linear_summary <- if (!is.null(iv_linear_model)) summary(iv_linear_model) else "IV (Linear) 模型估計失敗"

    iv_linear_robust_summary <- NULL
    if (!is.null(iv_linear_model) && !is.character(iv_linear_model)) {
      iv_linear_robust_summary <- tryCatch({
          coeftest(iv_linear_model, vcov. = vcovHC(iv_linear_model, type = "HC1"))
        }, error = function(e) { message("Robust SE (Linear) Error: ", e$message); NULL })
    } else {
        iv_linear_robust_summary <- "無法計算穩健標準誤 (Linear)。"
    }

    # --- Run RCS Models (if nk >= 3) ---
    first_stage_rcs_model <- NULL
    first_stage_rcs_summary <- "未執行 RCS 模型 (結數 < 3)。"
    reduced_form_rcs_model <- NULL
    reduced_form_rcs_summary <- "未執行 RCS 模型 (結數 < 3)。"
    ols_rcs_model <- NULL
    ols_rcs_summary <- "未執行 RCS 模型 (結數 < 3)。"
    iv_rcs_model <- NULL
    iv_rcs_summary <- "未執行 RCS 模型 (結數 < 3)。"
    iv_rcs_robust_summary <- "未執行 RCS 模型 (結數 < 3)。"

    if (nk >= 3) {
        first_stage_rcs_model <- tryCatch({ lm(formula_fs_rcs, data = ischemia_sim) }, error = function(e){ message("FS RCS Error: ", e$message); NULL})
        first_stage_rcs_summary <- if (!is.null(first_stage_rcs_model)) summary(first_stage_rcs_model) else "第一階段 (RCS) 模型估計失敗"

        reduced_form_rcs_model <- tryCatch({ lm(formula_rf_rcs, data = ischemia_sim) }, error = function(e){ message("RF RCS Error: ", e$message); NULL})
        reduced_form_rcs_summary <- if (!is.null(reduced_form_rcs_model)) summary(reduced_form_rcs_model) else "簡化式 (RCS) 模型估計失敗"

        ols_rcs_model <- tryCatch({ lm(formula_ols_rcs, data = ischemia_sim) }, error = function(e){ message("OLS RCS Error: ", e$message); NULL})
        ols_rcs_summary <- if (!is.null(ols_rcs_model)) summary(ols_rcs_model) else "OLS (RCS) 模型估計失敗"

        iv_rcs_model <- tryCatch({
            ivreg(formula_iv_rcs, data = ischemia_sim)
        }, error = function(e) { message("IV RCS Model Error: ", e$message); NULL })
        iv_rcs_summary <- if (!is.null(iv_rcs_model)) summary(iv_rcs_model) else "IV (RCS) 模型估計失敗"

        if (!is.null(iv_rcs_model) && !is.character(iv_rcs_model)) {
          iv_rcs_robust_summary <- tryCatch({
              coeftest(iv_rcs_model, vcov. = vcovHC(iv_rcs_model, type = "HC1"))
            }, error = function(e) { message("Robust SE (RCS) Error: ", e$message); NULL })
           if(is.null(iv_rcs_robust_summary)) iv_rcs_robust_summary <- "無法計算穩健標準誤 (RCS)。"
        } else {
            iv_rcs_robust_summary <- "無法計算穩健標準誤 (RCS)。"
        }
    }


    # Store ground truth parameters
    ground_truth_params <- list(
        n = n,
        intercept = intercept,
        beta_T_true_late = beta_T_true_late,
        beta_baseline = beta_baseline, # True effect is linear
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
      first_stage_linear_model = first_stage_linear_model,
      first_stage_linear_summary = first_stage_linear_summary,
      reduced_form_linear_model = reduced_form_linear_model,
      reduced_form_linear_summary = reduced_form_linear_summary,
      ols_linear_model = ols_linear_model,
      ols_linear_summary = ols_linear_summary,
      iv_linear_model = iv_linear_model,
      iv_linear_summary = iv_linear_summary,
      iv_linear_robust_summary = iv_linear_robust_summary,
      # RCS models
      first_stage_rcs_summary = first_stage_rcs_summary,
      reduced_form_rcs_summary = reduced_form_rcs_summary,
      ols_rcs_summary = ols_rcs_summary,
      iv_rcs_summary = iv_rcs_summary,
      iv_rcs_robust_summary = iv_rcs_robust_summary,
      # Params
      ground_truth_params = ground_truth_params,
      nk = nk # Pass nk to UI
    )
  }, ignoreNULL = FALSE)

  # --- 3. Render Outputs ---

  # Render data table preview
  output$data_head <- renderTable({
    results <- analysis_results()
    req(results$data_full)
    head(results$data_full)
  }, rownames = TRUE)

  # Render the FULL data table using DT
  output$full_data_table <- DT::renderDataTable({
    results <- analysis_results()
    req(results$data_full)
    DT::datatable(results$data_full,
                  options = list(pageLength = 10, scrollX = TRUE, searching = TRUE),
                  rownames = FALSE,
                  filter = 'top')
  })

  # --- Render Model Summaries (Linear and RCS) ---

  # Helper function to print summary and F-stat (if applicable)
  print_summary_and_fstat <- function(summary_obj, model_type = "Linear") {
      if (is.character(summary_obj)) { # Handle error messages
          cat(paste(model_type, ":", summary_obj, "\n"))
          return()
      }
      req(summary_obj)
      print(summary_obj)

      # Extract and display F-statistic for first stage models
      if (inherits(summary_obj, "summary.lm") && !is.null(summary_obj$fstatistic)) {
          f_stat <- summary_obj$fstatistic
          f_value <- f_stat[1]
          cat("\n---\n")
          cat(paste("第一階段 F 統計量 (檢定工具 Z 的聯合顯著性):", round(f_value, 2), "\n"))
          if (f_value < 10) {
              cat("警告：第一階段 F 統計量 < 10。可能存在弱工具變數問題。\n")
          } else {
              cat("第一階段 F 統計量 >= 10，工具變數強度尚可接受。\n")
          }
      }
  }

  # First Stage
  output$first_stage_summary_linear <- renderPrint({
      results <- analysis_results()
      req(results$first_stage_linear_summary)
      print_summary_and_fstat(results$first_stage_linear_summary, "Linear")
  })
  output$first_stage_summary_rcs <- renderPrint({
      results <- analysis_results()
      req(results$first_stage_rcs_summary)
      # Check if RCS was actually run
      if(is.character(results$first_stage_rcs_summary) || results$nk < 3) {
          cat(results$first_stage_rcs_summary)
      } else {
          print_summary_and_fstat(results$first_stage_rcs_summary, "RCS")
      }
  })

  # Reduced Form
  output$reduced_form_summary_linear <- renderPrint({
      results <- analysis_results()
      req(results$reduced_form_linear_summary)
      print(results$reduced_form_linear_summary)
  })
  output$reduced_form_summary_rcs <- renderPrint({
      results <- analysis_results()
      req(results$reduced_form_rcs_summary)
      if(is.character(results$reduced_form_rcs_summary) || results$nk < 3) {
          cat(results$reduced_form_rcs_summary)
      } else {
          print(results$reduced_form_rcs_summary)
      }
  })

  # IV / 2SLS
  output$iv_summary_linear <- renderPrint({
    results <- analysis_results()
    req(results$iv_linear_summary)
    if (is.character(results$iv_linear_summary)) {
        cat(results$iv_linear_summary)
    } else {
        print(results$iv_linear_summary)
    }
  })
  output$iv_summary_rcs <- renderPrint({
    results <- analysis_results()
    req(results$iv_rcs_summary)
    if (is.character(results$iv_rcs_summary) || results$nk < 3) {
        cat(results$iv_rcs_summary)
    } else {
        print(results$iv_rcs_summary)
    }
  })

  # IV / 2SLS Robust
  output$iv_robust_summary_linear <- renderPrint({
     results <- analysis_results()
     req(results$iv_linear_robust_summary)
     if (is.character(results$iv_linear_robust_summary)) {
         cat(results$iv_linear_robust_summary)
     } else {
         print(results$iv_linear_robust_summary)
     }
  })
   output$iv_robust_summary_rcs <- renderPrint({
     results <- analysis_results()
     req(results$iv_rcs_robust_summary)
     if (is.character(results$iv_rcs_robust_summary) || results$nk < 3) {
         cat(results$iv_rcs_robust_summary)
     } else {
         print(results$iv_rcs_robust_summary)
     }
  })

  # OLS
  output$ols_summary_linear <- renderPrint({
     results <- analysis_results()
     req(results$ols_linear_summary)
     print(results$ols_linear_summary)
  })
  output$ols_summary_rcs <- renderPrint({
     results <- analysis_results()
     req(results$ols_rcs_summary)
     if (is.character(results$ols_rcs_summary) || results$nk < 3) {
         cat(results$ols_rcs_summary)
     } else {
         print(results$ols_rcs_summary)
     }
  })

  # --- Render Formulas ---
  # (Formulas remain the same structurally, interpretation depends on model choice)
  output$ground_truth_formula_display <- renderUI({
      withMathJax(helpText(
        "$$ Y = \\beta_0 + \\beta_{T, LATE} T + \\beta_{baseline} \\text{baseline\\_SAQ} + ... + \\epsilon $$"
      ))
  })
  output$formula_first_stage <- renderUI({
      withMathJax(helpText(
        "$$ T = \\gamma_0 + \\gamma_Z Z + f(\\text{Covariates}) + \\nu $$"
      ))
  })
  output$formula_reduced_form <- renderUI({
      withMathJax(helpText(
        "$$ Y = \\pi_0 + \\pi_Z Z + f(\\text{Covariates}) + \\omega $$"
      ))
  })
  output$formula_2sls <- renderUI({
      withMathJax(helpText(
        "$$ Y = \\beta_0 + \\beta_{T, 2SLS} \\hat{T} + f(\\text{Covariates}) + \\epsilon $$"
      ))
  })
  output$formula_ols <- renderUI({
      withMathJax(helpText(
        "$$ Y = \\delta_0 + \\delta_T T + f(\\text{Covariates}) + \\mu $$"
      ))
  })

  # Render ground truth parameters
  output$ground_truth_params_display <- renderPrint({
      results <- analysis_results()
      req(results$ground_truth_params)
      params <- results$ground_truth_params
      prop_compliers = params$prob_revasc_if_Z1 - params$prob_revasc_if_Z0
      cat(paste("樣本數 (n):", params$n, "\n"))
      cat(paste("真實 LATE/CACE (β_T):", params$beta_T_true_late, "\n"))
      cat(paste("真實基線 SAQ 效應 (β_baseline, 線性):", params$beta_baseline, "\n")) # Clarify linear
      # ... other params ...
      cat(paste("P(T=1 | Z=1):", params$prob_revasc_if_Z1, "\n"))
      cat(paste("P(T=1 | Z=0):", params$prob_revasc_if_Z0, "\n"))
      cat(paste("遵從者比例估計:", round(prop_compliers, 3), "\n"))
  })

  # UPDATED: Comparison uses linear models, removed <strong> tags
  output$comparison_estimates <- renderUI({
      results <- analysis_results()
      req(results$reduced_form_linear_model, results$iv_linear_model, results$ols_linear_model, results$ground_truth_params)

      # Extract coefficients from LINEAR models safely
      itt_coef <- tryCatch(coef(results$reduced_form_linear_model)["Z"], error = function(e) NA)
      iv_coef <- tryCatch(coef(results$iv_linear_model)["T"], error = function(e) NA)
      ols_coef <- tryCatch(coef(results$ols_linear_model)["T"], error = function(e) NA)
      true_late <- results$ground_truth_params$beta_T_true_late

      if (is.na(itt_coef) || is.na(iv_coef) || is.na(ols_coef)) {
          return(p("無法提取所有線性模型係數進行比較。"))
      }

      # Create comparison text (without <strong>)
      tagList(
          p(paste0("在此模擬中，真實的遵從者平均因果效應 (LATE/CACE) 設定為: ", round(true_late, 3))),
          hr(),
          h5("不同分析方法的估計值 (基於線性 SAQ 模型)："),
          tags$ul(
              tags$li(paste0("ITT (意向治療) 效應: ", round(itt_coef, 3), ". 隨機『分派』對結果 Y 的平均影響。")),
              tags$li(paste0("OLS (依治療分析) 關聯性: ", round(ols_coef, 3), ". 比較實際『接受』治療者的結果差異，可能有偏誤。")),
              tags$li(paste0("2SLS (工具變數) 效應 (LATE/CACE): ", round(iv_coef, 3), ". 利用隨機分派 Z 估計實際『接受』治療 T 對『遵從者』的因果效應。"))
          ),
          hr(),
          h5("比較與解釋："),
          tags$ul(
              tags$li(paste0("2SLS 估計值 (", round(iv_coef, 3), ") 通常最接近真實 LATE (", round(true_late, 3), ")，因其校正選擇偏誤。")),
              tags$li(paste0("OLS 估計值 (", round(ols_coef, 3), ") 與真實 LATE 的差異反映了『遵從行為』與結果之間的混淆。")),
              tags$li(paste0("ITT 估計值 (", round(itt_coef, 3), ") 代表治療『策略』的平均效果，受遵從者比例影響 (ITT ≈ LATE × 遵從者比例)。"))
          ),
          hr(),
          p("注意：以上比較基於假設基線 SAQ 效應為線性的模型。您可以查看其他標籤頁中包含 RCS 的模型結果，以評估非線性假設是否影響估計。")
      )
  })

  # --- Download Handler ---
  output$download_data <- downloadHandler(
    filename = function() {
      paste0("simulated_ischemia_data_", Sys.Date(), "_nk", analysis_results()$nk, ".csv") # Include nk in filename
    },
    content = function(file) {
      results <- analysis_results()
      if (!is.null(results$data_full)) {
        write.csv(results$data_full, file, row.names = FALSE, fileEncoding = "UTF-8")
      } else {
        write.csv(data.frame(Error = "No simulation data generated yet."), file, row.names = FALSE)
      }
    }
  )

}

# --- Run the application ---
shinyApp(ui = ui, server = server)
