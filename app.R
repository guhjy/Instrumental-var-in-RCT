# --- Load Necessary Libraries ---
# shiny: Web application framework for R
# AER: Provides functions for applied econometrics, including ivreg() for IV estimation.
# sandwich: Provides functions for robust covariance matrix estimators.
# lmtest: Provides functions for hypothesis testing, including coeftest() for displaying coefficients with robust SEs.
# DT: Provides an interface to the DataTables library for interactive tables.

# --- References ---
# General IV Concepts in RCTs:
# Hernán MA, Robins JM. Instruments for Causal Inference: An Epidemiologist’s Dream? Epidemiology. 2006;17(4):360-372.
# Baiocchi M, Cheng J, Small DS. Instrumental Variable Methods for Causal Inference. Stat Med. 2014;33(13):2297-2340.

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

library(shiny)
library(AER)
library(sandwich)
library(lmtest)
library(DT) # Load DT for interactive tables

# --- UI Definition (使用者介面) ---
ui <- fluidPage(
  # Set language for the page
  tags$head(
    tags$meta(charset = "UTF-8"),
    tags$html(lang = "zh-Hant-TW"),
    # Include MathJax for rendering formulas
    tags$script(src = "https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.7/MathJax.js?config=TeX-MML-AM_CHTML")
  ),

  # App title
  titlePanel("工具變數分析 Shiny 應用程式 (模擬 RCT 中的遵從性問題)"),

  # Sidebar layout with input and output definitions
  sidebarLayout(
    # Sidebar panel for inputs
    sidebarPanel(
      h4("控制面板"),

      # --- Simulation Parameter Inputs ---
      h5("模擬參數設定:"),
      numericInput("n_sim", "樣本數 (Sample Size, n):", value = 5000, min = 100, max = 20000, step = 100),
      numericInput("beta_T_sim", "真實治療效應 (Complier Average Causal Effect, CACE/LATE):", value = 5.5, min = -10, max = 20, step = 0.5),
      sliderInput("compliance_inv", "遵從率 - 分派到介入組 (Proportion Receiving Treatment if Assigned Intervention):", min = 0, max = 1, value = 0.80, step = 0.01),
      sliderInput("compliance_con", "非遵從率 - 分派到對照組 (Proportion Receiving Treatment if Assigned Control):", min = 0, max = 1, value = 0.12, step = 0.01),
      hr(), # Horizontal line

      # Action button to trigger simulation and analysis
      actionButton("run_analysis", "重新模擬數據並執行分析", icon = icon("play"), class = "btn-primary"), # Added btn-primary for emphasis
      hr(), # Horizontal line

      # --- Download Button ---
      downloadButton("download_data", "下載模擬數據 (.csv)", class = "btn-success"),
      hr(), # Horizontal line


      h5("說明:"),
      p("此應用程式演示了在隨機對照試驗 (RCT) 中，當存在治療方案不遵從 (non-compliance) 時，如何使用工具變數 (Instrumental Variables, IV) 方法來估計治療的因果效應。"),
      p("模擬數據旨在反映類似 ISCHEMIA 試驗 (Maron et al., NEJM 2020) 的情況，其中隨機分派 (Z) 作為工具變數。"),
      p("調整上面的模擬參數，然後點擊按鈕以生成新的模擬數據集，並在右側查看分析結果。"),

    ),

    # Main panel for displaying outputs
    mainPanel(
      # Output: Tabset w/ different analysis results
      tabsetPanel(
        id = "results_tabs",
        tabPanel("分析說明與模擬設定",
                 h3("工具變數 (Instrumental Variables, IV) 與隨機試驗"),
                 p("在理想的隨機對照試驗 (RCT) 中，隨機分派 (Z) 保證了治療組間的基線可比性。然而，患者可能不完全遵從分配給他們的治療方案。例如，分配到介入治療 (Z=1) 的患者可能未接受治療，而分配到對照組 (Z=0) 的患者可能自行尋求了介入治療。這種情況稱為「不遵從性」(non-compliance)。"),
                 p("當存在不遵從性時，直接比較實際接受治療 (T=1) 和未接受治療 (T=0) 的患者（即『依治療分析』, as-treated analysis using OLS）可能會因為『選擇偏誤』而產生誤導，因為決定是否接受治療的因素可能也與結果 (Y) 相關。"),
                 p("工具變數 (IV) 分析利用最初的「隨機分派」(Z) 作為「工具」，來估計「實際接受治療」(T) 對「結果」(Y) 的因果效應。在 RCT 的背景下，IV 估計的是『局部平均處理效應』(Local Average Treatment Effect, LATE)，也稱為『遵從者平均因果效應』(Complier Average Causal Effect, CACE)。"),
                 p("LATE/CACE 是指治療對那些其治療狀態因隨機分派而改變的『遵從者』(Compliers) 的平均效應。遵從者是指如果被分派到介入組就會接受治療，如果被分派到對照組就不會接受治療的個體。"),
                 hr(),
                 h4("工具變數分析的核心假設:"),
                 tags$ul(
                   tags$li(strong("1. 關聯性 (Relevance):"), "隨機分派 (Z) 必須與實際接受的治療 (T) 相關。"),
                   tags$li(strong("2. 獨立性/隨機性 (Independence/Randomization):"), "隨機分派 (Z) 與任何影響結果 Y 的未觀察混淆因素無關 (由 RCT 設計保證)。"),
                   tags$li(strong("3. 排除限制 (Exclusion Restriction):"), "隨機分派 (Z) 只能透過影響實際接受的治療 (T) 來影響結果 (Y)，不能有其他直接影響結果的路徑。"),
                   tags$li(strong("4. 單調性 (Monotonicity):"), "不存在『反抗者』(Defiers)。")
                 ),
                 hr(),
                 h4("在此模擬中使用的變數："),
                 tags$ul(
                   # UPDATED: Explanation for Y (SAQ)
                   tags$li(strong("Y:"), " 結果變數。模擬一年後的西雅圖心絞痛問卷 (Seattle Angina Questionnaire, SAQ) 生活品質分數。範圍 0-100，分數越高代表健康狀況越好。"),
                   tags$li(strong("T:"), " 內生變數 (實際接受的血管重建治療，1=是, 0=否) - 受不遵從性影響"),
                   tags$li(strong("Z:"), " 工具變數 (隨機分派的治療策略，1=介入性, 0=保守性) - 不受患者行為影響"),
                   # UPDATED: Changed baseline_saq to baseline_SAQ
                   tags$li(strong("X:"), " 控制變數/共變數 (模擬基線 SAQ 分數 `baseline_SAQ`、地區 `region`)")
                 ),
                 hr(),
                 h4("數據生成模型 (Ground Truth Formula for Outcome Y)"),
                 p("模擬數據中的結果變數 (Y) 是根據以下公式生成的（注意：此公式描述了如果治療 T 被外生決定時的關係）："),
                 uiOutput("ground_truth_formula_display"), # Use uiOutput to render HTML formula
                 p("其中 ε 是隨機誤差項，服從常態分佈 N(0, σ²)。β_T 代表接受治療相對於未接受治療的真實因果效應（對於遵從者而言）。"),
                 h5("目前模擬使用的參數值:"),
                 verbatimTextOutput("ground_truth_params_display"), # Display parameters used
                 hr(),
                 h4("主要分析方法與公式："),
                 tags$ol(
                   tags$li(
                       strong("第一階段 (First Stage):"), " 迴歸實際接受的治療 (T) 對隨機分派 (Z) 和共變數 (X)。檢驗工具變數的『關聯性』假設（Z 的係數是否顯著不為零）和強度（F 統計量是否 > 10）。",
                       uiOutput("formula_first_stage")
                   ),
                   tags$li(
                       strong("簡化式 (Reduced Form) / 意向治療 (Intention-to-Treat, ITT):"), " 迴歸結果 (Y) 對隨機分派 (Z) 和共變數 (X)。Z 的係數估計了『分派』到介入組相對於對照組對結果的平均影響 (ITT 效應)。保留隨機分派優點，但評估的是『分派策略』而非『實際治療』效應。",
                       uiOutput("formula_reduced_form")
                   ),
                   tags$li(
                       strong("兩階段最小平方法 (2SLS / IV):"), " 結合第一階段和簡化式的資訊，估計實際治療 (T) 對結果 (Y) 的『局部平均處理效應』(LATE/CACE)。這是對『遵從者』的治療效應估計。LATE ≈ ITT 效應 / 第一階段效應。",
                       uiOutput("formula_2sls")
                   ),
                   tags$li(
                       strong("普通最小平方法 (OLS) / 依治療分析 (As-Treated):"), " 直接迴歸結果 (Y) 對實際接受的治療 (T) 和共變數 (X)。當存在不遵從性時，此分析通常因忽略 T 的內生性（選擇偏誤）而產生偏誤，其結果解釋為關聯性而非因果效應。",
                       uiOutput("formula_ols")
                   )
                 )
        ),
        tabPanel("數據預覽與遵從性檢查",
                 h4("模擬數據預覽 (前 6 筆)"),
                 tableOutput("data_head"),
                 hr(),
                 h4("第一階段迴歸 (遵從性檢查：Z 對 T 的影響)"),
                 p("此模型檢視隨機分派 (Z) 是否能有效預測實際接受的治療 (T)。Z 的係數估計了隨機分派到介入組 (Z=1) 相對於對照組 (Z=0)，使患者實際接受治療 (T=1) 的機率平均增加了多少。這個差異估計了樣本中『遵從者』的比例。"),
                 p("F 統計量檢驗了工具變數 Z 的顯著性。若 F 統計量遠大於 10，則表示工具變數的關聯性較強。"),
                 verbatimTextOutput("first_stage_summary")
        ),
        tabPanel("模型結果比較",
                 h4("簡化式迴歸 (ITT 效應：Z 對 Y 的影響)"),
                 p("此模型估計了『意向治療』(ITT) 效應。Z 的係數代表隨機『分派』到介入組相對於對照組，對最終結果 Y 的平均影響。這是對治療『策略』效應的無偏估計。"),
                 verbatimTextOutput("reduced_form_summary"),
                 hr(),
                 h4("工具變數迴歸 (2SLS 結果：T 對 Y 的因果效應估計 - LATE/CACE)"),
                 p("此模型使用 2SLS 估計了『局部平均處理效應』(LATE)，即實際接受治療 (T=1) 相對於未接受 (T=0) 對結果 Y 的平均因果效應，但此效應僅適用於『遵從者』群體。T 的係數即為 LATE 估計值。"),
                 p("理論上，LATE ≈ (簡化式中 Z 的係數) / (第一階段中 Z 的係數)。"),
                 verbatimTextOutput("iv_summary"),
                 hr(),
                 h4("工具變數迴歸 (含穩健標準誤)"),
                 p("使用穩健標準誤 (Robust Standard Errors) 可以使推論在異方差 (heteroskedasticity) 存在時更為可靠。"),
                 verbatimTextOutput("iv_robust_summary"),
                 hr(),
                 h4("普通最小平方法迴歸 (OLS / As-Treated：T 對 Y 的關聯性 - 可能有偏誤)"),
                 p("此模型直接比較實際接受治療 (T=1) 與未接受治療 (T=0) 患者的結果 Y，控制了共變數 X。由於 T 的選擇可能與影響 Y 的未觀察因素有關（內生性），T 的係數估計的是關聯性，通常不是無偏的因果效應估計。"),
                 verbatimTextOutput("ols_summary"),
                 hr(),
                 # UPDATED: Enhanced comparison section
                 h4("模型估計值比較 (ITT vs. OLS vs. 2SLS)"),
                 uiOutput("comparison_estimates") # Use uiOutput for dynamic text
        ),
        tabPanel("完整模擬數據",
                 h4("完整的模擬數據集"),
                 p("下方表格顯示了本次模擬生成的完整數據。您可以使用搜索框進行篩選，並點擊列標題進行排序。"),
                 DT::dataTableOutput("full_data_table")
        )
      )
    )
  )
)

# --- Server Logic (伺服器邏輯) ---
server <- function(input, output, session) {

  # Reactive expression for simulation and analysis, triggered by button
  analysis_results <- eventReactive(input$run_analysis, {

    showNotification("正在模擬數據並執行分析...", type = "message", duration = 3)

    # --- 1. Simulate Data ---
    set.seed(as.integer(Sys.time())) # Use time-based seed for variability

    n <- input$n_sim
    beta_T_true_late <- input$beta_T_sim # True LATE/CACE
    prob_revasc_if_Z1 <- input$compliance_inv # P(T=1|Z=1)
    prob_revasc_if_Z0 <- input$compliance_con # P(T=1|Z=0)

    # Z: Random assignment (Instrument)
    assignment_Z <- rbinom(n, 1, 0.5) # 1 = Intervention, 0 = Control

    # X: Covariates
    # UPDATED: Renamed baseline_saq to baseline_SAQ
    baseline_SAQ <- rnorm(n, mean = 60, sd = 15) # Baseline SAQ score
    region <- factor(sample(1:3, n, replace = TRUE), labels = c("地區A", "地區B", "地區C"))

    # T: Treatment Received (Endogenous Variable) - based on Z and compliance
    prob_revasc <- ifelse(assignment_Z == 1, prob_revasc_if_Z1, prob_revasc_if_Z0)
    treatment_revasc_T <- rbinom(n, 1, prob_revasc) # 1 = Received Revascularization, 0 = Did not

    # Y: Outcome (Final SAQ Score) - Ground Truth Formula
    # Y represents the final SAQ score (0-100, higher is better)
    intercept <- 76
    beta_baseline <- 0.1 # Effect of baseline SAQ
    beta_regionB <- 2
    beta_regionC <- -1
    error_sd <- 20

    # Calculate outcome Y
    # The effect beta_T_true_late applies *when* T=1 occurs due to assignment (implicitly for compliers)
    # UPDATED: Use baseline_SAQ in the formula
    outcome_Y <- intercept +
                 beta_T_true_late * treatment_revasc_T +
                 beta_baseline * baseline_SAQ +
                 ifelse(region == "地區B", beta_regionB, 0) +
                 ifelse(region == "地區C", beta_regionC, 0) +
                 rnorm(n, mean = 0, sd = error_sd) # Random error term (ε)

    # Create data frame
    ischemia_sim <- data.frame(
      Y = outcome_Y,
      T_revasc = treatment_revasc_T, # Actual treatment
      Z_assign = assignment_Z,      # Assigned treatment
      # UPDATED: Use baseline_SAQ
      baseline_SAQ = baseline_SAQ,
      region = region
    )
    # Rename columns for model formulas (use simple T, Z, and keep baseline_SAQ)
    # UPDATED: Column names reflect changes
    colnames(ischemia_sim) <- c("Y", "T", "Z", "baseline_SAQ", "region")


    # --- 2. Perform Analyses ---
    # UPDATED: Formulas use baseline_SAQ
    # First Stage Regression (T ~ Z + Covariates)
    first_stage_model <- lm(T ~ Z + baseline_SAQ + region, data = ischemia_sim)
    first_stage_summary_obj <- summary(first_stage_model)

    # Reduced Form Regression (Y ~ Z + Covariates) - ITT effect
    reduced_form_model <- lm(Y ~ Z + baseline_SAQ + region, data = ischemia_sim)
    reduced_form_summary_obj <- summary(reduced_form_model)

    # IV Regression (2SLS: Y ~ T + Covariates | Z + Covariates) - LATE
    iv_model <- tryCatch({
        ivreg(Y ~ T + baseline_SAQ + region | Z + baseline_SAQ + region, data = ischemia_sim)
      }, error = function(e) {
        showNotification(paste("IV 估計錯誤:", e$message), type = "error", duration = 5)
        return(NULL)
      })
    iv_summary_obj <- if (!is.null(iv_model)) summary(iv_model) else "IV 模型估計失敗"

    # Naive OLS Regression (Y ~ T + Covariates) - As-treated analysis
    naive_ols_model <- lm(Y ~ T + baseline_SAQ + region, data = ischemia_sim)
    ols_summary_obj <- summary(naive_ols_model)

    # Calculate Robust Standard Errors for IV model
    iv_robust_summary_output <- NULL
    if (!is.null(iv_model) && !is.character(iv_model)) {
      iv_robust_summary_output <- tryCatch({
          coeftest(iv_model, vcov. = vcovHC(iv_model, type = "HC1"))
        }, error = function(e) {
          showNotification(paste("計算穩健標準誤時發生錯誤:", e$message), type = "warning", duration = 5)
          return(NULL)
      })
    } else {
        iv_robust_summary_output <- "IV 模型不存在或估計失敗，無法計算穩健標準誤。"
    }

    # Store ground truth parameters
    ground_truth_params <- list(
        n = n,
        intercept = intercept,
        beta_T_true_late = beta_T_true_late, # True LATE/CACE
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
      first_stage_model = first_stage_model, # Return models for coef extraction
      first_stage_summary = first_stage_summary_obj,
      reduced_form_model = reduced_form_model, # Return models for coef extraction
      reduced_form_summary = reduced_form_summary_obj,
      iv_model = iv_model, # Return models for coef extraction
      iv_summary = iv_summary_obj,
      iv_robust_summary = iv_robust_summary_output,
      ols_model = naive_ols_model, # Return models for coef extraction
      ols_summary = ols_summary_obj,
      ground_truth_params = ground_truth_params
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

  # Render first stage model summary
  output$first_stage_summary <- renderPrint({
    results <- analysis_results()
    req(results$first_stage_summary, !is.character(results$first_stage_summary))
    print(results$first_stage_summary)
    f_stat <- results$first_stage_summary$fstatistic
    if(!is.null(f_stat) && length(f_stat) >= 1){
      f_value <- f_stat[1]
      cat("\n---\n")
      cat(paste("第一階段 F 統計量 (檢定工具 Z 的聯合顯著性):", round(f_value, 2), "\n"))
      if(f_value < 10){
         cat("警告：第一階段 F 統計量 < 10。可能存在弱工具變數問題。\n")
      } else {
         cat("第一階段 F 統計量 >= 10，工具變數強度尚可接受。\n")
      }
    } else {
        cat("\n---\n無法提取第一階段 F 統計量。\n")
    }
  })

  # Render reduced form model summary (ITT)
  output$reduced_form_summary <- renderPrint({
    results <- analysis_results()
    req(results$reduced_form_summary, !is.character(results$reduced_form_summary))
    print(results$reduced_form_summary)
  })

  # Render IV model summary (LATE)
  output$iv_summary <- renderPrint({
    results <- analysis_results()
    req(results$iv_summary)
    if (is.character(results$iv_summary)) {
        cat(results$iv_summary)
    } else {
        print(results$iv_summary)
    }
  })

  # Render IV model summary with robust SEs
  output$iv_robust_summary <- renderPrint({
     results <- analysis_results()
     req(results$iv_robust_summary)
     if (is.character(results$iv_robust_summary)) {
         cat(results$iv_robust_summary)
     } else {
         print(results$iv_robust_summary)
     }
  })

  # Render OLS model summary (As-Treated)
  output$ols_summary <- renderPrint({
     results <- analysis_results()
     req(results$ols_summary, !is.character(results$ols_summary))
     print(results$ols_summary)
  })

  # --- Render Formulas using MathJax ---
  output$ground_truth_formula_display <- renderUI({
      # UPDATED: Use baseline_SAQ in formula display
      withMathJax(helpText(
        "$$ Y = \\beta_0 + \\beta_{T, LATE} T + \\beta_{baseline} \\text{baseline\\_SAQ} + \\beta_{regionB} I(\\text{region=B}) + \\beta_{regionC} I(\\text{region=C}) + \\epsilon $$"
      ))
  })

  output$formula_first_stage <- renderUI({
      # UPDATED: Use baseline_SAQ in formula display
      withMathJax(helpText(
        "$$ T = \\gamma_0 + \\gamma_Z Z + \\gamma_{baseline} \\text{baseline\\_SAQ} + \\gamma_{region} \\text{region} + \\nu $$"
      ))
  })

  output$formula_reduced_form <- renderUI({
      # UPDATED: Use baseline_SAQ in formula display
      withMathJax(helpText(
        "$$ Y = \\pi_0 + \\pi_Z Z + \\pi_{baseline} \\text{baseline\\_SAQ} + \\pi_{region} \\text{region} + \\omega $$"
      ))
  })

  output$formula_2sls <- renderUI({
      # UPDATED: Use baseline_SAQ in formula display
      withMathJax(helpText(
        "概念上，第一階段預測 \\( \\hat{T} \\)。第二階段使用 \\( \\hat{T} \\) 取代 T:",
        "$$ Y = \\beta_0 + \\beta_{T, 2SLS} \\hat{T} + \\beta_{baseline} \\text{baseline\\_SAQ} + \\beta_{region} \\text{region} + \\epsilon $$",
        "其中 \\( \\beta_{T, 2SLS} \\approx \\frac{\\hat{\\pi}_Z}{\\hat{\\gamma}_Z} \\) (LATE 估計值)。"
      ))
  })

  output$formula_ols <- renderUI({
      # UPDATED: Use baseline_SAQ in formula display
      withMathJax(helpText(
        "$$ Y = \\delta_0 + \\delta_T T + \\delta_{baseline} \\text{baseline\\_SAQ} + \\delta_{region} \\text{region} + \\mu $$"
      ))
  })


  # Render the ground truth parameters used in the simulation
  output$ground_truth_params_display <- renderPrint({
      results <- analysis_results()
      req(results$ground_truth_params)
      params <- results$ground_truth_params
      prop_compliers = params$prob_revasc_if_Z1 - params$prob_revasc_if_Z0

      cat(paste("樣本數 (n):", params$n, "\n"))
      cat(paste("截距項 (β₀):", params$intercept, "\n"))
      cat(paste("真實 LATE/CACE (β_T):", params$beta_T_true_late, "\n"))
      # UPDATED: Reflect baseline_SAQ name change
      cat(paste("基線 SAQ 效應 (β_baseline):", params$beta_baseline, "\n"))
      cat(paste("地區 B 效應 (β_regionB):", params$beta_regionB, "\n"))
      cat(paste("地區 C 效應 (β_regionC):", params$beta_regionC, "\n"))
      cat(paste("誤差標準差 (σ_ε):", params$error_sd, "\n"))
      cat(paste("P(T=1 | Z=1) (介入組治療率):", params$prob_revasc_if_Z1, "\n"))
      cat(paste("P(T=1 | Z=0) (對照組治療率):", params$prob_revasc_if_Z0, "\n"))
      cat(paste("遵從者比例估計 (P(T=1|Z=1) - P(T=1|Z=0)):", round(prop_compliers, 3), "\n"))
  })

  # UPDATED: Enhanced comparison of ITT, OLS, 2SLS estimates
  output$comparison_estimates <- renderUI({
      results <- analysis_results()
      # Ensure all required models and parameters are available
      req(results$reduced_form_model, results$iv_model, results$ols_model, results$ground_truth_params)

      # Extract coefficients safely
      itt_coef <- tryCatch(coef(results$reduced_form_model)["Z"], error = function(e) NA)
      iv_coef <- tryCatch(coef(results$iv_model)["T"], error = function(e) NA)
      ols_coef <- tryCatch(coef(results$ols_model)["T"], error = function(e) NA)
      true_late <- results$ground_truth_params$beta_T_true_late

      if (is.na(itt_coef) || is.na(iv_coef) || is.na(ols_coef)) {
          return(p("無法提取所有模型係數進行比較。"))
      }

      # Create comparison text
      tagList(
          p(paste0("在此模擬中，真實的遵從者平均因果效應 (LATE/CACE) 設定為: ", strong(round(true_late, 3)))),
          hr(),
          h5("不同分析方法的估計值："),
          tags$ul(
              tags$li(paste0(strong("ITT (意向治療) 效應: "), round(itt_coef, 3), "。這是隨機『分派』到介入組 (Z=1) 相對於對照組 (Z=0) 對結果 Y 的平均影響。它保留了隨機化的優點，估計了『策略』的效果，但因不遵從性而被稀釋。")),
              tags$li(paste0(strong("OLS (依治療分析) 關聯性: "), round(ols_coef, 3), "。這是比較實際『接受』治療 (T=1) 與未接受治療 (T=0) 患者的平均結果差異。由於 T 的選擇不是隨機的（存在選擇偏誤/內生性），此估計值通常是有偏誤的因果效應估計。")),
              tags$li(paste0(strong("2SLS (工具變數) 效應 (LATE/CACE): "), round(iv_coef, 3), "。這是利用隨機分派 Z 作為工具變數，估計實際『接受』治療 T 對結果 Y 的因果效應，但此效應僅適用於『遵從者』。理論上，LATE ≈ ITT / (第一階段 Z 的係數)。"))
          ),
          hr(),
          h5("比較與解釋："),
          tags$ul(
              tags$li(paste0("2SLS 估計值 (", round(iv_coef, 3), ") 通常最接近真實的 LATE (", round(true_late, 3), ")，因為它校正了 OLS 中的選擇偏誤，並針對『遵從者』估計了治療本身的因果效應。")),
              tags$li(paste0("OLS 估計值 (", round(ols_coef, 3), ") 與真實 LATE 的差異 (偏誤) 反映了『遵從行為』與結果之間的混淆。例如，如果更傾向於接受治療的患者本身預後就更好（即使不治療），OLS 可能會高估治療效果。")),
              tags$li(paste0("ITT 估計值 (", round(itt_coef, 3), ") 通常介於 0 和 LATE 之間（如果 LATE 為正）。它代表了在現實世界中實施該治療『策略』的平均效果，考慮到了預期的不遵從性。其大小受遵從者比例的影響 (ITT ≈ LATE × 遵從者比例)。"))
          )
      )
  })

  # --- Download Handler ---
  output$download_data <- downloadHandler(
    filename = function() {
      paste0("simulated_ischemia_data_", Sys.Date(), ".csv")
    },
    content = function(file) {
      results <- analysis_results()
      if (!is.null(results$data_full)) {
        # Ensure column names in the CSV match the final internal names
        df_to_save <- results$data_full
        # Optional: Rename columns back to more descriptive names if desired for the CSV
        # colnames(df_to_save) <- c("Final_SAQ", "Received_Revasc", "Assigned_Intervention", "Baseline_SAQ", "Region")
        write.csv(df_to_save, file, row.names = FALSE, fileEncoding = "UTF-8")
      } else {
        write.csv(data.frame(Error = "No simulation data generated yet."), file, row.names = FALSE)
      }
    }
  )

}

# --- Run the application ---
# Keep this line commented out when providing the code block.
# shinyApp(ui = ui, server = server)

