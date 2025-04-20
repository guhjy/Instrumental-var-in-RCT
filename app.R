# --- Load Necessary Libraries ---
# shiny: Web application framework for R
# AER: Provides functions for applied econometrics, including ivreg() for IV estimation.
# sandwich: Provides functions for robust covariance matrix estimators.
# lmtest: Provides functions for hypothesis testing, including coeftest() for displaying coefficients with robust SEs.
# DT: Provides an interface to the DataTables library for interactive tables.

# Reference for IV concepts in RCTs:
# Hernán MA, Robins JM. Instruments for Causal Inference: An Epidemiologist’s Dream? Epidemiology. 2006;17(4):360-372. (General reference)
# Baiocchi M, Cheng J, Small DS. Instrumental Variable Methods for Causal Inference. Stat Med. 2014;33(13):2297-2340. (General reference)
# Note: The user also provided a specific reference: https://www.dropbox.com/scl/fi/lsnmm2ge2q1ircugjgqsf/nejm_iv.pdf?rlkey=8xbouk44mnryh2o98fdjbvlfx&raw=1

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
      # NEW: Added download button for the simulated data
      downloadButton("download_data", "下載模擬數據 (.csv)", class = "btn-success"),
      hr(), # Horizontal line


      h5("說明:"),
      p("此應用程式演示了在隨機對照試驗 (RCT) 中，當存在治療方案不遵從 (non-compliance) 時，如何使用工具變數 (Instrumental Variables, IV) 方法來估計治療的因果效應。"),
      p("模擬數據旨在反映類似 ISCHEMIA 試驗的情況，其中隨機分派 (Z) 作為工具變數。"),
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
                 p("當存在不遵從性時，直接比較實際接受治療 (T=1) 和未接受治療 (T=0) 的患者（即『依治療分析』, as-treated analysis）可能會因為『選擇偏誤』而產生誤導，因為決定是否接受治療的因素可能也與結果 (Y) 相關。"),
                 p("工具變數 (IV) 分析利用最初的「隨機分派」(Z) 作為「工具」，來估計「實際接受治療」(T) 對「結果」(Y) 的因果效應。在 RCT 的背景下，IV 估計的是『局部平均處理效應』(Local Average Treatment Effect, LATE)，也稱為『遵從者平均因果效應』(Complier Average Causal Effect, CACE)。"),
                 p("LATE/CACE 是指治療對那些其治療狀態因隨機分派而改變的『遵從者』(Compliers) 的平均效應。遵從者是指如果被分派到介入組就會接受治療，如果被分派到對照組就不會接受治療的個體。"),
                 hr(),
                 h4("工具變數分析的核心假設:"),
                 tags$ul(
                   tags$li(strong("1. 關聯性 (Relevance):"), "隨機分派 (Z) 必須與實際接受的治療 (T) 相關。也就是說，分派到不同組別會影響患者接受治療的可能性。"),
                   tags$li(strong("2. 獨立性/隨機性 (Independence/Randomization):"), "隨機分派 (Z) 與任何影響結果 Y 的未觀察混淆因素無關。這是由 RCT 的設計所保證的。"),
                   tags$li(strong("3. 排除限制 (Exclusion Restriction):"), "隨機分派 (Z) 只能透過影響實際接受的治療 (T) 來影響結果 (Y)，不能有其他直接影響結果的路徑（例如，分派本身引起的安慰劑效應，除非該效應也透過實際治療實現）。"),
                   tags$li(strong("4. 單調性 (Monotonicity):"), "不存在『反抗者』(Defiers)，即沒有患者會在被分派到介入組時反而不接受治療，而在被分派到對照組時反而接受治療。")
                 ),
                 hr(),
                 h4("在此模擬中使用的變數："),
                 tags$ul(
                   # UPDATED: Changed "血運重建" to "血管重建"
                   tags$li(strong("Y:"), " 結果變數 (模擬一年後的 SAQ 生活品質分數)"),
                   tags$li(strong("T:"), " 內生變數 (實際接受的血管重建治療，1=是, 0=否) - 受不遵從性影響"),
                   tags$li(strong("Z:"), " 工具變數 (隨機分派的治療策略，1=介入性, 0=保守性) - 不受患者行為影響"),
                   tags$li(strong("X:"), " 控制變數/共變數 (模擬基線 SAQ 分數 `baseline_saq`、地區 `region`)")
                 ),
                 hr(),
                 h4("數據生成模型 (Ground Truth Formula for Outcome Y)"),
                 p("模擬數據中的結果變數 (Y) 是根據以下公式生成的（注意：此公式描述了如果治療 T 被外生決定時的關係）："),
                 uiOutput("ground_truth_formula_display"), # Use uiOutput to render HTML formula
                 p("其中 ε 是隨機誤差項，服從常態分佈 N(0, σ²)。β_T 代表接受治療相對於未接受治療的真實因果效應（對於遵從者而言）。"),
                 h5("目前模擬使用的參數值:"),
                 verbatimTextOutput("ground_truth_params_display"), # Display parameters used
                 hr(),
                 h4("主要分析方法與公式："), # UPDATED: Added "與公式"
                 tags$ol(
                   tags$li(
                       strong("第一階段 (First Stage):"), " 迴歸實際接受的治療 (T) 對隨機分派 (Z) 和共變數 (X)。Z 的係數估計了隨機分派對實際接受治療機率的影響，這代表了『遵從者』在樣本中的比例（假設單調性）。此階段用於檢驗工具變數的『關聯性』假設（Z 的係數是否顯著不為零）和強度（F 統計量是否足夠大，通常 > 10）。",
                       # NEW: Added formula
                       uiOutput("formula_first_stage")
                   ),
                   tags$li(
                       strong("簡化式 (Reduced Form) / 意向治療 (Intention-to-Treat, ITT):"), " 迴歸結果 (Y) 對隨機分派 (Z) 和共變數 (X)。Z 的係數估計了『分派』到介入組相對於對照組對結果的平均影響。這是 ITT 效應，它保留了隨機分派的優點，但評估的是『分派策略』而非『實際治療』的效應。",
                       # NEW: Added formula
                       uiOutput("formula_reduced_form")
                   ),
                   tags$li(
                       strong("兩階段最小平方法 (2SLS):"), " 結合第一階段和簡化式的資訊，估計實際治療 (T) 對結果 (Y) 的『局部平均處理效應』(LATE/CACE)。這是對『遵從者』的治療效應估計。計算上，LATE ≈ ITT 效應 / 第一階段效應。",
                       # NEW: Added formula
                       uiOutput("formula_2sls")
                   ),
                   tags$li(
                       strong("普通最小平方法 (OLS) / 依治療分析 (As-Treated):"), " 直接迴歸結果 (Y) 對實際接受的治療 (T) 和共變數 (X)。當存在不遵從性時，此分析通常因忽略了 T 的內生性（選擇偏誤）而產生偏誤，其結果解釋為關聯性而非因果效應。",
                       # NEW: Added formula
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
                 p("F 統計量檢驗了工具變數 Z（在此例中只有 Z）的顯著性。若 F 統計量遠大於 10，則表示工具變數的關聯性較強，較不易受弱工具變數問題影響。"),
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
                 p("使用穩健標準誤 (Robust Standard Errors) 可以使推論在異方差 (heteroskedasticity) 存在時更為可靠，這在 IV 估計中是常見的建議做法。"),
                 verbatimTextOutput("iv_robust_summary"),
                 hr(),
                 h4("普通最小平方法迴歸 (OLS / As-Treated：T 對 Y 的關聯性 - 可能有偏誤)"),
                 p("此模型直接比較實際接受治療 (T=1) 與未接受治療 (T=0) 患者的結果 Y，控制了共變數 X。然而，由於 T 的選擇可能與影響 Y 的未觀察因素有關（內生性），T 的係數估計的是關聯性，通常不是無偏的因果效應估計。比較此 OLS 係數與 IV (LATE) 係數，可以看出不遵從性可能導致的偏誤大小和方向。"),
                 verbatimTextOutput("ols_summary"),
                 hr(), # NEW: Added separator
                 # NEW: Added section for comparing estimates to ground truth
                 h4("模型估計值與真實效應比較"),
                 uiOutput("comparison_to_truth") # Use uiOutput for dynamic text
        ),
        tabPanel("完整模擬數據", # New Tab for full data
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
    # This is the true effect for compliers (LATE/CACE) used in data generation
    beta_T_true_late <- input$beta_T_sim
    # Probability of receiving T=1 if assigned Z=1 (Intervention)
    prob_revasc_if_Z1 <- input$compliance_inv
    # Probability of receiving T=1 if assigned Z=0 (Control)
    prob_revasc_if_Z0 <- input$compliance_con

    # Z: Random assignment (Instrument)
    assignment_Z <- rbinom(n, 1, 0.5) # 1 = Intervention, 0 = Control

    # X: Covariates
    baseline_saq <- rnorm(n, mean = 60, sd = 15)
    region <- factor(sample(1:3, n, replace = TRUE), labels = c("地區A", "地區B", "地區C"))

    # T: Treatment Received (Endogenous Variable) - based on assignment Z and compliance probabilities
    # Probability of receiving treatment T=1 depends on assigned group Z
    prob_revasc <- ifelse(assignment_Z == 1, prob_revasc_if_Z1, prob_revasc_if_Z0)
    # UPDATED: Changed variable name to reflect "Revascularization"
    treatment_revasc_T <- rbinom(n, 1, prob_revasc) # 1 = Received Revascularization, 0 = Did not

    # Y: Outcome (e.g., SAQ QoL Score at 1 year) - Ground Truth Formula
    # Note: The simulation generates data consistent with the LATE being beta_T_true_late.
    # The actual generation process might implicitly define effects for always-takers/never-takers,
    # but the IV analysis targets the effect among compliers.
    intercept <- 76      # Baseline mean for control group (Z=0, T=0) if covariates were 0 and no treatment effect
    beta_baseline <- 0.1 # Effect of baseline SAQ
    beta_regionB <- 2    # Effect of Region B relative to A
    beta_regionC <- -1   # Effect of Region C relative to A
    error_sd <- 20       # Standard deviation of the error term

    # Calculate outcome Y based on the ground truth formula
    # The effect beta_T_true_late applies *when* T=1 occurs due to assignment (implicitly for compliers)
    # A simplified way to generate data consistent with a target LATE:
    # We assume the effect beta_T_true_late applies to T.
    # The endogeneity comes from the fact that T is correlated with unobservables influencing Y,
    # but IV uses Z to isolate the variation in T that is exogenous.
    outcome_Y <- intercept +
                 beta_T_true_late * treatment_revasc_T +  # Applying the effect to the observed T
                 beta_baseline * baseline_saq +
                 ifelse(region == "地區B", beta_regionB, 0) +
                 ifelse(region == "地區C", beta_regionC, 0) +
                 rnorm(n, mean = 0, sd = error_sd) # Random error term (ε)
                 # In a more complex simulation, Y might depend on latent compliance type.
                 # This simpler model still allows IV to recover beta_T_true_late under standard assumptions.

    # Create data frame
    ischemia_sim <- data.frame(
      Y = outcome_Y,
      # UPDATED: Renamed T variable
      T_revasc = treatment_revasc_T,
      Z_assign = assignment_Z, # Renamed for clarity
      baseline_saq = baseline_saq,
      region = region
    )
    # Rename columns for model formulas (use simple T and Z)
    colnames(ischemia_sim) <- c("Y", "T", "Z", "baseline_saq", "region")


    # --- 2. Perform Analyses ---
    # First Stage Regression (T ~ Z + Covariates) - Assess relevance of Z for T
    first_stage_model <- lm(T ~ Z + baseline_saq + region, data = ischemia_sim)
    first_stage_summary_obj <- summary(first_stage_model) # Store the summary object

    # Reduced Form Regression (Y ~ Z + Covariates) - Estimate ITT effect
    reduced_form_model <- lm(Y ~ Z + baseline_saq + region, data = ischemia_sim)
    reduced_form_summary_obj <- summary(reduced_form_model) # Store the summary object

    # IV Regression (2SLS: Y ~ T + Covariates | Z + Covariates) - Estimate LATE
    iv_model <- tryCatch({
        # Formula: Outcome ~ Endogenous Regressor(s) + Exogenous Regressors | Instrument(s) + Exogenous Regressors
        ivreg(Y ~ T + baseline_saq + region | Z + baseline_saq + region, data = ischemia_sim)
      }, error = function(e) {
        showNotification(paste("IV 估計錯誤:", e$message), type = "error", duration = 5)
        return(NULL) # Return NULL if IV fails
      })

    iv_summary_obj <- if (!is.null(iv_model)) summary(iv_model) else "IV 模型估計失敗" # Store summary or error message

    # Naive OLS Regression (Y ~ T + Covariates) - As-treated analysis (likely biased)
    naive_ols_model <- lm(Y ~ T + baseline_saq + region, data = ischemia_sim)
    ols_summary_obj <- summary(naive_ols_model) # Store the summary object

    # Calculate Robust Standard Errors for IV model if it exists
    iv_robust_summary_output <- NULL
    if (!is.null(iv_model) && !is.character(iv_model)) { # Check it's a model object
      iv_robust_summary_output <- tryCatch({
          coeftest(iv_model, vcov. = vcovHC(iv_model, type = "HC1")) # HC1 is common
        }, error = function(e) {
          showNotification(paste("計算穩健標準誤時發生錯誤:", e$message), type = "warning", duration = 5)
          return(NULL) # Return NULL if robust SE calculation fails
      })
    } else {
        iv_robust_summary_output <- "IV 模型不存在或估計失敗，無法計算穩健標準誤。"
    }


    # Store ground truth parameters used in this simulation run
    ground_truth_params <- list(
        n = n,
        intercept = intercept,
        beta_T_true_late = beta_T_true_late, # True LATE/CACE
        beta_baseline = beta_baseline,
        beta_regionB = beta_regionB,
        beta_regionC = beta_regionC,
        error_sd = error_sd,
        prob_revasc_if_Z1 = prob_revasc_if_Z1, # P(T=1|Z=1)
        prob_revasc_if_Z0 = prob_revasc_if_Z0  # P(T=1|Z=0)
    )

    # Return all results as a list
    list(
      data_full = ischemia_sim,
      first_stage_summary = first_stage_summary_obj,
      reduced_form_summary = reduced_form_summary_obj,
      iv_model = iv_model, # Return the model object itself for coefficient extraction
      iv_summary = iv_summary_obj,
      iv_robust_summary = iv_robust_summary_output,
      ols_model = naive_ols_model, # Return the model object
      ols_summary = ols_summary_obj,
      ground_truth_params = ground_truth_params
    )
  }, ignoreNULL = FALSE) # ignoreNULL = FALSE makes it run once automatically

  # --- 3. Render Outputs ---

  # Render data table preview (first 6 rows)
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

    # Extract F-statistic for the instrument(s) relevance check
    f_stat <- results$first_stage_summary$fstatistic
    # Check if f_stat is available (it should be for lm summary)
    if(!is.null(f_stat) && length(f_stat) >= 1){
      f_value <- f_stat[1] # The F value is the first element
      cat("\n---\n")
      cat(paste("第一階段 F 統計量 (檢定工具 Z 的聯合顯著性):", round(f_value, 2), "\n"))
      if(f_value < 10){
         cat("警告：第一階段 F 統計量 < 10。可能存在弱工具變數問題，IV 估計可能不可靠。\n")
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
        cat(results$iv_summary) # Print error message if estimation failed
    } else {
        print(results$iv_summary) # Print the summary.ivreg object
    }
  })

  # Render IV model summary with robust SEs
  output$iv_robust_summary <- renderPrint({
     results <- analysis_results()
     req(results$iv_robust_summary)
     if (is.character(results$iv_robust_summary)) {
         cat(results$iv_robust_summary) # Print error/info message
     } else {
         print(results$iv_robust_summary) # Print the coeftest object
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
      withMathJax(helpText(
        # Updated formula explanation
        "$$ Y = \\beta_0 + \\beta_{T, LATE} T + \\beta_{baseline} \\text{baseline\\_saq} + \\beta_{regionB} I(\\text{region=B}) + \\beta_{regionC} I(\\text{region=C}) + \\epsilon $$"
      ))
  })

  # NEW: Render 1st Stage Formula
  output$formula_first_stage <- renderUI({
      withMathJax(helpText(
        "$$ T = \\gamma_0 + \\gamma_Z Z + \\gamma_X X + \\nu $$"
      ))
  })

  # NEW: Render Reduced Form (ITT) Formula
  output$formula_reduced_form <- renderUI({
      withMathJax(helpText(
        "$$ Y = \\pi_0 + \\pi_Z Z + \\pi_X X + \\omega $$"
      ))
  })

  # NEW: Render 2SLS Formula Explanation
  output$formula_2sls <- renderUI({
      withMathJax(helpText(
        "概念上，第一階段預測 \\( \\hat{T} = \\hat{\\gamma}_0 + \\hat{\\gamma}_Z Z + \\hat{\\gamma}_X X \\)。第二階段使用 \\( \\hat{T} \\) 取代 T:",
        "$$ Y = \\beta_0 + \\beta_{T, 2SLS} \\hat{T} + \\beta_X X + \\epsilon $$",
        "其中 \\( \\beta_{T, 2SLS} \\approx \\frac{\\hat{\\pi}_Z}{\\hat{\\gamma}_Z} \\) (LATE 估計值)。"
      ))
  })

  # NEW: Render OLS Formula
  output$formula_ols <- renderUI({
      withMathJax(helpText(
        "$$ Y = \\delta_0 + \\delta_T T + \\delta_X X + \\mu $$"
      ))
  })


  # Render the ground truth parameters used in the simulation
  output$ground_truth_params_display <- renderPrint({
      results <- analysis_results()
      req(results$ground_truth_params)
      params <- results$ground_truth_params
      # Calculate proportion of compliers from input probabilities
      prop_compliers = params$prob_revasc_if_Z1 - params$prob_revasc_if_Z0

      cat(paste("樣本數 (n):", params$n, "\n"))
      cat(paste("截距項 (β₀):", params$intercept, "\n"))
      cat(paste("真實 LATE/CACE (β_T):", params$beta_T_true_late, "\n"))
      cat(paste("基線 SAQ 效應 (β_baseline):", params$beta_baseline, "\n"))
      cat(paste("地區 B 效應 (β_regionB):", params$beta_regionB, "\n"))
      cat(paste("地區 C 效應 (β_regionC):", params$beta_regionC, "\n"))
      cat(paste("誤差標準差 (σ_ε):", params$error_sd, "\n"))
      cat(paste("P(T=1 | Z=1) (介入組治療率):", params$prob_revasc_if_Z1, "\n"))
      cat(paste("P(T=1 | Z=0) (對照組治療率):", params$prob_revasc_if_Z0, "\n"))
      cat(paste("遵從者比例估計 (P(T=1|Z=1) - P(T=1|Z=0)):", round(prop_compliers, 3), "\n"))
  })

  # NEW: Render comparison of estimates to ground truth
  output$comparison_to_truth <- renderUI({
      results <- analysis_results()
      req(results$iv_model, results$ols_model, results$ground_truth_params)

      # Extract coefficients safely
      iv_coef <- tryCatch(coef(results$iv_model)["T"], error = function(e) NA)
      ols_coef <- tryCatch(coef(results$ols_model)["T"], error = function(e) NA)
      true_late <- results$ground_truth_params$beta_T_true_late

      if (is.na(iv_coef) || is.na(ols_coef)) {
          return(p("無法提取模型係數進行比較。"))
      }

      # Create comparison text
      tagList(
          p(paste0("在此模擬中，真實的遵從者平均因果效應 (LATE/CACE) 設定為: ", round(true_late, 3))),
          p(paste0("2SLS 估計的治療 (T) 對結果 (Y) 的效應 (LATE 估計值) 為: ", round(iv_coef, 3))),
          p(paste0("OLS (依治療分析) 估計的治療 (T) 對結果 (Y) 的關聯性為: ", round(ols_coef, 3))),
          hr(),
          p("比較："),
          tags$ul(
              tags$li(paste0("2SLS 估計值 (", round(iv_coef, 3), ") 通常較接近真實的 LATE (", round(true_late, 3), ")，因為它利用隨機分派 (Z) 作為工具變數來處理因不遵從性引起的內生性問題。")),
              tags$li(paste0("OLS 估計值 (", round(ols_coef, 3), ") 可能與真實 LATE 相差較大。這是因為 OLS 比較的是實際接受治療與未接受治療的兩組，而這兩組可能因自我選擇或其他未觀察因素而在結果上存在差異（選擇偏誤），OLS 未能校正此偏誤。在此模擬中，OLS 估計值與真實 LATE 的差異反映了這種偏誤的大小和方向。"))
          )
      )
  })

  # --- Download Handler ---
  # NEW: Added download handler for the data
  output$download_data <- downloadHandler(
    filename = function() {
      paste0("simulated_ischemia_data_", Sys.Date(), ".csv")
    },
    content = function(file) {
      results <- analysis_results()
      if (!is.null(results$data_full)) {
        write.csv(results$data_full, file, row.names = FALSE, fileEncoding = "UTF-8") # Ensure UTF-8 encoding for compatibility
      } else {
        # Create an empty file or write an error message if data is not available
        write.csv(data.frame(Error = "No simulation data generated yet."), file, row.names = FALSE)
      }
    }
  )

}

# --- Run the application ---
# Keep this line commented out when providing the code block.
shinyApp(ui = ui, server = server)
