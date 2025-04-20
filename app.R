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

# --- UI Definition (使用者介面) ---
ui <- fluidPage(
  tags$head(
    tags$meta(charset = "UTF-8"),
    tags$html(lang = "zh-Hant-TW"),
    tags$script(src = "https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.7/MathJax.js?config=TeX-MML-AM_CHTML")
  ),
  titlePanel("有許多非遵從者之隨機分配對照臨床試驗 RCT 的工具變數分析"),
  sidebarLayout(
    sidebarPanel(
      h4("控制面板"),
      h5("模擬參數設定:"),
      numericInput("n_sim", "樣本數 (Sample Size, n):", value = 5000, min = 100, max = 20000, step = 100),
      numericInput("beta_T_sim", "真實治療效應 (LATE/CACE):", value = 5.5, min = -10, max = 20, step = 0.5),
      sliderInput("compliance_inv", "遵從率 - 介入組 (P(T=1|Z=1)):", min = 0, max = 1, value = 0.80, step = 0.01),
      sliderInput("compliance_con", "非遵從率 - 對照組 (P(T=1|Z=0)):", min = 0, max = 1, value = 0.12, step = 0.01),
      hr(),
      h5("模型設定:"),
      numericInput("rcs_knots", "基線 SAQ 的 RCS 結數 (Number of Knots for RCS):", value = 4, min = 3, max = 7, step = 1),
      helpText("設為 0 或 1 或 2 則使用線性項。建議 3-5 個結。"), # Updated help text slightly
      hr(),
      actionButton("run_analysis", "重新模擬數據並執行分析", icon = icon("play"), class = "btn-primary"),
      hr(),
      downloadButton("download_data", "下載模擬數據 (.csv)", class = "btn-success"),
      hr(),
      h5("說明:"),
      p("此應用程式演示 RCT 中存在不遵從性時，如何使用 IV 方法估計治療效應。"),
      p("模擬數據旨在反映類似 ISCHEMIA 試驗的情況。"),
      p("可選擇使用受限三次樣條 (Restricted Cubic Splines, RCS) 來模擬基線 SAQ 的非線性效應。"),
      p("IV 模型結果現在包含診斷檢定。"),
      p("OLS 模型結果包含比較線性與 RCS 模型的 ANOVA 檢定，以評估非線性。") # Added explanation
    ),
    mainPanel(
      tabsetPanel(
        id = "results_tabs",
        tabPanel("分析說明與模擬設定",
                 h3("工具變數 (IV) 與隨機試驗"),
                 p("當存在不遵從性時，直接比較實際接受治療 (T=1) 和未接受治療 (T=0) 的患者（即『依治療分析』, OLS）可能會因為『選擇偏誤』而產生誤導，因為決定接受治療的因素可能也與結果相關。"),
                 p("工具變數 (IV) 分析利用最初的「隨機分派」(Z) 作為「工具」，來估計「實際接受治療」(T) 對「結果」(Y) 的因果效應。在存在不遵從性的 RCT 中，這通常估計的是『遵從者的平均因果效應』(Complier Average Causal Effect, CACE)，也稱為『局部平均治療效應』(Local Average Treatment Effect, LATE)。"),
                 hr(),
                 h4("工具變數分析的核心假設:"),
                 tags$ul(
                   tags$li("1. 關聯性: 工具變數 Z 與內生變數（跟誤差項相關的右側變項）T 相關 (在控制共變數後)。這是可以檢驗的 (見『遵從性檢查』分頁)。"),
                   tags$li("2. 獨立性/隨機性: 工具變數 Z 與結果 Y 的誤差項無關，亦即 Z 是外生性變數。在 RCT 中，隨機分派通常能滿足此假設。"),
                   tags$li("3. 排除限制（Z →X→Y）: 工具變數 Z 只能通過內生變數 T 影響結果 Y，沒有直接經由或其他間接路徑 (在控制共變數後) 影響 Y。這個假設無法直接檢驗，但在 RCT 中通常被認為是合理的 (亦即分派本身不影響結果)。"),
                   tags$li("4. 單調性：沒有『反抗者』，即沒有人會因為被分派到介入組反而不接受治療，同時如果被分派到對照組反而接受治療。這意味著 Z 對 T 的影響方向對所有人都是一致的 (或者為零)。")
                 ),
                 hr(),
                 h4("在此模擬中使用的變數："),
                 tags$ul(
                   tags$li("Y: 結果變數 (一年後的 SAQ 生活品質分數, 0-100, 高分=好)"),
                   tags$li("T: 內生變數 (實際接受的血管重建治療，1=是, 0=否)"),
                   tags$li("Z: 工具變數 (隨機分派的治療策略，1=介入性, 0=保守性)"),
                   tags$li("X: 控制變數/共變數 (模擬基線 SAQ 分數 `baseline_SAQ`、地區 `region`)")
                 ),
                 h4("受限三次樣條 (Restricted Cubic Splines, RCS)"),
                 p("您可以選擇使用 RCS 來對基線 SAQ (`baseline_SAQ`) 與結果 Y 或治療 T 之間的關係進行建模。RCS 是一種靈活的方式，可以在不過度擬合數據的情況下捕捉潛在的非線性關係。結 (knots) 的數量決定了曲線的靈活性（通常選擇 3 到 5 個結）。"),
                 hr(),
                 h4("數據生成模型 (Ground Truth Formula for Outcome Y)"),
                 uiOutput("ground_truth_formula_display"),
                 h5("目前模擬使用的參數值:"),
                 verbatimTextOutput("ground_truth_params_display"),
                 hr(),
                 h4("主要分析方法與公式："),
                 tags$ol(
                   tags$li("第一階段: T ~ Z + Covariates", uiOutput("formula_first_stage")),
                   tags$li("簡化式/ ITT: Y ~ Z + Covariates, ITT ≈ LATE × (P(T=1|Z=1) - P(T=1|Z=0))", uiOutput("formula_reduced_form")),
                   tags$li("兩階段最小平方法 (2SLS / IV): Y ~ T + Covariates | Z + Covariates", uiOutput("formula_2sls")),
                   tags$li("普通最小平方法 (OLS) / As-Treated: Y ~ T + Covariates", uiOutput("formula_ols"))
                 )
        ),
        tabPanel("數據預覽與遵從性檢查",
                 h4("模擬數據預覽 (前 6 筆)"),
                 tableOutput("data_head"),
                 hr(),
                 h4("第一階段迴歸 (遵從性檢查：Z 對 T 的影響)"),
                 h5("工具變數強度與弱工具變數問題"),
                 p("工具變數的『關聯性』假設要求 Z 必須與 T 相關。我們通常使用第一階段迴歸中工具變數 Z 的 F 統計量來評估工具變數的強度。"),
                 tags$ul(
                    tags$li("計算方式：檢驗第一階段迴歸中所有工具變數（在此例中只有 Z）係數是否聯合為零的 F 檢定統計量。"),
                    tags$li("經驗法則：常用的判斷標準是 F 統計量 > 10。這表示工具變數與內生變數有足夠強的關聯。"),
                    tags$li("弱工具變數：若 F < 10，則認為工具變數 Z 與 T 的關聯性較弱，稱為『弱工具變數』。"),
                    tags$li("弱工具變數的後果："),
                    tags$ul(
                        tags$li("2SLS 估計量會產生偏誤 (Bias)，且偏誤方向會朝向 OLS 估計量。"),
                        tags$li("2SLS 估計量的標準誤會變得不可靠，導致信賴區間和假說檢定 (如 t 檢定) 的結果不準確。"),
                        tags$li("即使樣本數很大，弱工具變數問題仍然可能存在。")
                    ),
                    tags$li("如何處理弱工具變數："),
                    tags$ul(
                        tags$li("尋找更強的工具變數 (如果可能)。"),
                        tags$li("使用對弱工具變數較不敏感的估計方法，例如『有限資訊最大概似法』LIML)，但 AER 套件中的 `ivreg` 不直接提供 LIML。"),
                        tags$li("使用對弱工具變數穩健的推論方法，例如 Anderson-Rubin 檢定。"),
                        tags$li("承認弱工具變數的存在，並在解釋結果時更加謹慎。")
                    ),
                    tags$li("重要性：檢查第一階段 F 統計量是 IV 分析中非常重要的一步。`ivreg` 的 `summary` 函數加上 `diagnostics = TRUE` 會自動報告此 F 統計量。")
                 ),
                 h5("線性基線 SAQ 模型結果:"),
                 verbatimTextOutput("first_stage_summary_linear"),
                 h5("RCS 基線 SAQ 模型結果:"),
                 verbatimTextOutput("first_stage_summary_rcs")
        ),
        tabPanel("模型結果比較",
                 h4("簡化式迴歸 (ITT 效應：Z 對 Y 的影響)"),
                 h5("線性基線 SAQ 模型結果:"),
                 verbatimTextOutput("reduced_form_summary_linear"),
                 h5("RCS 基線 SAQ 模型結果:"),
                 verbatimTextOutput("reduced_form_summary_rcs"),
                 hr(),
                 h4("工具變數迴歸 (2SLS 結果：T 對 Y 的 LATE/CACE)"),
                 p("注意：以下結果包含使用 `summary(..., diagnostics = TRUE)` 產生的診斷檢定。"),
                 tags$ul(
                    tags$li(strong("弱工具變數檢定:"), "報告第一階段的 F 統計量 (同『遵從性檢查』分頁)。F > 10 通常表示工具變數強度足夠。"),
                    tags$li(strong("內生性檢定 (Wu-Hausman Test):"), "檢定內生變數 (T) 是否真的具有內生性 (即 OLS 是否有偏誤)。虛無假設 (H0) 是 T 為外生。若 p 值顯著 (例如 < 0.05)，則拒絕 H0，支持 T 具有內生性，使用 IV 是適當的。"),
                    tags$li(strong("過度識別檢定 (Sargan Test of Overidentifying Restrictions):"), "僅在工具變數數量 > 內生變數數量時適用 (在此模擬中不適用，因為只有一個工具 Z 和一個內生變數 T)。檢定『額外的』工具變數是否有效 (即是否滿足排除限制且與誤差項無關)。虛無假設 (H0) 是所有工具變數都有效。若 p 值不顯著，則無法拒絕 H0，支持工具變數的有效性。")
                 ),
                 h5("線性基線 SAQ 模型結果 (含診斷):"),
                 verbatimTextOutput("iv_summary_linear"),
                 h5("RCS 基線 SAQ 模型結果 (含診斷):"),
                 verbatimTextOutput("iv_summary_rcs"),
                 hr(),
                 h4("工具變數迴歸 (含穩健標準誤)"),
                 p("注意：使用 `coeftest` 搭配 `vcovHC` 計算穩健標準誤時，不會顯示 `summary` 的診斷檢定。診斷檢定應參考上面的標準 IV 結果。"),
                 h5("線性基線 SAQ 模型結果 (穩健 SE):"),
                 verbatimTextOutput("iv_robust_summary_linear"),
                 h5("RCS 基線 SAQ 模型結果 (穩健 SE):"),
                 verbatimTextOutput("iv_robust_summary_rcs"),
                 hr(),
                 h4("普通最小平方法迴歸 (OLS / As-Treated)"),
                 h5("線性基線 SAQ 模型結果:"),
                 verbatimTextOutput("ols_summary_linear"),
                 h5("RCS 基線 SAQ 模型結果:"),
                 verbatimTextOutput("ols_summary_rcs"),
                 h5("OLS 中 Baseline SAQ 非線性檢定 (ANOVA: 線性 vs RCS): 即使 p >0.05, 也要保留 RCS"), # New Title for ANOVA test
                 verbatimTextOutput("ols_nonlinearity_test_output"),          # New Output for ANOVA test
                 hr(),
                 h4("模型估計值比較 (基於線性 SAQ 模型)"),
                 uiOutput("comparison_estimates")
        ),
        tabPanel("完整模擬數據",
                 h4("完整的模擬數據集"),
                 DT::dataTableOutput("full_data_table")
        ),
        tabPanel("參考文獻",
                 h4("主要參考文獻"),
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
      )
    )
  )
)


# --- Server Logic (伺服器邏輯) ---
server <- function(input, output, session) {

  analysis_results <- eventReactive(input$run_analysis, {

    showNotification("正在模擬數據並執行分析...", type = "message", duration = 3)

    # --- 1. Simulate Data ---
    set.seed(as.integer(Sys.time()) %% 10000) # Use modulo to keep seed within reasonable integer range
    n <- input$n_sim
    beta_T_true_late <- input$beta_T_sim
    prob_revasc_if_Z1 <- input$compliance_inv
    prob_revasc_if_Z0 <- input$compliance_con
    nk <- input$rcs_knots # Get number of knots for RCS
    use_rcs <- nk >= 3 # Determine if RCS should be used based on knots

    # Z: Random assignment
    assignment_Z <- rbinom(n, 1, 0.5)

    # X: Covariates
    baseline_SAQ <- rnorm(n, mean = 60, sd = 15)
    region <- factor(sample(1:3, n, replace = TRUE), labels = c("地區A", "地區B", "地區C"))

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
    # confounding_effect <- 0.05 * (baseline_SAQ - 60) # Confounding removed for clarity, OLS bias comes from T correlation

    outcome_Y <- intercept +
                 beta_T_true_late * treatment_revasc_T +
                 beta_baseline * baseline_SAQ + # Effect of baseline SAQ
                 ifelse(region == "地區B", beta_regionB, 0) +
                 ifelse(region == "地區C", beta_regionC, 0) +
                 # confounding_effect + # Removed confounding effect
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
    formula_rcs_cov_terms <- if (use_rcs) paste("rcs(baseline_SAQ, nk=", nk, ") + region") else formula_linear_cov_terms

    # Define formulas (Linear and RCS)
    formula_fs_linear <- as.formula(paste("T ~ Z +", formula_linear_cov_terms))
    formula_rf_linear <- as.formula(paste("Y ~ Z +", formula_linear_cov_terms))
    formula_ols_linear <- as.formula(paste("Y ~ T +", formula_linear_cov_terms))
    formula_iv_linear <- as.formula(paste("Y ~ T +", formula_linear_cov_terms, "| Z +", formula_linear_cov_terms))

    formula_fs_rcs <- if (use_rcs) as.formula(paste("T ~ Z +", formula_rcs_cov_terms)) else NULL
    formula_rf_rcs <- if (use_rcs) as.formula(paste("Y ~ Z +", formula_rcs_cov_terms)) else NULL
    formula_ols_rcs <- if (use_rcs) as.formula(paste("Y ~ T +", formula_rcs_cov_terms)) else NULL
    formula_iv_rcs <- if (use_rcs) as.formula(paste("Y ~ T +", formula_rcs_cov_terms, "| Z +", formula_rcs_cov_terms)) else NULL


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

    iv_linear_summary <- if (!is.null(iv_linear_model) && inherits(iv_linear_model, "ivreg")) {
        tryCatch({
            summary(iv_linear_model, diagnostics = TRUE)
        }, error = function(e) { message("IV Linear Summary Error: ", e$message); "IV (Linear) 模型摘要失敗 (含診斷)"})
    } else {
        "IV (Linear) 模型估計失敗"
    }


    iv_linear_robust_summary <- NULL
    if (!is.null(iv_linear_model) && inherits(iv_linear_model, "ivreg")) {
      iv_linear_robust_summary <- tryCatch({
          coeftest(iv_linear_model, vcov. = vcovHC(iv_linear_model, type = "HC1"))
        }, error = function(e) { message("Robust SE (Linear) Error: ", e$message); NULL })
      if(is.null(iv_linear_robust_summary)) iv_linear_robust_summary <- "無法計算穩健標準誤 (Linear)。"
    } else {
        iv_linear_robust_summary <- "無法計算穩健標準誤 (Linear)。"
    }

    # --- Run RCS Models (if use_rcs is TRUE) ---
    first_stage_rcs_model <- NULL
    first_stage_rcs_summary <- "未執行 RCS 模型 (結數 < 3)。"
    reduced_form_rcs_model <- NULL
    reduced_form_rcs_summary <- "未執行 RCS 模型 (結數 < 3)。"
    ols_rcs_model <- NULL
    ols_rcs_summary <- "未執行 RCS 模型 (結數 < 3)。"
    iv_rcs_model <- NULL
    iv_rcs_summary <- "未執行 RCS 模型 (結數 < 3)。"
    iv_rcs_robust_summary <- "未執行 RCS 模型 (結數 < 3)。"
    nonlinearity_test_ols <- "未執行 RCS 模型 (結數 < 3)，無法進行非線性檢定。" # Initialize nonlinearity test result

    if (use_rcs) { # Check if RCS should be used
        first_stage_rcs_model <- tryCatch({ lm(formula_fs_rcs, data = ischemia_sim) }, error = function(e){ message("FS RCS Error: ", e$message); NULL})
        first_stage_rcs_summary <- if (!is.null(first_stage_rcs_model)) summary(first_stage_rcs_model) else "第一階段 (RCS) 模型估計失敗"

        reduced_form_rcs_model <- tryCatch({ lm(formula_rf_rcs, data = ischemia_sim) }, error = function(e){ message("RF RCS Error: ", e$message); NULL})
        reduced_form_rcs_summary <- if (!is.null(reduced_form_rcs_model)) summary(reduced_form_rcs_model) else "簡化式 (RCS) 模型估計失敗"

        ols_rcs_model <- tryCatch({ lm(formula_ols_rcs, data = ischemia_sim) }, error = function(e){ message("OLS RCS Error: ", e$message); NULL})
        ols_rcs_summary <- if (!is.null(ols_rcs_model)) summary(ols_rcs_model) else "OLS (RCS) 模型估計失敗"

        # *** ADDED: Perform ANOVA test for nonlinearity in OLS ***
        if (!is.null(ols_linear_model) && !is.null(ols_rcs_model)) {
             nonlinearity_test_ols <- tryCatch({
                 anova(ols_linear_model, ols_rcs_model) # Compare linear vs RCS OLS models
             }, error = function(e) {
                 message("ANOVA comparison error (OLS): ", e$message)
                 "無法執行 OLS 線性 vs RCS 模型比較。"
             })
        } else {
             nonlinearity_test_ols <- "OLS 線性或 RCS 模型之一未成功擬合，無法比較。"
        }
        # *** END OF ADDED PART ***

        iv_rcs_model <- tryCatch({
            ivreg(formula_iv_rcs, data = ischemia_sim)
        }, error = function(e) { message("IV RCS Model Error: ", e$message); NULL })

        iv_rcs_summary <- if (!is.null(iv_rcs_model) && inherits(iv_rcs_model, "ivreg")) {
             tryCatch({
                summary(iv_rcs_model, diagnostics = TRUE)
             }, error = function(e) { message("IV RCS Summary Error: ", e$message); "IV (RCS) 模型摘要失敗 (含診斷)"})
        } else {
             "IV (RCS) 模型估計失敗"
        }


        if (!is.null(iv_rcs_model) && inherits(iv_rcs_model, "ivreg")) {
          iv_rcs_robust_summary <- tryCatch({
              coeftest(iv_rcs_model, vcov. = vcovHC(iv_rcs_model, type = "HC1"))
            }, error = function(e) { message("Robust SE (RCS) Error: ", e$message); NULL })
           if(is.null(iv_rcs_robust_summary)) iv_rcs_robust_summary <- "無法計算穩健標準誤 (RCS)。"
        } else {
            iv_rcs_robust_summary <- "無法計算穩健標準誤 (RCS)。"
        }
    } # End of if(use_rcs)


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
      first_stage_rcs_model = first_stage_rcs_model,
      first_stage_rcs_summary = first_stage_rcs_summary,
      reduced_form_rcs_model = reduced_form_rcs_model,
      reduced_form_rcs_summary = reduced_form_rcs_summary,
      ols_rcs_model = ols_rcs_model,
      ols_rcs_summary = ols_rcs_summary,
      iv_rcs_model = iv_rcs_model,
      iv_rcs_summary = iv_rcs_summary,
      iv_rcs_robust_summary = iv_rcs_robust_summary,
      nonlinearity_test_ols = nonlinearity_test_ols, # Include ANOVA test result
      # Params
      ground_truth_params = ground_truth_params,
      nk = nk,
      use_rcs = use_rcs # Pass whether RCS was used
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
  print_first_stage_summary <- function(summary_obj, model_type = "Linear", use_rcs_flag = FALSE, knots = NULL) {
      if (is.character(summary_obj)) { # Handle error messages or "not run" messages
          cat(paste("第一階段 (", model_type, "): ", summary_obj, "\n", sep=""))
          return()
      }
      # Check if summary object is valid before printing
      req(summary_obj, inherits(summary_obj, "summary.lm"))
      cat(paste("--- 第一階段模型 (", model_type, ") ---\n", sep=""))
      if (use_rcs_flag && !is.null(knots)) {
          cat(paste("使用 rcs(baseline_SAQ, nk=", knots, ")\n", sep=""))
      } else {
          cat("使用 baseline_SAQ (線性)\n")
      }
      print(summary_obj)

      # Extract and display F-statistic for the instrument Z
      if ("Z" %in% rownames(summary_obj$coefficients)) {
          t_stat_Z <- summary_obj$coefficients["Z", "t value"]
          f_stat_Z <- t_stat_Z^2
          cat("\n---\n")
          cat(paste("工具變數 Z 的 F 統計量 (近似值 t^2):", round(f_stat_Z, 2), "\n"))
          if (f_stat_Z < 10) {
              cat("警告：第一階段 F 統計量 < 10。可能存在弱工具變數問題。\n")
          } else {
              cat("第一階段 F 統計量 >= 10，工具變數強度尚可接受。\n")
          }
      } else {
          cat("\n---\n無法提取工具變數 Z 的 F 統計量。\n")
      }
  }

  # First Stage Output
  output$first_stage_summary_linear <- renderPrint({
      results <- analysis_results()
      req(results$first_stage_linear_summary)
      print_first_stage_summary(results$first_stage_linear_summary, "Linear", use_rcs_flag = FALSE)
  })
  output$first_stage_summary_rcs <- renderPrint({
      results <- analysis_results()
      req(results$first_stage_rcs_summary) # Requires the summary object or the message string
      print_first_stage_summary(results$first_stage_rcs_summary, "RCS", use_rcs_flag = results$use_rcs, knots = results$nk)
  })

  # Helper function to print model summaries safely
  print_model_summary <- function(summary_obj, model_name) {
      cat(paste("--- ", model_name, " ---\n", sep=""))
      if(is.character(summary_obj)) {
          cat(summary_obj)
      } else if (inherits(summary_obj, c("summary.lm", "summary.ivreg", "coeftest"))) {
          print(summary_obj)
      } else {
          cat("無法顯示模型摘要。")
      }
      cat("\n") # Add a newline for spacing
  }

  # Reduced Form Output
  output$reduced_form_summary_linear <- renderPrint({
      results <- analysis_results(); req(results$reduced_form_linear_summary)
      print_model_summary(results$reduced_form_linear_summary, "簡化式模型 (Linear)")
  })
  output$reduced_form_summary_rcs <- renderPrint({
      results <- analysis_results(); req(results$reduced_form_rcs_summary)
      print_model_summary(results$reduced_form_rcs_summary, paste("簡化式模型 (RCS, nk=", results$nk, ")", sep=""))
  })

  # IV / 2SLS Output (Now includes diagnostics)
  output$iv_summary_linear <- renderPrint({
    results <- analysis_results(); req(results$iv_linear_summary)
    print_model_summary(results$iv_linear_summary, "IV/2SLS 模型 (Linear, 含診斷)")
  })
  output$iv_summary_rcs <- renderPrint({
    results <- analysis_results(); req(results$iv_rcs_summary)
    print_model_summary(results$iv_rcs_summary, paste("IV/2SLS 模型 (RCS, nk=", results$nk, ", 含診斷)", sep=""))
  })

  # IV / 2SLS Robust Output (No diagnostics here)
  output$iv_robust_summary_linear <- renderPrint({
     results <- analysis_results(); req(results$iv_linear_robust_summary)
     print_model_summary(results$iv_linear_robust_summary, "IV/2SLS 模型 (Linear, 穩健 SE)")
  })
   output$iv_robust_summary_rcs <- renderPrint({
     results <- analysis_results(); req(results$iv_rcs_robust_summary)
     print_model_summary(results$iv_rcs_robust_summary, paste("IV/2SLS 模型 (RCS, nk=", results$nk, ", 穩健 SE)", sep=""))
  })

  # OLS Output
  output$ols_summary_linear <- renderPrint({
     results <- analysis_results(); req(results$ols_linear_summary)
     print_model_summary(results$ols_linear_summary, "OLS 模型 (Linear)")
  })
  output$ols_summary_rcs <- renderPrint({
     results <- analysis_results(); req(results$ols_rcs_summary)
     print_model_summary(results$ols_rcs_summary, paste("OLS 模型 (RCS, nk=", results$nk, ")", sep=""))
  })

  # *** ADDED: Render ANOVA test output for OLS nonlinearity ***
  output$ols_nonlinearity_test_output <- renderPrint({
      results <- analysis_results()
      req(results$nonlinearity_test_ols) # Requires the test result or message

      # Check if the result is the message string indicating RCS wasn't run
      if(is.character(results$nonlinearity_test_ols) && grepl("未執行 RCS 模型", results$nonlinearity_test_ols)) {
          cat(results$nonlinearity_test_ols)
          return() # Stop here if RCS wasn't run
      }

      # If we got here, RCS was run, proceed to display test info
      cat("--- OLS: Baseline SAQ 非線性檢定 ---\n")
      if(results$use_rcs) { # Ensure nk was >= 3
        cat("比較模型：\n")
        cat(paste("  模型 1 (線性): Y ~ baseline_SAQ + region + T\n"))
        cat(paste("  模型 2 (RCS):   Y ~ rcs(baseline_SAQ, nk=", results$nk, ") + region + T\n", sep=""))
        cat("虛無假設 (H0): 線性模型已足夠 (非線性項係數為 0)\n")
        cat("對立假設 (H1): 非線性項顯著改善模型\n\n")

        # Print the anova result (which could be the table or an error message)
        print(results$nonlinearity_test_ols)

        # Add interpretation hint if the anova object is valid
        if(inherits(results$nonlinearity_test_ols, "anova")) {
            # Extract P-value (Pr(>F)) from the second row if it exists
            p_value <- tryCatch(results$nonlinearity_test_ols$"Pr(>F)"[2], error = function(e) NA)
            if(!is.na(p_value)) {
                cat("\n---\n")
                cat(paste("非線性項的 P 值 (來自 ANOVA F 檢定):", format.pval(p_value, digits = 3, eps = 0.001), "\n"))
                if(p_value < 0.05) {
                    cat("結論：有顯著證據表明 Baseline SAQ 與結果 Y 之間存在非線性關係 (在 OLS 模型中，p < 0.05)。\n")
                } else {
                    cat("結論：沒有足夠證據拒絕線性關係假設 (在 OLS 模型中，p >= 0.05)。\n")
                }
            } else {
                cat("\n---\n無法從 ANOVA 結果中提取 P 值。\n")
            }
        } else if (is.character(results$nonlinearity_test_ols)) {
            # If it's a character string, it's likely an error message from tryCatch
             cat("\n---\n")
             cat(paste("ANOVA 比較時發生錯誤:", results$nonlinearity_test_ols, "\n"))
        }
      } else {
           cat("未執行 RCS 模型 (結數 < 3)，無法進行非線性檢定。\n")
      }
  })
  # *** END OF ADDED RENDER BLOCK ***


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
      cat(paste("樣本數 (n):", params$n, "\n"))
      cat(paste("真實 LATE/CACE (β_T):", params$beta_T_true_late, "\n"))
      cat(paste("真實基線 SAQ 效應 (β_baseline, 線性):", params$beta_baseline, "\n"))
      cat(paste("地區 B 效應 (β_regionB):", params$beta_regionB, "\n"))
      cat(paste("地區 C 效應 (β_regionC):", params$beta_regionC, "\n"))
      cat(paste("誤差標準差 (σ_ε):", params$error_sd, "\n"))
      cat(paste("P(T=1 | Z=1) (輸入設定):", params$prob_revasc_if_Z1, "\n")) # Note: This is the input slider value
      cat(paste("P(T=1 | Z=0) (輸入設定):", params$prob_revasc_if_Z0, "\n")) # Note: This is the input slider value
      cat(paste("遵從者比例估計 (基於輸入):", round(prop_compliers, 3), "\n"))
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
          return(p("錯誤：無法從一個或多個線性模型中提取係數進行比較。請檢查模型結果。"))
      }

      # Create comparison text
      tagList(
          p(paste0("在此模擬中，真實的遵從者平均因果效應 (LATE/CACE) 設定為: ", round(true_late, 3))),
          hr(),
          h5("不同分析方法的估計值 (基於線性 SAQ 模型)："),
          tags$ul(
              tags$li(paste0("ITT (意向治療) 效應 (Z 對 Y): ", round(itt_coef, 3), ". 隨機『分派』對結果 Y 的平均影響。")),
              tags$li(paste0("OLS (依治療分析) 關聯性 (T 對 Y): ", round(ols_coef, 3), ". 比較實際『接受』治療者的結果差異，可能有偏誤。")),
              tags$li(paste0("2SLS (工具變數) 效應 (T 對 Y 的 LATE/CACE): ", round(iv_coef, 3), ". 利用隨機分派 Z 估計實際『接受』治療 T 對『遵從者』的因果效應。"))
          ),
          hr(),
          h5("比較與解釋："),
          tags$ul(
              tags$li(paste0("2SLS 估計值 (", round(iv_coef, 3), ") 通常最接近真實 LATE (", round(true_late, 3), ")，因其試圖校正由 T 的內生性（選擇偏誤）引起的偏誤。")),
              tags$li(paste0("OLS 估計值 (", round(ols_coef, 3), ") 與真實 LATE 的差異反映了『遵從行為』與結果之間的混淆。在此模擬中，基線 SAQ 可能同時影響治療選擇和結果，導致 OLS 偏誤。")),
              tags$li(paste0("ITT 估計值 (", round(itt_coef, 3), ") 代表治療『策略』的平均效果，其大小受遵從者比例影響。理論上，ITT ≈ LATE × (P(T=1|Z=1) - P(T=1|Z=0))。"))
          ),
          hr(),
          p("注意：以上比較基於假設基線 SAQ 效應為線性的模型。您可以查看包含 RCS 的模型結果，以及新增的 OLS 非線性檢定結果，以評估非線性假設是否影響估計。同時，請檢查 IV 模型的診斷檢定 (如弱工具變數檢定)。") # Updated note
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
