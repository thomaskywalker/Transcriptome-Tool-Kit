library(shiny)
library(shinyWidgets)

ui <- fluidPage(
  titlePanel("YAML Configuration for DESeq2 Analysis"),

  tabsetPanel(
    # Data Parameters Main Tab
    tabPanel("Data Parameters",
      tabsetPanel(
        tabPanel("Basic Settings",
          sidebarLayout(
            sidebarPanel(
              textInput("EXP_DESIGN_FORMULA", "Experimental Design Formula", "~group"),
              textInput("REF_GROUP", "Reference Group", "Normal"),
              selectInput("NORMALIZED_FUNCTION", "Normalization Function", choices = c("vst", "rln")),
              numericInput("pvalueThreshold", "P-value Threshold", 0.05),
              numericInput("lfcThreshold", "Log Fold Change Threshold", 1),
              # Add other basic data_parameters fields here
            ),
            mainPanel()
          )
        ),
        tabPanel("Enrichment Analyses",
          sidebarLayout(
            sidebarPanel(
              checkboxInput("enrichment_ctrl", "Control", TRUE),
              checkboxInput("GO_CC", "GO Cellular Component", TRUE),
              checkboxInput("GO_MF", "GO Molecular Function", TRUE),
              # Add other enrichment analyses fields here
            ),
            mainPanel()
          )
        )
        # Add more nested tabs for other lists under data_parameters if needed
      )
    ),
    # Plot Parameters Main Tab
    tabPanel("Plot Parameters",
      tabsetPanel(
        tabPanel("Expression Heatmap",
          sidebarLayout(
            sidebarPanel(
              numericInput("heatmap_gene_num", "Number of Genes for Heatmap", 30),
              selectInput("heatmap_scale", "Heatmap Scale", choices = c("row", "none")),
              # Add other expression_heatmap fields here
            ),
            mainPanel()
          )
        )
        # Add more nested tabs for other lists under plot_parameters if needed
      )
    ),
    # Flow Controller Main Tab
    tabPanel("Flow Controller",
      sidebarLayout(
        sidebarPanel(
          checkboxInput("ctrl", "Control", TRUE),
          checkboxInput("dds_transformation", "DDS Transformation", TRUE),
          checkboxInput("normalized_data_generation", "Normalized Data Generation", TRUE),
          # Add other flow_controller fields here
        ),
        mainPanel()
      )
    )
    # Add more main tabs as needed
  )
)


server <- function(input, output, session) {
    observeEvent(input$runButton, {
    # 当按钮被点击时执行的代码
    # 例如，更新一个文本输出来显示一些信息
    output$text <- renderText({
      paste("按钮已点击", input$runButton, "次。")
    })
  })
  observeEvent(input$saveButton, {
    yaml_list <- list(
      data_parameters = list(
        EXP_DESIGN_FORMULA = input$EXP_DESIGN_FORMULA,
        REF_GROUP = list(group = input$REF_GROUP),
        NORMALIZED_FUNCTION = input$NORMALIZED_FUNCTION,
        pvalueThreshold = as.numeric(input$pvalueThreshold),
        lfcThreshold = as.numeric(input$lfcThreshold),
        enrichment_analyses = list(
          Ctrl = input$enrichment_ctrl,
          GO_CC = input$GO_CC,
          GO_MF = input$GO_MF
          # Add other enrichment analyses inputs
        )
        # Add other data_parameters
      ),
      plot_parameters = list(
        expression_heatmap = list(
          heatmap_gene_num = input$heatmap_gene_num,
          heatmap_scale = input$heatmap_scale
          # Add other expression_heatmap inputs
        )
        # Add other plot_parameters
      ),
      flow_controller = list(
        ctrl = input$ctrl,
        dds_transformation = input$dds_transformation,
        normalized_data_generation = input$normalized_data_generation
        # Add other flow_controller fields
      )
      # Add additional sections as needed
    )
    
    # Convert to YAML and handle the output
    yaml_str <- yaml::as.yaml(yaml_list)
    # Save or display the YAML string
    session$onSessionEnded(function() {
    # 當用戶關閉窗口時停止app
    stopApp()
     })
  })
}



shinyApp(ui = ui, server = server)
