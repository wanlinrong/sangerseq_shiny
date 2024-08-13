#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#

################################################
# setRepositories(ind = 1:7)                   #
# options(rsconnect.packrat = TRUE)            #
# rsconnect::deployApp('../sangerseq_shiny/')  #
################################################


library(shiny)
library(tidyverse)
library(Biostrings)
library(shinytoastr)
library(bslib)





# Define UI for application that draws a histogram
ui <- fluidPage(#theme = shinythemes::shinytheme("simplex"),
                #shinythemes::themeSelector(),
                title = "wanlinrong:sangerseq_shiny",
                waiter::autoWaiter(),
                shinytoastr::useToastr(),#设置消息提示
    # Application title
    titlePanel("Analysis and separation of sequencing files[*.ab1]"),
    
    p(em("Maintained by Linrong Wan & Bowen Yin [3S team]")),
    
    # Sidebar with a slider input for number of bins 
    sidebarLayout(
      sidebarPanel(width = 3,
        p('Select the ab1 file that needs to be parsed, allowing multiple files to be selected(required).'),
        
        fileInput('ab1_file_path','select your ab1 files',accept = c('.ab1')),
     
        #tableOutput('file_selected'),
        hr(),
        p("In general, the sequencing results will have certain noise data, which can be corrected by setting the signal-to-noise ratio parameters."),
        numericInput('signal_noise_ratio', label= 'signal / noise ratio', value= 0.33,min = 0.1, max = 1.0, step = 0.1),
        hr(),
        p("In general, the front[trim5] and back[trim3] parts of the sequencing results are inaccurate,and it is best to discard them."),
        numericInput('ab1_trim5', label= "ab1_trim5'", value= 50,min = 1, max = 1000, step = 1),
        numericInput('ab1_trim3', label= "ab1_trim3'", value= 20,min = 1, max = 1000, step = 1),
        hr(),
        p('The reference sequence is used to correct the sequencing results and analyze the results of the variation (required)'),
        textAreaInput('ref_seq', 'Reference seq', value = 'paste seq here',placeholder =T),
        verbatimTextOutput("ref",placeholder = T),
        hr(),
        actionButton(inputId = 'run','submit your mission'),
        ),
        
#-----------------------------------------------------------------------------------------------------------------------      
      # Show a plot of the generated distribution
        mainPanel(width = 9,
                  navset_tab(
                   nav_panel(title = "Chromatograms_Plot", plotOutput('Chromatograms')),
                   nav_panel(title = "Aligemnt_Method",
                             verbatimTextOutput("aligemnt_result")),
                   nav_panel(title = "Matching_Method",
                             verbatimTextOutput("seq_matching"))
                  ),
           
           #msaROutput("aligemnt_result"),
          
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output,session) {

  output$file_selected <- renderTable(colnames = F,bordered = T,
    {input$ab1_file_path$name }) #输出已选择文件
 
  ref_seq <- reactive({input$ref_seq %>% str_extract_all('[ATCG]+') %>% unlist() %>% paste0(collapse = '') })
  
  output$ref <- renderText({   ref_seq()   })
 
  
  #Chromatograms_Plot
  output$Chromatograms <- renderPlot({
    req(input$ab1_file_path$datapath,input$ab1_trim5,input$ab1_trim3,input$signal_noise_ratio)
    sangerseqR::readsangerseq(input$ab1_file_path$datapath) %>% 
      makeBaseCalls(ratio = input$signal_noise_ratio ) %>% 
      chromatogram(width = 100,
                   height = 2, 
                   trim5 = input$ab1_trim5, 
                   trim3 = input$ab1_trim3,
                   showcalls = "both")  })   #Chromatograms_Plot
  
  #Aligemnt_Method
  observeEvent(input$run,{
    req(input$ab1_file_path,input$signal_noise_ratio,input$ab1_trim5,input$ab1_trim3,input$ref_seq)
    
    seq <- core_sanggerseq_split_new(ab1 = input$ab1_file_path$datapath,
                               signal_noise_ratio = input$signal_noise_ratio,
                               ab1_trim5 = input$ab1_trim5,
                               ab1_trim3 = input$ab1_trim3,
                               ref_seq = ref_seq()) %>% 
      tryCatch(error = function(e) shinytoastr::toastr_error(title = "aligemnt_error", conditionMessage(e),showDuration = 10000),
               warning = function(w) shinytoastr::toastr_warning(title = "aligemnt_warning", conditionMessage(w),showDuration = 10000),
               finally = shinytoastr::toastr_success('aligemnt_完成'))
                          
    #Sys.sleep(2)
    # aln <- msa::msaClustalW(p,type ="dna" ) %>% msa::msaConvert(type = "ape::DNAbin")
    # 
    # output$aligemnt_result <- renderMsaR( msaR(aln,rowheight = 25,alignmentHeight = 500,) )
    
    
    output$aligemnt_result <- renderPrint({ 
      
      req(seq,secondarySeq(seq))
      
      pwalign::pairwiseAlignment(ref_seq(),secondarySeq(seq),gapExtension = 6.66) %>% writePairwiseAlignments(block.width = 150) 
      
      })
    
    }) #Aligemnt_Method
  
  
  #Matching_Method
  observeEvent(input$run,{
    req(input$ab1_file_path,input$signal_noise_ratio,input$ab1_trim5,input$ab1_trim3,input$ref_seq)
    
    seq <- core_sanggerseq_split(ab1 = input$ab1_file_path$datapath,
                                     signal_noise_ratio = input$signal_noise_ratio,
                                     ab1_trim5 = input$ab1_trim5,
                                     ab1_trim3 = input$ab1_trim3,
                                     ref_seq = ref_seq()) %>% 
      tryCatch(error = function(e) shinytoastr::toastr_error(title = "matching_error", conditionMessage(e),showDuration = 10000),
               warning = function(w) shinytoastr::toastr_warning(title = "matching_warning", conditionMessage(w),showDuration = 10000),
               finally = shinytoastr::toastr_success('matching_完成'))
    
    #Sys.sleep(2)
    # aln <- msa::msaClustalW(p,type ="dna" ) %>% msa::msaConvert(type = "ape::DNAbin")
    # 
    # output$aligemnt_result <- renderMsaR( msaR(aln,rowheight = 25,alignmentHeight = 500,) )
    
    
    output$seq_matching <- renderPrint({ 
      
      req(seq)
      
      pwalign::pairwiseAlignment(ref_seq(),seq,gapExtension = 6.66) %>% writePairwiseAlignments(block.width = 150) 
      
    })
    
  }) #Matching_Method

  
}

# Run the application 
shinyApp(ui = ui, server = server)

























