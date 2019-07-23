library(shiny)
library(scales)
library(ggplot2)
library(plotly)
# library(seewave)
library(plot.matrix)

options(shiny.maxRequestSize = 100*1024^2)



# Function to create sine waves famyly  

comSin <- function(Wtime,freq){
  f1 = rep(freq,each = length(Wtime))
  Csin <- exp( 1i*2*pi*f1*Wtime )
  Csin1 <- matrix(Csin, nrow = length(Wtime), byrow = F)
  return(Csin1)
}

# Function to create Gausian windows. 

comGau <- function(Wtime,s){
  s1 = rep(s,each = length(Wtime))
  gaus_win  = exp( (-Wtime^2) / (2*s1^2) )
  gaus_win1 <- matrix(gaus_win, nrow = length(Wtime), byrow = F)
  return(gaus_win1)
}
  
# Set frequency and sine parameters for wavelets

min_frex = 4
max_frex = 40
num_frex = 20

freq = seq(min_frex, max_frex, length.out = num_frex)

# g1 and g2 are the parameters
g1 = 20
g2 = 25
gW   = seq(g1,g2, length.out = length(freq))
s = gW / (2*pi*freq); 


  #####################################_________UI____________###################################

ui <- navbarPage("BackAv",
                 
                 
                 #####################_______GENERAL_VIEW_____####################
                 tabPanel("General view", 
                          
                          fluidRow(column = 12,
                                   plotlyOutput("plot"), 
                                   hr(),         
                                   fluidRow(
                                     column( width = 2, offset = 2,
                                             # channel inputs
                                             numericInput(inputId = "EEG_channel",
                                                          label   = "EEG channel",
                                                          value   = 1  ),
                                             numericInput(inputId = "EMG_channel",
                                                          label   = "EMG channel",
                                                          value   = 1 )),
                                     column( width = 5, offset = 2,
                                             # sampling rate
                                             numericInput(inputId = "srate",
                                                          label   = "Sampling Rate",
                                                          value   = 1000 ),
                            
                              
                                           # sampling rate
                                            checkboxInput(inputId = "rectify ", 
                                                          label   = "rectify EMG channels",
                                                          value   = T  ),
                                            
                                            br(),
                                            fileInput(inputId = "file1", 
                                                      label = "Choose File",
                                                      accept = c(
                                                        "text/plain",
                                                        "text/csv",
                                                        "text/comma-separated-values,text/plain",
                                                        ".csv"))
                                     )
                                   
                                     )
                                     
                                   )),
                                   
                                   
                          
                 ###########################_________TIME_DOMAIN______#####################
                 tabPanel("Time domain",
                          
                          fluidRow(column = 12,
                                   
                                   plotlyOutput("plot_time"), 
                                   br(),
                                   br(),
                                   hr(),                  
                                   fluidRow(
                                            column(width = 4, offset = 1,
                                                   h4("Parameters for EMG markers"),
                                                   numericInput(inputId = "threshold",
                                                                label   = "Threshold",
                                                                step = 0.01,
                                                                value   = 0.2 ),
                                                   numericInput(inputId = "time_after",
                                                                label   = "Time After (in sec)",
                                                                step = 0.01,
                                                                value   = 0.02 ),
                                                   numericInput(inputId = "time_before",
                                                                label   = "Time Before (in sec)",
                                                                step = 0.01,
                                                                value   = 0.02)
                                                   ),
                                             column(width = 3, 
                                                    br(),
                                                    br(),
                                                   numericInput(inputId = "after_a",
                                                                label   = "Amp After > than",
                                                                step = 0.01,
                                                                value   = 0.15 ),
                                                                
                                                   numericInput(inputId = "before_a",
                                                                label   = "Amp Before < than",
                                                                step = 0.01,
                                                                value   = 0.20 )
                                            ),
                                            column(width = 2, 
                                                    br(),
                                                    br(),
                                                    numericInput(inputId = "window",
                                                                 label   = "Window (in sec) ",
                                                                 step = 0.01,
                                                                 value   = 1 ),
                                                    numericInput(inputId = "onset",
                                                                 label   = "Onset (in sec)",
                                                                 step = 0.01,
                                                                 value   = .8 ),
                                                   actionButton("RUN1", "RUN")
                                                   )
                                            )
                                   ) 
                          
                 ),
                 ###########################_________Average______#####################
                 tabPanel("Average",
                          column(12,
                                
                                plotlyOutput("plot_ave"), 
                                br(),
                                hr(),
                                numericInput(inputId = "baseline",
                                             label   = "baseline correction in sec",
                                             step = 0.01,
                                             value   = .1 ),
                                actionButton("RUN2","RUN")
                                
                          )
                 ),
                 
                 
                 ###########################_________Spectrogram______#####################
                tabPanel("Spectrogram",
                         column(12,
                              plotlyOutput("pspec"),
                              # plotlyOutput("pow"),
                              # br(),
                              # hr(),
                              actionButton("RUN3","RUN")
                              # sliderInput("frexpow", "Power at freqency", min = 4, max = 40, value = 4, step = 1.9)
                         )
                         )
                # ,
                # 
                # ###################test########################
                # 
                # tabPanel("test",
                #          column(12,
                #               plotOutput("wav")  ,
                #               tableOutput("test1")
                #          )
                # )
                
)

#####################################________Server___________###################################

server <- function(input, output) {
  
  #### load data
  data <- reactive({
    req(input$file1)
    inFile    <- input$file1
    read.delim(inFile$datapath, header = F, sep = "")
  })
  
  ##### General information
  ext             <- reactive({dim(data())[1]}) # extention of data in pnts
  srate           <- reactive({input$srate})   # to be specify
  tsec            <- reactive({ext()/srate()}) # Extention in seconds
  ttime           <- reactive({seq(1/srate(),tsec(), by = 1/srate())}) #time vector
  threshold       <- reactive({input$threshold})
  time_after      <- reactive({input$time_after * srate()})
  time_before     <- reactive({input$time_before * srate()})
  after_a         <- reactive({input$after_a})
  before_a        <- reactive({input$before_a})
  wind            <- reactive({input$window * srate()})  # window for the average
  onset           <- reactive({input$onset * srate()})   # onset to put the 0   
  timeW           <- reactive({seq(onset()*-1,wind() - onset())}) # time window for the average
  
   # creating wavelets family for the spectrogram 
  
  Wtime           <- reactive(seq(-2,2,by = 1/srate()))  # wavelet time
  half_wav        <- reactive((length(Wtime()) - 1)/2)
  nKern           <- reactive(length(Wtime()));
  
 
  sin_f           <- reactive(t(comSin(Wtime(),freq)))
  gau_f           <- reactive(t(comGau(Wtime(),s)))
  wavelets        <- reactive(sin_f() * gau_f())
  
  zpadD           <- reactive(as.vector( matrix(0L, nrow = 1, ncol = nKern() - 1 ))) # zero padding for data
  

  
  
  
  # select ACC and EMG channels
  EEGch <- reactive({input$EEG_channel})
  EMGch <- reactive({input$EMG_channel})
  
  EEG   <- reactive({data()[,EEGch()]})
  EMG   <- reactive({data()[,EMGch()]})
  
  # 
  t1 <- reactive({plot_ly(x = ttime(), y = EEG()) %>% add_lines( name = "EEG", line = list(width = .5)) %>% rangeslider()  })
  t2 <- reactive({plot_ly(x = ttime(), y = EMG()) %>% add_lines( name = "EMG", line = list(width = .5))  })
  
  output$plot <- renderPlotly({subplot(t1(),t2(),nrows = 2, shareX = T, shareY = T, titleX = T) })
  
  #####Time Domains#####
  # creating markers
  
  datM <- eventReactive(input$RUN1,{
    
  emg     <- scales::rescale(abs(EMG()))
  
  emg1    <- emg > threshold() 
  emgtrl  <- diff(emg1)
  emgon   <- which(emgtrl %in%  1)  # beggening of the burst
  emgoff  <- which(emgtrl %in% -1) #finish of the burst
  ons     <- emgon[which(diff(emgon) > (srate()/3))] # for real data has to be >40
  
  on_after  <- 0;
  on_before <- 0;
  
  for (li in 1:length(ons)) {
    if ((ons[li] - 20) < 1) next
    on_after[li]  <- mean(emg[ons[li]:(ons[li] + time_after())]) # should adjust according to sampling rate?
    on_before[li] <- mean(emg[(ons[li] - time_before()):ons[li]])
  }
  
  ons2 <- ons[which(on_after > after_a() & on_before < before_a() )] # original numbers are 0.04 and 0.03
  beg    <- ons2 - onset()
  ending <- ons2 + (wind() - onset())
  trial <- data.frame( beg[which(beg > 0)], ending[which( beg > 0)])
  trial
  })
  
  p <- reactive({plot_ly(x = ttime(), y = scales::rescale(abs(EMG()))) %>% add_lines( name = "EMG", line = list(width = .5)) %>% rangeslider()   })

  # initiate a line shape object (for each marker)
 
  lines <- reactive({
  
    line <- list(
    type = "line",
    line = list(color = "red"),
    xref = "x",
    yref = "y"
  
  )
    
  l1 <- list()

  for (i in 1:dim(datM())[1]) {
    line[["y0"]] <- 0
    line[["y1"]] <- 0.4
    v1 <- datM()[i,1]
    line[c("x0", "x1")] <- v1/srate()
    l1 <- c(l1, list(line))
    }
  
  l1
  }
  

  )

   p1 <- reactive(layout(p(), title = 'EMG plus markers', shapes = lines()))
   
  output$plot_time <- renderPlotly({p1()})
  
  ## segmentation ##
  
  segDat <- eventReactive(input$RUN2,{
    
    s1 <- matrix(nrow = dim(datM())[1], ncol = wind() + 1)
    s2 <- matrix(nrow = dim(datM())[1], ncol = wind() + 1)
    
    for (ti in 1:dim(datM())[1]) {
      
      s1[ti,] <- EEG()[datM()[ti,1]:datM()[ti,2]]
      s2[ti,] <- EMG()[datM()[ti,1]:datM()[ti,2]]
      
    }
    s3 <- list(s1,s2)
    s3
  })
      
    # baseline correction
  
  EEGpre <- reactive({
    dat <- segDat()[[1]]
  for (ti in 1:dim(datM())[1]) {
    d1 <- mean(dat[ti, 1:(input$baseline*srate())])
    dat[ti,] <- dat[ti,] - d1
  }
    dat
  })
  
 EEGave <- reactive({colMeans(EEGpre())})
 EMGave <- reactive({colMeans(segDat()[[2]])})

 EEGp   <- reactive({plot_ly(x = timeW(), y = EEGave()) %>% add_lines( name = "EEG", line = list(width = .8)) %>%  layout(yaxis = list(autorange = "reversed"))  })
 EMGp   <- reactive({plot_ly(x = timeW(), y = EMGave()) %>% add_lines( name = "EMG", line = list(width = .8)) %>%  layout(yaxis = list(autorange = "reversed"))}) 


 output$plot_ave <- renderPlotly({subplot(EEGp(),EMGp(), nrows = 2,shareX = T ,shareY = F)})

 #  #### Spectrogram #####
 
  # s1  <- eventReactive(input$RUN3,{segDat()[[1]]})
  # len <- reactive(dim(segDat()[[1]])[1]*dim(segDat()[[1]])[2]) # length of apended segments
  # s2  <- reactive(array(s1(),dim(c(len()[1],1) ))  )           # reshape to long segment
  # 
  # ps1 <- reactive(spectro(s2(),f = srate(),plot = F)) # spectrum of new segment
  # 
  # m1  <- reactive((ps1()[["amp"]])) # extract matrix
  # lp  <- reactive(length(ps1()[["freq"]])) # number of frequencies
  # ns  <- reactive(dim(segDat()[[1]])[1])
  # m2  <- reactive(array(m1(),    c( dim(segDat()[[1]])[1],dim(segDat()[[1]])[2], lp() )))  # reshape matrix
  # m3  <- reactive(apply(m2(),c(2,3),mean))
  # 
  # spct  <- eventReactive(input$RUN3,{
  #   lspec <- array(NA, dim = c(256,length(timeW()),dim(segDat()[[1]])[1]))
  #                      for (ti in 1:dim(segDat()[[1]])[1]) {
  #                        d1 <- segDat()[[1]]
  #                        d2 <- d1[ti,]
  #                        s1 <- spectro(d2,f = srate(),plot = F)
  #                        s2 <- s1[["amp"]]
  #                        f1 <- s1[["freq"]]
  #                        t1 <- s1[["time"]]
  #                        lspec[,,ti] <- s2
  #                      }
  #   final <- t((apply(lspec,c(1,2),mean)))
  #   final
  #   })
 
 nSeg    <- reactive(dim(datM())[1]) # number of segments
 nApp    <- reactive(nSeg() * length(timeW()) ) # length of apended segments
 
 
 
 nConv   <- reactive(nApp() + nKern() - 1)
 zpad    <- reactive(as.vector( matrix(0L, nrow = 1, ncol =   nApp() - 1 ))) # zero padding for wavelet
 
 data2  <- reactive({
   dat  <- segDat()[[1]]
   ndat <- matrix(0, nrow = nSeg(),ncol = length(timeW()) )
   
   for (ti in 1:nSeg()){
     ndat[ti,] <- dat[ti,] - EEGave()
     
   }
   ndat
   })
 
 
 
 
 spct  <- eventReactive(input$RUN3,{
 
      conMat  <- matrix(0,nrow = length(freq), ncol = nApp())     
      s1      <- data2()
      s2      <- (array(s1, dim = nApp()  ))           # reshape to long segment
      data1   <- append(zpadD(),s2)
      dataX   <- fft(data1)
   

    for (li in 1:length(freq)) {
      cmw1      <- append(zpad(),wavelets()[li,])
      tempX     <- fft(cmw1)
      cmwX      <- tempX / max(abs(tempX))

      as1       <- fft((dataX*cmwX),inverse = T)
      as2       <- as1[(half_wav() - 1 ):(nConv() - half_wav() - 2)]

      conMat[li,] <- as2

   }
   
     mat2 <- array(abs(conMat),c(num_frex,nSeg(),length(timeW()))) # reshaping
     mat3 <- apply(mat2,c(1,3),mean) # averaging trials
     bas  <- apply(mat3[,100:200],2,mean) #baseline normalization
     
     mat4 <- mat3 /bas
     mat4
 })
 
   f2p <- reactive({
    rfreq <- round(freq,1)
    indF  <- match(input$frexpow, rfreq)
    indF
    })
  # 
  output$pspec <- renderPlotly(plot_ly(x = timeW()[100:(max(timeW()) - 100)], y = freq, z = (spct()[,100:(length(timeW()) - 100)]), type = "contour",contours = list(showlines = F), ncontours = 50))

  # Wave to plot

  w2p <- reactive({
    smooth(spct()[f2p(),])
  })
  t2p <- reactive({
    timeW()[100:(max(timeW()) - 100)]
  })
  #
  output$pow  <- renderPlotly(plot_ly( x = timeW(),y = spct()[f2p(),]) %>% add_lines( name = "Power", line = list(width = .8)))

  
  
  output$test1  <- renderTable({(dim(spct()))})
  #output$linesN <- renderText({input$baseline})
  output$wav    <- renderPlot({plot(   (Re(wavelets()[50,]))  , type="l")})
   
}


##########################_______RUN________################################
shinyApp(ui, server)