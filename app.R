## app.R ##
library(shiny)
library(shinydashboard)
library(xlsx)
library("DT")
library(limma)
library(DESeq2)
library(edgeR)
library(caret)
library(mclust)
library(geneplotter)
library(cluster)
library(fossil)
library(SummarizedExperiment)
library(xlsx)
library(clValid)
library(genefilter)
library(survival)
library(lubridate)
library(foreach)
library(PoiClaClu)
library(dplyr)
library(NbClust)
library(MBCluster.Seq)
library(reticulate)
library(pheatmap)
library(plotly)
library(heatmaply)
library(ggfortify)
library(FactoMineR)
options(shiny.maxRequestSize = 30*1024^2)

use_python('C:/Users/Lenovo/AppData/Local/Programs/Python/Python36/python.exe')
sns <- import('seaborn')
plt <- import('matplotlib.pyplot')
pd <- import('pandas')

######## FUNCTIONS ########
NullModelAhu = function (x, type = c("mle", "deseq", "quantile", "none")) 
{
  type <- match.arg(type)
  rowsum <- rowSums(x)
  colsum <- colSums(x)
  mle <- outer(rowsum, colsum, "*")/sum(rowsum)
  if (type == "mle") {
    return(list(n = mle, sizes = rowSums(x)/sum(x)))
  }
  else if (type == "quantile") {
    sample.qts <- apply(x, 1, quantile, 0.75)
    sample.qts <- pmax(sample.qts, 1)
    sample.qts <- sample.qts/sum(sample.qts)
    fit <- outer(sample.qts, colsum, "*")
    return(list(n = fit, sizes = sample.qts))
  }
  else if (type == "deseq") {
    counts <- t(x)
    geomeans <- exp(rowMeans(log(counts)))
    sizes <- apply(counts, 2, function(cnts) median((cnts/geomeans)[geomeans > 
                                                                      0]))
    rawsizestr <- sizes
    sizes <- sizes/sum(sizes)
    fit <- outer(sizes, colsum, "*")
    return(list(n = fit, sizes = sizes, geomeans = geomeans, 
                rawsizestr = rawsizestr))
  }
  else if (type == "none") {
    return(list(n = x))
  }
}

PoissonDistanceAhu = function (x, beta = 1, type = c("mle", "deseq", "quantile", "none"), 
                               transform = TRUE, alpha = NULL, perfeature = FALSE) 
{
  type <- match.arg(type)
  if (!transform && !is.null(alpha)) 
    stop("You have asked for NO transformation but have entered alpha.")
  if (transform && !is.null(alpha)) {
    if (alpha > 0 && alpha <= 1) 
      x <- x^alpha
    if (alpha <= 0 || alpha > 1) 
      stop("alpha must be between 0 and 1")
  }
  if (transform && is.null(alpha)) {
    alpha <- FindBestTransform(x)
    x <- x^alpha
  }
  dd <- matrix(0, nrow = nrow(x), ncol = nrow(x))
  ddd <- NULL
  if (perfeature) 
    ddd <- array(0, dim = c(nrow(x), nrow(x), ncol(x)))
  for (i in 2:nrow(dd)) {
    xi <- x[i, ]
    for (j in 1:(i - 1)) {
      xj <- x[j, ]
      n <- NullModelAhu(x[c(i, j), ], type = type)$n
      ni <- n[1, ]
      nj <- n[2, ]
      di <- (xi + beta)/(ni + beta)
      dj <- (xj + beta)/(nj + beta)
      dd[i, j] <- sum(ni + nj - ni * di - nj * dj + xi * 
                        log(di) + xj * log(dj))
      if (perfeature) 
        ddd[j, i, ] <- ddd[i, j, ] <- ni + nj - ni * 
        di - nj * dj + xi * log(di) + xj * log(dj)
    }
  }
  save <- list(dd = as.dist(dd + t(dd)), alpha = alpha, x = x, 
               ddd = ddd, alpha = alpha, type = type)
  class(save) <- "poidist"
  return(save)
}


###############################  DISTANCE MATRIX FUNCTIONS  ############################################
weighted.distance = function (datam, method = c("euclidean","squared euclidean",
                                                "manhattan","canberra","chebychev",
                                                "bray curtis","cosine correlation",
                                                "pearson correlation","uncentered pearson correlation"),
                              cluster = c("row","column"))
{
  
  
  E=voomWithQualityWeights(datam, normalize="none", plot=FALSE)$E
  w=voomWithQualityWeights(datam, normalize="none", plot=FALSE)$target[,2] #weights ##sample.weights veya target[,2]
  
  
  if (cluster == "row")
  {E = E}
  else if (cluster == "column")
  {E = t(E)
  w = t(w)}
  
  x <- list()
  #if (cluster == "row")
  #{
  for (i in 1:nrow(E)){
    x[[i]] <- E[i,]
  }
  #}
  #else if (cluster == "column")
  #{for (i in 1:ncol(E)){
  # x[[i]] <- E[,i]
  # }}
  
  y <- NULL
  y <- matrix(nrow=nrow(E),ncol=nrow(E))
  
  #foreach(i=1:nrow(E)) %:%
  for (i in 1:nrow(E)){ 
    #foreach(j=(i+1):nrow(E)) %dopar% {
    y[i,i] <- 0
    for (j in i+1:nrow(E)){
      if (i<=nrow(E) && j<=nrow(E)){
        if (method == "euclidean")
        {
          y[i,j] <- sqrt(sum(w[i]*(x[[i]]-x[[j]])^2))
          y[j,i] <- y[i,j]
        }
        else if (method == "squared euclidean")
        {
          y[i,j] <- sum(w[i]*(x[[i]]-x[[j]])^2)
          y[j,i] <- y[i,j]
        }
        else if (method == "manhattan")
        {
          y[i,j] <- sum(w[i]*abs(x[[i]]-x[[j]]))
          y[j,i] <- y[i,j]
        }
        else if (method == "canberra")
        {
          y[i,j] <- sum(w[i]*(abs(x[[i]]-x[[j]])/abs(x[[i]]+x[[j]])))
          y[j,i] <- y[i,j]
        }
        else if (method == "chebychev")
        {
          y[i,j] <- max(w[i]*abs(x[[i]]-x[[j]]))
          y[j,i] <- y[i,j]
        }
        else if (method == "bray curtis")
        {
          #y[i,j] <- sum(w*abs(x[[i]]-x[[j]]))/sum(w*x[[i]]+x[[j]])
          y[i,j] <- sum(w[i]*abs(x[[i]]-x[[j]]))/sum(w[i]*(x[[i]]+x[[j]]))
          y[j,i] <- y[i,j]
        }
        else if (method == "cosine correlation")
        {
          #y[i,j] <- sum(w*(x[[i]]*x[[j]]))/sqrt(sum(x[[i]])*sum(x[[j]]))
          y[i,j] <- sum(w[i]*(x[[i]]*x[[j]]))/(sqrt(sum(w[i]*x[[i]]^2))*sqrt(sum(w[i]*x[[j]]^2)))
          y[j,i] <- y[i,j]
        }
        else if (method == "pearson correlation")
        {
          #xi_ort = sum(w*x[[i]])/sum(w)
          #xj_ort = sum(w*x[[j]])/sum(w)
          
          #sxj = sum(w*(x[[j]]-xj_ort)^2)/sum(w)
          #sxixj = sum(w*(x[[i]]-xi_ort)*(x[[j]]-xj_ort))
          #y[i,j] <- sxixj / sqrt(sxi*sxj)
          #y[j,i] <- y[i,j]
          y[i,j] <- sum(w[i]*(x[[i]]-mean(x[[i]]))*(x[[j]]-mean(x[[j]])))/(sqrt(sum(w[i]*(x[[i]]-mean(x[[i]]))^2))*sqrt(sum(w[i]*(x[[j]]-mean(x[[j]]))^2)))
          y[j,i] <- y[i,j]
        }
        else if (method == "uncentered pearson correlation")
        {
          #xi_ort = sum(w*x[[i]])/sum(w)
          #xj_ort = sum(w*x[[j]])/sum(w)
          #sxi = sum(w*(x[[i]]-xi_ort)^2)/sum(w)
          #sxj = sum(w*(x[[j]]-xj_ort)^2)/sum(w)
          #sxixj = sum(w*(x[[i]]*x[[j]]))
          #y[i,j] <- sxixj / sqrt(sxi*sxj)
          #y[j,i] <- y[i,j]
          y[i,j]  <- sum(w[i]*x[[i]]*x[[j]])/(sqrt(sum(w[i]*(x[[i]]-mean(x[[i]]))^2))*sqrt(sum(w[i]*(x[[j]]-mean(x[[j]]))^2)))
          y[j,i] <- y[i,j]
        }
      }
    }
  }
  
  if (method == "cosine correlation" || method == "pearson correlation" || method == "uncentered pearson correlation")
  {y=1-y}
  return(y)
}

weighted.distance2 = function (datam, method = c("euclidean","squared euclidean",
                                                 "manhattan","canberra","chebychev",
                                                 "bray curtis","cosine correlation",
                                                 "pearson correlation","uncentered pearson correlation"),
                               cluster = c("row","column"))
{
  
  
  E=voomWithQualityWeights(datam, normalize="none", plot=FALSE)$E
  w=voomWithQualityWeights(datam, normalize="none", plot=FALSE)$weights #weights ##sample.weights veya target[,2]
  
  
  if (cluster == "row")
  {E = E}
  else if (cluster == "column")
  {E = t(E)
  w = t(w)}
  
  x <- list()
  #if (cluster == "row")
  #{
  for (i in 1:nrow(E)){
    x[[i]] <- E[i,]
  }
  #}
  #else if (cluster == "column")
  #{for (i in 1:ncol(E)){
  # x[[i]] <- E[,i]
  # }}
  
  y <- NULL
  y <- matrix(nrow=nrow(E),ncol=nrow(E))
  
  #foreach(i=1:nrow(E)) %:%
  for (i in 1:nrow(E)){ 
    #foreach(j=(i+1):nrow(E)) %dopar% {
    y[i,i] <- 0
    for (j in i+1:nrow(E)){
      if (i<=nrow(E) && j<=nrow(E)){
        if (method == "euclidean")
        {
          y[i,j] <- sqrt(sum(w[i,]*(x[[i]]-x[[j]])^2))
          y[j,i] <- y[i,j]
        }
        else if (method == "squared euclidean")
        {
          y[i,j] <- sum(w[i,]*(x[[i]]-x[[j]])^2)
          y[j,i] <- y[i,j]
        }
        else if (method == "manhattan")
        {
          y[i,j] <- sum(w[i,]*abs(x[[i]]-x[[j]]))
          y[j,i] <- y[i,j]
        }
        else if (method == "canberra")
        {
          y[i,j] <- sum(w[i,]*(abs(x[[i]]-x[[j]])/abs(x[[i]]+x[[j]])))
          y[j,i] <- y[i,j]
        }
        else if (method == "chebychev")
        {
          y[i,j] <- max(w[i,]*abs(x[[i]]-x[[j]]))
          y[j,i] <- y[i,j]
        }
        else if (method == "bray curtis")
        {
          #y[i,j] <- sum(w*abs(x[[i]]-x[[j]]))/sum(w*x[[i]]+x[[j]])
          y[i,j] <- sum(w[i,]*abs(x[[i]]-x[[j]]))/sum(w[i,]*(x[[i]]+x[[j]]))
          y[j,i] <- y[i,j]
        }
        else if (method == "cosine correlation")
        {
          #y[i,j] <- sum(w*(x[[i]]*x[[j]]))/sqrt(sum(x[[i]])*sum(x[[j]]))
          y[i,j] <- sum(w[i,]*(x[[i]]*x[[j]]))/(sqrt(sum(w[i,]*x[[i]]^2))*sqrt(sum(w[i,]*x[[j]]^2)))
          y[j,i] <- y[i,j]
        }
        else if (method == "pearson correlation")
        {
          #xi_ort = sum(w*x[[i]])/sum(w)
          #xj_ort = sum(w*x[[j]])/sum(w)
          
          #sxj = sum(w*(x[[j]]-xj_ort)^2)/sum(w)
          #sxixj = sum(w*(x[[i]]-xi_ort)*(x[[j]]-xj_ort))
          #y[i,j] <- sxixj / sqrt(sxi*sxj)
          #y[j,i] <- y[i,j]
          y[i,j] <- sum(w[i,]*(x[[i]]-mean(x[[i]]))*(x[[j]]-mean(x[[j]])))/(sqrt(sum(w[i,]*(x[[i]]-mean(x[[i]]))^2))*sqrt(sum(w[i,]*(x[[j]]-mean(x[[j]]))^2)))
          y[j,i] <- y[i,j]
        }
        else if (method == "uncentered pearson correlation")
        {
          #xi_ort = sum(w*x[[i]])/sum(w)
          #xj_ort = sum(w*x[[j]])/sum(w)
          #sxi = sum(w*(x[[i]]-xi_ort)^2)/sum(w)
          #sxj = sum(w*(x[[j]]-xj_ort)^2)/sum(w)
          #sxixj = sum(w*(x[[i]]*x[[j]]))
          #y[i,j] <- sxixj / sqrt(sxi*sxj)
          #y[j,i] <- y[i,j]
          y[i,j]  <- sum(w[i,]*x[[i]]*x[[j]])/(sqrt(sum(w[i,]*(x[[i]]-mean(x[[i]]))^2))*sqrt(sum(w[i,]*(x[[j]]-mean(x[[j]]))^2)))
          y[j,i] <- y[i,j]
        }
      }
    }
  }
  if (method == "cosine correlation" || method == "pearson correlation" || method == "uncentered pearson correlation")
  {y=1-y}
  return(y)
}


raw.distance = function (datam, method = c("euclidean","squared euclidean",
                                           "manhattan","canberra","chebychev",
                                           "bray curtis","cosine correlation",
                                           "pearson correlation","uncentered pearson correlation"),
                         cluster = c("row","column"))
{
  if (cluster == "column")
  {datam = t(datam)}
  else if (cluster == "row")
  {datam = datam}
  
  x <- list()
  #if (cluster == "row")
  #{
  for (i in 1:nrow(datam)){
    x[[i]] <- datam[i,]
  }
  #}
  #else if (cluster == "column")
  #{for (i in 1:ncol(datam)){
  # x[[i]] <- datam[,i]
  #}}
  
  y <- NULL
  y <- matrix(nrow=nrow(datam),ncol=nrow(datam))
  
  #foreach(i=1:nrow(datam)) %:%
  for (i in 1:nrow(datam)){ 
    #foreach(j=(i+1):nrow(datam)) %dopar% {
    for (j in i+1:nrow(datam)){
      y[i,i] <- 0
      if (i<=nrow(datam) && j<=nrow(datam)){
        if (method == "euclidean")
        {
          y[i,j] <- sqrt(sum((x[[i]]-x[[j]])^2))
          y[j,i] <- y[i,j]
        }
        else if (method == "squared euclidean")
        {
          y[i,j] <- sum((x[[i]]-x[[j]])^2)
          y[j,i] <- y[i,j]
        }
        else if (method == "manhattan")
        {
          y[i,j] <- sum(abs(x[[i]]-x[[j]]))
          y[j,i] <- y[i,j]
        }
        else if (method == "canberra")
        {
          y[i,j] <- sum((abs(x[[i]]-x[[j]])/abs(x[[i]]+x[[j]])))
          y[j,i] <- y[i,j]
        }
        else if (method == "chebychev")
        {
          y[i,j] <- max(abs(x[[i]]-x[[j]]))
          y[j,i] <- y[i,j]
        }
        else if (method == "bray curtis")
        {
          y[i,j] <- sum(abs(x[[i]]-x[[j]]))/sum(x[[i]]+x[[j]])
          y[j,i] <- y[i,j]
        }
        else if (method == "cosine correlation")
        {
          #y[i,j] <- sum((x[[i]]*x[[j]]))/(sqrt(sum(x[[i]]^k))*sqrt(sum(x[[j]]^k)))
          y[i,j] <- sum(x[[i]]*x[[j]])/(sqrt(sum(x[[i]]^2))*sqrt(sum(x[[j]]^2)))
          y[j,i] <- y[i,j]
        }
        else if (method == "pearson correlation")
        {
          y[i,j] <- sum((x[[i]]-mean(x[[i]]))*(x[[j]]-mean(x[[j]])))/(sqrt(sum((x[[i]]-mean(x[[i]]))^2))*sqrt(sum((x[[j]]-mean(x[[j]]))^2)))
          y[j,i] <- y[i,j]                                         
        }
        else if (method == "uncentered pearson correlation")
        {
          y[i,j]  <- sum(x[[i]]*x[[j]])/(sqrt(sum((x[[i]]-mean(x[[i]]))^2))*sqrt(sum((x[[j]]-mean(x[[j]]))^2)))
          y[j,i] <- y[i,j]
        }
      }
    }
  }
  if (method == "cosine correlation" || method == "pearson correlation" || method == "uncentered pearson correlation")
  {y=1-y}
  return(y)
} 


 ui <- dashboardPage(
  #dashboardHeader(title = "voomQW & voomPW:Clustering of RNA-Seq Data",titleWidth = 450),
  dashboardHeader(
    dropdownMenu(
      type = "notifications", 
      icon = icon("question-circle"),
      badgeStatus = NULL,
      headerText = "See also:",
      
      notificationItem("shiny", icon = icon("file"),
                       href = "http://shiny.rstudio.com/"),
      notificationItem("shinydashboard", icon = icon("file"),
                       href = "https://rstudio.github.io/shinydashboard/"),
      notificationItem(
        text = "5 new users today",
        icon("users")
      ),
      notificationItem(
        text = "12 items delivered",
        icon("truck"),
        status = "success"
      ),
      notificationItem(
        text = "Server load at 86%",
        icon = icon("exclamation-triangle"),
        status = "warning"
      )
    ),
    dropdownMenu(
      type = "messages"
    ),
    dropdownMenu(type = "tasks", badgeStatus = "success",
                 taskItem(value = 90, color = "green",
                          "Documentation"
                 ),
                 taskItem(value = 17, color = "aqua",
                          "Project X"
                 ),
                 taskItem(value = 75, color = "yellow",
                          "Server deployment"
                 ),
                 taskItem(value = 80, color = "red",
                          "Overall project"
                 )
    )
  ),
  dashboardSidebar(
    width = 250,
    sidebarMenu(
      tags$head(tags$style(HTML('.shiny-server-account { display: none; }'))),
      uiOutput("userpanel"),
      menuItem("Introduction", tabName = "Introduction",icon = icon("info-circle")),
      menuItem("Data Upload", tabName = "dataupload", icon = icon("download")),
      menuItem("Clustering", tabName = "clustering",icon = icon("sitemap")),
      menuItem("Heatmaps", tabName = "heatmaps", icon = icon("th")),
      menuItem("PCA", tabName = "PCA", icon = icon("circle")),
      menuItem("Manual", icon = icon ("comment")),
      menuItem("Authors & News",tabName="authors",icon = icon("users")),
      menuItem("Citation", icon = icon ("file-alt"))
    )
  ),
  dashboardBody(
    tags$style(type='text/css',"p { font-size: 18px; }"),
    tags$style(type='text/css',"p { font-family: Helvetica; }"),
    tabItems(
    tabItem(tabName = "dataupload",
      box(
        width = 2, status = "primary",h4("Input Data"),
        radioButtons("dataInput", "", list("Load example data"=1, "Upload a file"=2), selected=1),
        conditionalPanel(condition="input.dataInput=='1'",
                         h5("Datasets:"),
                         radioButtons("sampleData", "", list("Cervical data (n=58, p=714)"=1, "Alzheimer data (n=70, p=2801)"=2), selected=1),
                         h5("n: number of observations"),
                         h5("p: number of variables")
        ),
        
        conditionalPanel(condition="input.dataInput=='2'",
                         h5("Upload a delimited text file: "),
                         radioButtons("radio3", "File Format:", list(".xlsx"=1, ".txt"=2, ".sav"=3, ".csv"=4), selected=1),
                         #HTML('<i class="fa fa-beer fa-lg"></i>'),
                         fileInput("upload", "", multiple = FALSE),
                         radioButtons("fileSepDF", "Delimiter:", list("Comma"=1, "Tab"=2, "Semicolon"=3, "Space"=4),selected=2),
                         
                         conditionalPanel(condition="input.fileSepDF!='1'",
                                          checkboxInput(inputId = "decimal", label = "Use comma as decimal", value = FALSE)
                         ),
                         
                         HTML('<br>'),
                         HTML('<p>You can upload your data as separated by comma, tab, semicolon or space.</p>'),
                         HTML('<p><b>Note</b>: First row must be header.</p>')
        )
      ),
      box(width = 10, status = "primary",h4("Data"),
          div(style = 'overflow-x: scroll', DT::dataTableOutput('RawData')))
    ),
    tabItem(tabName = "authors",
            h4("Authors"),
            br(),
            HTML('<p><a href="http://bit.do/ahudurmuscelebi" target="_blank"> <b>Ahu Durmuscelebi, MSc</b></a><p>'),
            HTML('<p>Erciyes University Faculty of Medicine <a href="http://http://hastaneler.erciyes.edu.tr/Sayfa/Temeltip/3102" target="_blank"> Department of Biostatistics</a><p>'),
            HTML('<p><a href="mailto:ahudurmuscelebi.87@gmail.com" target="_blank">ahudurmuscelebi.87@gmail.com</a><p>'),
            
            br(),
            HTML('<p><a href="http://bit.do/gokmenzararsiz" target="_blank"> <b>Gokmen Zararsiz, PhD</b></a><p>'),
            HTML('<p>Hacettepe University Faculty of Medicine <a href="http://www.biostatistics.hacettepe.edu.tr" target="_blank"> Department of Biostatistics</a><p>'),
            HTML('<p><a href="mailto:gokmen.zararsiz@hacettepe.edu.tr" target="_blank">gokmen.zararsiz@hacettepe.edu.tr</a><p>'),
            
            br(),
            HTML('<p><a href="http://bit.do/gozdeezararsiz" target="_blank"> <b>Gozde E. Zararsiz, PhD</b></a><p>'),
            HTML('<p>Hacettepe University Faculty of Medicine <a href="http://www.biostatistics.hacettepe.edu.tr" target="_blank"> Department of Biostatistics</a><p>'),
            HTML('<p><a href="mailto:gozdeerturk@erciyes.edu.tr " target="_blank">gozdeerturk@erciyes.edu.tr </a><p>'),
            
            br(),
            
            HTML('<p><a href="http://www.biostatistics.hacettepe.edu.tr/cv/Dincer_Goksuluk_CV_Eng.pdf" target="_blank"> <b>Dincer Goksuluk, PhD</b></a><p>'),
            HTML('<p>Hacettepe University Faculty of Medicine <a href="http://www.biostatistics.hacettepe.edu.tr" target="_blank"> Department of Biostatistics</a><p>'),
            HTML('<p><a href="mailto:dincer.goksuluk@hacettepe.edu.tr" target="_blank">dincer.goksuluk@hacettepe.edu.tr</a><p>'),
            br(),
            
            HTML('<p><a href="http://www-huber.embl.de/users/klaus/" target="_blank"> <b>Bernd Klaus, PhD</b></a><p>'),
            HTML('<p>EMBL Heidelberg<p>'),
            HTML('<p><a href="bernd.klaus@embl.de" target="_blank">bernd.klaus@embl.de</a><p>'),
            br(),
            
            HTML('<p><a href="http://yunus.hacettepe.edu.tr/~selcuk.korkmaz/" target="_blank"> <b>Selcuk Korkmaz, PhD</b></a><p>'),
            HTML('<p>Hacettepe University Faculty of Medicine <a href="http://www.biostatistics.hacettepe.edu.tr" target="_blank"> Department of Biostatistics</a><p>'),
            HTML('<p><a href="mailto:selcuk.korkmaz@hacettepe.edu.tr" target="_blank">selcuk.korkmaz@hacettepe.edu.tr</a><p>'),
            br(),
            
            HTML('<p><a href="http://aves.istanbul.edu.tr/vahap.eldem/" target="_blank"> <b>Vahap Eldem, PhD</b></a><p>'),
            HTML('<p>Istanbul University Faculty of Science <a href="http://fen.istanbul.edu.tr/biyoloji/#" target="_blank"> Department of Biology</a><p>'),
            HTML('<p><a href="mailto:vahap.eldem@istanbul.edu.tr" target="_blank">vahap.eldem@istanbul.edu.tr</a><p>'),
            br(),

            HTML('<p><a href="http://aves.erciyes.edu.tr/ahmetozturk/" target="_blank"> <b>Ahmet Ozturk, PhD</b></a><p>'),
            HTML('<p>Erciyes University Faculty of Medicine <a href="http://biyoistatistik.erciyes.edu.tr" target="_blank"> Department of Biostatistics</a><p>'),
            HTML('<p><a href="mailto:ahmets67@hotmail.com" target="_blank">ahmets67@hotmail.com</a><p>'),
            
            
            HTML('<br>'),
            
            
            h4("News"),
            h5("Version 1.5 (November 25, 2016)"),
            HTML('<p>(2) Lung cancer data added as an example dataset <p>'),
            HTML('<p>(2) Bug fixes and improvements <p>'),
            
            h5("Version 1.4 (November 14, 2016)"),
            HTML('<p>(1) Bug fixes and improvements <p>'),
            
            h5("Version 1.3 (November 5, 2016)"),
            HTML('<p>(1) Gene ontology results added <p>'),
            HTML('<p>(2) Gene ontology plot added <p>'),
            HTML('<p>(3) Bug fixes and improvements <p>'),
            
            h5("Version 1.2 (August 20, 2016)"),
            HTML('<p>(1) Heatmap added <p>'),
            HTML('<p>(2) Network plot added <p>'),
            HTML('<p>(3) Bug fixes and improvements <p>'),
            
            h5("Version 1.1 (July 18, 2016)"),
            HTML('<p>(1) Upgraded to shiny version 0.14 <p>'),
            
            h5("Version 1.0 (June 18, 2015)"),
            HTML('<p>(1) VoomDDA web application has been released. <p>'),
            
            
            
            
            
            
            
            
            HTML('<br>'),
            
            h5("Other Tools"),
            
            HTML('<p><a href="http://www.bioconductor.org/packages/release/bioc/html/MLSeq.html" target="_blank"> <b>MLSeq: Machine learning interface for RNA-Seq data </b></a><p>'),
            HTML('<p><a href="http://www.biosoft.hacettepe.edu.tr/MLViS/" target="_blank"> <b>MLViS: machine learning-based virtual screening tool </b></a><p>'),
            HTML('<p><a href="http://www.biosoft.hacettepe.edu.tr/easyROC" target="_blank"> <b>easyROC: a web-tool for ROC curve analysis </b></a><p>'),
            HTML('<p><a href="http://www.biosoft.hacettepe.edu.tr/MVN" target="_blank"> <b>MVN: a web-tool for assessing multivariate normality </b></a><p>'),
            HTML('<p><a href="http://www.biosoft.hacettepe.edu.tr/DDNAA/" target="_blank"> <b>DDNAA: Decision support system for differential diagnosis of nontraumatic acute abdomen </b></a><p>'),
            HTML('<br>'),
            
            
            
            h6("Please feel free to send us bugs and feature requests.")      
    ),
    tabItem(tabName = "Introduction",
            #fluidRow(column(3,
             #               imageOutput("images12")),
                #     column(3,imageOutput("images15")),
                 #    column(3,imageOutput("images13"))
                  #   
            #),
            fluidRow(column(9,
            HTML('<br>'),
            HTML('<p>voomQW and voomPW is a web-tool that is developed for clustering analises for RNA-sequencing data. </p>'),
            HTML('<p>In addition to the user can upload their own data to the program with .xlsx, .txt, .sav, .csv etc. extensions of file, they can use data within the program. 
                Before clustering process, it can be done do various filtering and normalizastion steps for pre-processing. Then, the user 
                do cluster analises with vst,rlog,Poisson, model-based and EdgeR methods that previously recommended or voomQW and voomPW, which newly created.
                Both voomQW and voomPW use voom transformation methods; but voomQW obtains distance matrix by combining log-cpm values at the individuals 
                observations level and their associated quality weights quality weights which are combined with sample-specific weights; while 
                voomPW obtains distance matris by combining log-cpm values and their precisions weights which are calculated from weights at observational
                levels. After clustering, the user can see which observationals in which cluster is estimated with numbers. Moreover, heatmaps and 
                the other charts gives information visually.</p>'),
            HTML('<br>'),            
            HTML('<p>The tool is maked by using several R packages such as shiny, shinydasboard, heatmaply, limma, etc. The source code of the tool is available
               in Github.</p>')
            ),
            column(3,HTML('<br>'),HTML('<br>'),HTML('<br>'),imageOutput("images4"))),
            
            
            fluidRow(
              column(3,imageOutput("images7"))
              ,
              (column(6,
                      HTML('<p> [1] Zararsiz, G., Goksuluk, G., Korkmaz, S., et al. (2015). VoomDDA: Discovery of Diagnostic Biomarkers and Classification of RNA-Seq Data.</p>'),
                      
                      HTML('<p> [2] Law, C.W., Chen, Y., Shi, W., et al. (2014). <a href="http://www.genomebiology.com/2014/15/2/R29">voom: Precision weights unlock linear model analysis tools for RNA-Seq read counts.</a> Genome Biology; 15:R29.</p>'),
                      
                      HTML('<p> [3] Tibshirani, R., Hastie, T., Narasimhan, B., et al. (2002). <a href="http://www.pnas.org/content/99/10/6567.abstract" target="_blank">Diagnosis of multiple cancer types by shrunken centroids of gene expression.</a> PNAS; 99(10): 6567-72. </p>'),
                      
                      HTML('<p> [4] Dudoit, S., Fridlyand, J. and Speed, T.P. (2002). <a href="http://amstat.tandfonline.com/doi/abs/10.1198/016214502753479248" target="_blank">Comparison of Discrimination Methods for the Classification of Tumors Using Gene Expression Data.</a> Journal of the American Statistical Association; 97(457): 77-87.</p>')
                      
                      
              )),
              column(3,imageOutput("images13")))
            
            #fluidRow(
             # column(3,imageOutput("images6")),
              #       column(3,imageOutput("images4")),
              #column(3,imageOutput("images14")),
              #column(3,imageOutput("images10"))
            #+)
             
    ),
    tabItem(tabName = "clustering",box(
      width = 3, status = "primary",
                                                    h4(strong("1. Pre-processing")),
                                                    h5(strong("a) Filtering")),
                                                    fluidRow(column(1),
                                                             column(10,
                                                                    checkboxInput(inputId = "nearZeroF", label = "Near-zero variance filtering", value = TRUE),
                                                                    checkboxInput(inputId = "varF", label = "Variance filtering", value = FALSE),
                                                                    checkboxInput(inputId = "stand", label = "Standardization", value = FALSE))
                                                    ),
                                                    conditionalPanel(condition = "input.varF",
                                                                     textInput(inputId = "maxVar", label = "Number of genes with maximum variance", value = "2000")
                                                    ),
                                                    h5(strong("b) Normalization")),
                                                    fluidRow(column(1),
                                                             column(10,
                                                                    radioButtons(inputId = "normMeth", label = "", choices = list("None" = "none", "DESeq median ratio" = "deseq", "TMM" = "TMM"), 
                                                                                 selected = "none",inline = TRUE,))
                                                    ),
                                                    h4(strong("2. Clustering for Distance Matrix")),
                                                     fluidRow(column(1),
                                                      column(11,
                                                      radioButtons(inputId = "clustering", label = "", selected = "Raw", 
                                                      choices = c("Raw" = "Raw", "Vst" = "Vst", "Rlog" = "Rlog",
                                                                  "VoomQW" = "VoomQW", "VoomPW" = "VoomPW", "Poisson" = "Poisson"
                                                                  ,"Model-based" = "modelbased", "EdgeR" = "EdgeR"))
                                                    )),   
                                                    #HTML('<br>'),
                                                    conditionalPanel(condition="input.clustering=='Raw'",
                                                                     fluidRow(column(6,selectInput("distancemeasure", "Distance Matrix Measures:", choices = c("Euclidean" = "euclidean", "Squared Euclidean" = "squared euclidean", "Manhattan" = "manhattan",
                                                                                                                                                               "Canberra" = "canberra", "Chebychev" = "chebychev", "Bray Curtis" = "bray curtis"
                                                                                                                                                               ,"Cosine Correlation" = "cosine correlation", "Pearson Correlation" = "pearson correlation"),
                                                                                                   selected = "Euclidean",width = '100%')),
                                                                              column(6,selectInput("methods", "Clustering Methods:", choices = c("Complete" = "complete", "Ward.D2" = "ward.D2", "Average" = "average",
                                                                                                                                                 "Mcquitty" = "mcquitty", "Kmeans" = "kmeans", "Kmedoids" = "kmedoids"),
                                                                                                   selected = "Complete",width = '100%'))
                                                                              
                                                   )),
      conditionalPanel(condition="input.clustering=='Vst'",
                       fluidRow(column(6,selectInput("distancemeasure", "Distance Matrix Measures:", choices = c("Euclidean" = "euclidean", "Squared Euclidean" = "squared euclidean", "Manhattan" = "manhattan",
                                                                                                                 "Canberra" = "canberra", "Chebychev" = "chebychev", "Bray Curtis" = "bray curtis"
                                                                                                                 ,"Cosine Correlation" = "cosine correlation", "Pearson Correlation" = "pearson correlation"),
                                                     selected = "Euclidean",width = '100%')),
                                column(6,selectInput("methods", "Clustering Methods:", choices = c("Complete" = "complete", "Ward.D2" = "ward.D2", "Average" = "average",
                                                                                                   "Mcquitty" = "mcquitty", "Kmeans" = "kmeans", "Kmedoids" = "kmedoids"),
                                                     selected = "Complete",width = '100%'))
                                
                       )),
      conditionalPanel(condition="input.clustering=='Rlog'",
                       fluidRow(column(6,selectInput("distancemeasure", "Distance Matrix Measures:", choices = c("Euclidean" = "euclidean", "Squared Euclidean" = "squared euclidean", "Manhattan" = "manhattan",
                                                                                                                 "Canberra" = "canberra", "Chebychev" = "chebychev", "Bray Curtis" = "bray curtis"
                                                                                                                 ,"Cosine Correlation" = "cosine correlation", "Pearson Correlation" = "pearson correlation"),
                                                     selected = "Euclidean",width = '100%')),
                                column(6,selectInput("methods", "Clustering Methods:", choices = c("Complete" = "complete", "Ward.D2" = "ward.D2", "Average" = "average",
                                                                                                   "Mcquitty" = "mcquitty", "Kmeans" = "kmeans", "Kmedoids" = "kmedoids"),
                                                     selected = "Complete",width = '100%'))
                                
                       )),
      conditionalPanel(condition="input.clustering=='VoomQW'",
                       fluidRow(column(6,selectInput("distancemeasure", "Distance Matrix Measures:", choices = c("Euclidean" = "euclidean", "Squared Euclidean" = "squared euclidean", "Manhattan" = "manhattan",
                                                                                                                 "Canberra" = "canberra", "Chebychev" = "chebychev", "Bray Curtis" = "bray curtis"
                                                                                                                 ,"Cosine Correlation" = "cosine correlation", "Pearson Correlation" = "pearson correlation"),
                                                     selected = "Euclidean",width = '100%')),
                                column(6,selectInput("methods", "Clustering Methods:", choices = c("Complete" = "complete", "Ward.D2" = "ward.D2", "Average" = "average",
                                                                                                   "Mcquitty" = "mcquitty", "Kmeans" = "kmeans", "Kmedoids" = "kmedoids"),
                                                     selected = "Complete",width = '100%'))
                                
                       )),
      conditionalPanel(condition="input.clustering=='VoomPW'",
                       fluidRow(column(6,selectInput("distancemeasure", "Distance Matrix Measures:", choices = c("Euclidean" = "euclidean", "Squared Euclidean" = "squared euclidean", "Manhattan" = "manhattan",
                                                                                                                 "Canberra" = "canberra", "Chebychev" = "chebychev", "Bray Curtis" = "bray curtis"
                                                                                                                 ,"Cosine Correlation" = "cosine correlation", "Pearson Correlation" = "pearson correlation"),
                                                     selected = "Euclidean",width = '100%')),
                                column(6,selectInput("methods", "Clustering Methods:", choices = c("Complete" = "complete", "Ward.D2" = "ward.D2", "Average" = "average",
                                                                                                   "Mcquitty" = "mcquitty", "Kmeans" = "kmeans", "Kmedoids" = "kmedoids"),
                                                     selected = "Complete",width = '100%'))
                                
                       )),
      conditionalPanel(condition="input.clustering=='Poisson'",
                       checkboxInput(inputId = "powertrans", label = "Power Transformation", value = TRUE)
      ),
      conditionalPanel(condition="input.clustering=='EdgeR'",
                       fluidRow(column(1),
                                column(11,
                                       selectInput("methods", "Clustering Methods:", choices = c("Complete" = "complete", "Ward.D2" = "ward.D2", "Average" = "average",
                                                                                                 "Mcquitty" = "mcquitty", "Kmeans" = "kmeans", "Kmedoids" = "kmedoids"),
                                                   selected = "Complete",width = '50%')
                                ))),
                                                    #textInput(inputId = "numbercluster", label = "Number of clusters", value = "2"),
      h4(strong("3. Number of clusters:")),
                                                    sliderInput("numbercluster", "", 1, 10, 2),
      actionButton("reset", "Run")
                                                    
    ),
    box(width = 9, status = "primary",
        verbatimTextOutput("trainConsole"),h4("Distance Matrix:"),div(style = 'overflow-x: scroll', DT::dataTableOutput('DistanceMatrix')),
        verbatimTextOutput("datamcutree")
    )
    ),
    tabItem(tabName = "heatmaps", fluidRow(column(3, box(weight=4,h4("Heatmap Features:"),
                                                        textInput(inputId = "mainname", label = "Main Name:", value = ""),
                                                        textInput(inputId = "xaxisname", label = "x-axis name:", value = ""),
                                                        textInput(inputId = "yaxisname", label = "y-axis name:", value = ""),
                                                        selectInput("dendrogram", "Dendrogram:", choices = c("none" = "none", "row" = "row", "column" = "column",
                                                           "both" = "both"),selected = "both",width = '50%'),
                                                        checkboxInput(inputId = "colsidecolor", label = "Column side color", value = FALSE),
                                                        checkboxInput(inputId = "rowsidecolor", label = "Row side color", value = FALSE),
                                                        selectInput("seriate", "Seriate:", choices = c("none" = "none", "mean" = "mean", "OLO" = "OLO",
                                                                                                             "GW" = "GW"),selected = "none",width = '50%'),
                                                        actionButton("runheatmap", "Run")
      )),
      column(9,plotlyOutput("heatmap"))
  )),
  tabItem(tabName = "PCA", fluidRow(column(3, box(weight=4,h4("PCA Features:"))),
                                    column(9,plotOutput("plot1"))))
  )
  ))


server <- function(input, output,session) { 
  output$images12 <- renderImage({
    return(list(src = "C:/Users/Lenovo/Documents/TEZ_WEBTOOL/images/images12.png",contentType = "image/png", width=250, height=250))
  }, deleteFile = FALSE) #where the src is wherever you have the picture
  
  output$images10 <- renderImage({
    return(list(src = "C:/Users/Lenovo/Documents/TEZ_WEBTOOL/images/images10.png",contentType = "image/png",width=300, height=300))
  }, deleteFile = FALSE) #where the src is wherever you have the picture
  
  output$images13 <- renderImage({
    return(list(src = "C:/Users/Lenovo/Documents/TEZ_WEBTOOL/images/images13.png",contentType = "image/png", width=300, height=300))
  }, deleteFile = FALSE) #where the src is wherever you have the picture
  
  output$images4 <- renderImage({
    return(list(src = "C:/Users/Lenovo/Documents/TEZ_WEBTOOL/images/images4.jpg",contentType = "image/png"))
  }, deleteFile = FALSE) #where the src is wherever you have the picture
  
  output$images15 <- renderImage({
    return(list(src = "C:/Users/Lenovo/Documents/TEZ_WEBTOOL/images/images15.png",contentType = "image/png", width=300, height=300))
  }, deleteFile = FALSE) #where the src is wherever you have the picture
  
  output$images6 <- renderImage({
    return(list(src = "C:/Users/Lenovo/Documents/TEZ_WEBTOOL/images/images6.jpg",contentType = "image/png",width=300, height=300))
  }, deleteFile = FALSE) #where the src is wherever you have the picture
  
  output$images14 <- renderImage({
    return(list(src = "C:/Users/Lenovo/Documents/TEZ_WEBTOOL/images/images14.png",contentType = "image/png", width=300, height=300))
  }, deleteFile = FALSE) #where the src is wherever you have the picture
  
  output$images1 <- renderImage({
    return(list(src = "C:/Users/Lenovo/Documents/TEZ_WEBTOOL/images/images1.jpg",contentType = "image/png", width=300, height=300))
  }, deleteFile = FALSE) #where the src is wherever you have the picture
  
  output$images7 <- renderImage({
    return(list(src = "C:/Users/Lenovo/Documents/TEZ_WEBTOOL/images/images7.jpg",contentType = "image/png", width=300, height=300))
  }, deleteFile = FALSE) #where the src is wherever you have the picture
  
  
  dataM <- reactive({  ## Data input.
    if(input$dataInput==1){  ## Load example data.
      if(input$sampleData==1){
        #verim <- read.xlsx("nafld.xlsx", sheetIndex = 1)
        verim <- read.table("Cervical_cancer.txt",header =TRUE,sep="")
      }
      
      else if(input$sampleData==2){
        verim <- read.table("alzheimer.txt",header =TRUE,sep="")
      }
      
    } 
    
    else if(input$dataInput==2){  ## Upload data.
      
      inFile <- input$upload
      mySep <- switch(input$fileSepDF, '1'=",",'2'="\t",'3'=";", '4'="")
      
      if (is.null(input$upload))  {return(NULL)}
      
      if (file.info(inFile$datapath)$size <= 10485800){
        if (input$radio3 == 1){
          verim <- read.xlsx(inFile$datapath, sheetIndex = 1)
          read.
        }
        else if (input$radio3 == 2){
          verim <- read.table(inFile$datapath, sep=mySep, header=TRUE, fill=TRUE, dec = ifelse(input$decimal, ",", "."))
        rownames(verim) <- verim[,1]
          verim <- verim[,-1]
        }
        else if (input$radio3 == 3){
          verim <- as.data.set(spss.system.file(inFile$datapath))
        }
        else if (input$radio3 == 4){
          verim <- read.csv(inFile$datapath, sep=mySep, header=TRUE, fill=TRUE, dec = ifelse(input$decimal, ",", "."))
        }
      }
      else print("File is bigger than 10MB and will not be uploaded.") 
    } 
    
    verim <- data.frame(verim)
    return(verim)   
  })
  
  output$RawData <- renderDataTable({
    
    #if (input$tabs1 == "Data Upload"){ 
      dataM()
   # }
  #}, options = list(iDisplayLength = 10)
  })
  
  datam_filtered <- reactive({  
    datam_filtered <- dataM()
    nZ_datam = caret::nearZeroVar(t(dataM()),saveMetrics = FALSE)
    if(length(nZ_datam) != 0) {
      datam_filtered = dataM()[-nZ_datam,]
    }
    ss = apply(datam_filtered,1,sd)
    ort = apply(datam_filtered,1,mean)
    cv = ss/ort
    siralama = order(cv, decreasing = TRUE)
    if (input$maxVar != 0){
      son = nrow(datam_filtered)*1
    }
    else{
      son = input$maxVar
    }
    datam_filtered <- datam_filtered[siralama[1:son],]
    
    if (input$clustering == "Raw" || input$clustering == "Vst" || input$clustering == "Rlog" || input$clustering == "VoomQW" || input$clustering == "VoomPW" 
        || input$clustering == "Poisson" || input$clustering == "modelbased" || input$clustering == "EdgeR")
    {
      if (input$normMeth == "none")
        datam_filtered = datam_filtered
      else if (input$normMeth == "deseq"){
        datam_filtered = DGEList(counts = as.matrix(dataM()))
        datam_filtered = calcNormFactors(datam_filtered, method = "RLE")
        datam_filtered = estimateCommonDisp(datam_filtered)$pseudo.counts
        datam_filtered <- datam_filtered[siralama[1:son],]}
      else if (input$normMeth == "TMM"){
        datam_filtered = DGEList(counts = as.matrix(dataM()))
        datam_filtered = calcNormFactors(datam_filtered, method = "TMM")
        datam_filtered = estimateCommonDisp(datam_filtered)$pseudo.counts
        datam_filtered <- datam_filtered[siralama[1:son],]}
    }
    
    
    return(datam_filtered)
  })
  
  output$trainConsole <- renderPrint({
    if(input$reset){ 
    
    cat("Model Summary:", "\n")
    cat("-----------------------------------------","\n")
    cat(paste("Raw Data"),"\n")
    cat(paste("   Data includes the read counts of ", dim(dataM())[1], " genes belong to ", dim(dataM())[2], " observations.", sep=""),"\n\n")
    
    nZ_datam = caret::nearZeroVar(t(dataM()),saveMetrics = FALSE)
    if(length(nZ_datam) != 0) {
      datam_filtered = dataM()[-nZ_datam,]
    }
    
    if (input$nearZeroF){
      cat(paste("Near-zero filtering"), "\n")
      cat(paste("   ", length(nZ_datam), " out of ", dim(dataM())[1]," genes are filtered.", sep=""), "\n\n")
    }
    
    if (input$varF){
      cat(paste("Variance filtering"), "\n")
      cat(paste("   ", input$maxVar, " genes are selected based on their maximum variance.",sep=""), "\n\n")
    }
    
    if (input$stand){
      datam_stand <- NULL
      datam_stand <- matrix(nrow=nrow(datam_filtered),ncol=nrow(datam_filtered))
      for(i in 1:ncol(datam_filtered)){
        for (j in 1:nrow(datam_filtered)) {
          datam_stand[j,i] <- (datam_filtered[j,i] - mean(datam_filtered[,i]))/sd(datam_filtered[,i])
        }
      }
      datam_filtered <- datam_stand[j,i]
    }
    
    if (input$stand){
      cat(paste("Standardization is applied."), "\n")
    }
    cat("-----------------------------------------","\n\n")
    

  
  datam_filtered <- reactive({  
    datam_filtered <- dataM()
    nZ_datam = caret::nearZeroVar(t(dataM()),saveMetrics = FALSE)
    if(length(nZ_datam) != 0) {
      datam_filtered = dataM()[-nZ_datam,]
    }
    ss = apply(datam_filtered,1,sd)
    ort = apply(datam_filtered,1,mean)
    cv = ss/ort
    siralama = order(cv, decreasing = TRUE)
    if (input$maxVar != 0){
      son = nrow(datam_filtered)*1
    }
    else{
      son = input$maxVar
    }
    datam_filtered <- datam_filtered[siralama[1:son],]
    
    if (input$clustering == "Raw" || input$clustering == "Vst" || input$clustering == "Rlog" || input$clustering == "VoomQW" || input$clustering == "VoomPW" 
        || input$clustering == "Poisson" || input$clustering == "modelbased" || input$clustering == "EdgeR")
    {
    if (input$normMeth == "none")
      datam_filtered = datam_filtered
    else if (input$normMeth == "deseq"){
      datam_filtered = DGEList(counts = as.matrix(dataM()))
      datam_filtered = calcNormFactors(datam_filtered, method = "RLE")
      datam_filtered = estimateCommonDisp(datam_filtered)$pseudo.counts
      datam_filtered <- datam_filtered[siralama[1:son],]}
    else if (input$normMeth == "TMM"){
      datam_filtered = DGEList(counts = as.matrix(dataM()))
      datam_filtered = calcNormFactors(datam_filtered, method = "TMM")
      datam_filtered = estimateCommonDisp(datam_filtered)$pseudo.counts
      datam_filtered <- datam_filtered[siralama[1:son],]}
    }
    
    
    return(datam_filtered)
  })
  
  datam_dist <- reactive({  
    if (input$clustering != "modelbased"){
      if (input$clustering == "Vst"){
        datam_dist <- DESeq2::varianceStabilizingTransformation(round(as.matrix(datam_filtered()+1)), blind=FALSE,fitType="local") }
      if (input$clustering == "Rlog"){
        datam_dist <-  DESeq2::rlog(round(as.matrix(datam_filtered())+1), blind=FALSE,fitType="local")}
      
      if (input$clustering == "Raw" || input$clustering == "Vst" || input$clustering == "Rlog")
        datam_dist = raw.distance(datam_filtered(),input$distancemeasure,"column")
      else if (input$clustering == "VoomQW")
        datam_dist = weighted.distance (datam_filtered(),input$distancemeasure,"column")
      else if (input$clustering == "VoomPW")
        datam_dist = weighted.distance2 (datam_filtered(),input$distancemeasure,"column")
      else if (input$clustering == "Poisson")
      {
        if (input$powertrans)
          datam_dist = as.matrix(PoissonDistanceAhu(t(datam_filtered()),type="none", transform=TRUE)$dd)
        else
          datam_dist = as.matrix(PoissonDistanceAhu(t(datam_filtered()),type="none", transform=FALSE)$dd)
      }
      else if (input$clustering == "EdgeR")
      {
        dge_veri <- DGEList(as.matrix(datam_filtered()))
        dge_veri <- calcNormFactors(dge_veri, method = "none")
        mds_veri <- plotMDS(dge_veri, top = 100, plot = FALSE)
        datam_dist <- mds_veri$distance.matrix
      }
      return(datam_dist)}
  })
  
  output$DistanceMatrix <- DT::renderDataTable({
  
  #output$DistanceMatrix <- DT::renderDataTable({
    #if(nrow(datam_dist()) == 0)
      #return("No data to show")
    
    if (input$clustering != "modelbased"){
      datatable(datam_dist(), extensions = c('Buttons','KeyTable', 'Responsive'), options = list(
       dom = 'Bfrtip',
        buttons = c('copy', 'csv', 'excel', 'pdf', 'print'), keys = FALSE, pageLength = 10
      ))
      
      
    }else{matrix(1,nrow=1,ncol=1)}
  })
  
  
  
  output$datamcutree <- renderPrint({
  if (input$methods == "complete" || input$methods == "ward.D2" || input$methods == "average" || input$methods == "mcquitty"){
    if (input$clustering == "Poisson")
      {datammethod = hclust(datam_dist()$dd, method = input$methods, members = NULL)}
    else
      {datammethod = hclust(as.dist(datam_dist()), method = input$methods, members = NULL)}
   
  datamcutree <- cutree(datammethod, input$numbercluster)}
  else if (input$methods == "kmeans")
  {
    if (input$clustering == "Poisson")
      datammethod = kmeans(datam_dist()$dd, input$numbercluster)
    else
    datamcutree = kmeans(as.dist(datam_dist()), input$numbercluster)
  }
  else if (input$methods == "kmedoids")
  {
    if (input$clustering == "Poisson")
      datammethod = pam(datam_dist()$dd, input$numbercluster)
    else
      datamcutree = pam(as.dist(datam_dist()), input$numbercluster)
  }
    
  if (input$clustering == "modelbased")
  {
    modelbased_matrix <- datam_filtered()
    GeneID_modelbased=1:nrow(modelbased_matrix)
    mb_RNAsEQ = RNASeq.Data(t(modelbased_matrix),Normalize = NULL, GeneID_modelbased)
    c0_KmeansPlus=KmeansPlus.RNASeq(mb_RNAsEQ, nK=input$numbercluster)$centers
    datamcutree=Cluster.RNASeq(data=mb_RNAsEQ, model="nbinom",centers=c0_KmeansPlus, method="EM")$cluster
    
  }
  cat("Cluster Predictions:", "\n")
  cat("-----------------------------------------","\n")
  cat(paste(datamcutree, sep=""), "\n\n")
  })
    }   

  })
  
  output$heatmap <- renderPlotly({
    if (input$runheatmap){
    HM_data <- t(log2(t(datam_filtered() + 0.5)/(apply(datam_filtered(), 2, sum) + 1) * 1e+06))
    rc <- colorspace::rainbow_hcl(nrow(HM_data))
    cc <- colorspace::rainbow_hcl(ncol(HM_data))
    if (input$rowsidecolor == TRUE && input$colsidecolor == FALSE)
    {rc <- colorspace::rainbow_hcl(nrow(HM_data))
    heatmaply(HM_data, xlab = as.character(input$xaxisname), ylab = as.character(input$yaxisname), 
              main = as.character(input$mainname),k_col = 3,grid_size = 1,colors = viridis,RowSideColors=rc,ColSideColors = NULL,
              dendrogram = input$dendrogram,seriate = input$seriate)  %>% layout(height=900,width=1000)}
    #if (input$colsidecolor == TRUE)
    #{cc <- colorspace::rainbow_hcl(ncol(HM_data))}

    else if (input$rowsidecolor == TRUE && input$colsidecolor == TRUE) {

        heatmaply(HM_data, xlab = as.character(input$xaxisname), ylab = as.character(input$yaxisname), 
              main = as.character(input$mainname),k_col = 3,grid_size = 1,colors = viridis,RowSideColors=rc,ColSideColors = cc,
              dendrogram = input$dendrogram,seriate = input$seriate)  %>% layout(height=900,width=1000)}
    
    else if(input$rowsidecolor == FALSE && input$colsidecolor == TRUE)
    {heatmaply(HM_data, xlab = as.character(input$xaxisname), ylab = as.character(input$yaxisname), 
               main = as.character(input$mainname),k_col = 3,grid_size = 1,colors = viridis,RowSideColors=NULL,ColSideColors = cc,
               dendrogram = input$dendrogram,seriate = input$seriate)  %>% layout(height=900,width=1000)}
    
    else if(input$rowsidecolor == FALSE && input$colsidecolor == FALSE)
    {heatmaply(HM_data, xlab = as.character(input$xaxisname), ylab = as.character(input$yaxisname), 
               main = as.character(input$mainname),k_col = 3,grid_size = 1,colors = viridis,RowSideColors=NULL,ColSideColors = NULL,
               dendrogram = input$dendrogram,seriate = input$seriate,plot_method = "plotly")  %>% layout(height=900,width=1000)}
    #heatmaply(mtcars, xlab = "Features", ylab = "Cars", 
      #        scale = "column",
       #       main = "Data transformation using 'scale'",yaxis_width=600)
    }
  })
  
  output$plot1 <- renderPlot({
    df <- dataM()
    df <- t(df)
    c_mydata <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
                  2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2)
    #autoplot(prcomp(t(df),data = t(df), colour = c_mydata,label = TRUE,frame = TRUE))
    iris.pca <- PCA(df#, quali.sup=5
                    )
    plot(iris.pca,axes = c(1, 2), #habillage = 2, 
         col.hab = c("black","blue", "red")
         #title = "Dataset projected onto PC1-2 Subspace"
         )
  })
  

}
shinyApp(ui, server)
