# librerias ------------------------------------------------------------------------------
library(readxl)
library(shiny)
library(shinyWidgets)
library(bslib)

library(plotly)
library(dplyr)
library(DT)

library(heatmaply)
library(adegenet)
library(pegas)
library(shinyalert)
library(shinydashboard)
library(shinyjs)

options(warn = -1)
shiny::addResourcePath("www", "www")

traducciones <- read.csv("www/traducciones.csv",stringsAsFactors = FALSE, fileEncoding = "latin1")
sample <- read.csv("www/sample.txt",sep = '\t')

# funciones adicionales ------------------------------------------------------------------
# opciones del pickerInput
options_picker <- pickerOptions(actionsBox = TRUE,hideDisabled = TRUE,liveSearch = TRUE)
# formato de tabla de salida
return_table <- function(datatble, col_2formt, idioma, formato_valor_combinado = FALSE) {
  lenguaje <- ifelse(idioma == "en", "en-GB.json", "es-MX.json")
  dt <- datatable(
    datatble,
    rownames = FALSE,
    extensions = 'Buttons',
    options = list(
      columnDefs = list(list(className = 'dt-center', targets = "_all")),
      dom = 'frtiB',
      language = list(url = lenguaje),
      paging = FALSE,
      fixedHeader = TRUE,
      scrollX = TRUE,
      scrollY = '50vh',
      scrollCollapse = TRUE,
      buttons = list(
        list(extend = "copy", text = 'Copiar'),
        list(extend = "csv", text = 'Descargar TXT', action = JS("exportToTxt"))
      )
    ),
    class = "tablashiny",
    selection = "none"
  ) %>%
    formatStyle(columns = col_2formt, 'text-align' = 'center') %>%
    formatStyle(columns = 1:ncol(datatble), fontSize = "12px")
  
  if (formato_valor_combinado && ("Valor combinado" %in% colnames(datatble) || "Combined value" %in% colnames(datatble))) {
    idx <- which(colnames(datatble) %in% c("Valor combinado", "Combined value"))
    dt <- DT::formatStyle(dt, columns = idx,
                          target = 'cell')
    
  }
  return(dt)
}

# funcion para llamar a las distintas funciones de estadígrafos
call_function <- function(name, values, size) {
  do.call(name, list(values, size))
}
# función para combinar columnas diploides
calculate_tabtotal <- function(tab,colindex){
  # tab = tabla totalleída, colindex = columnas de individuo y población
  cols <- names(tab[ , !(names(tab) %in% colindex)])
  cols1 <- cols[grep("\\.1$",cols,invert = TRUE)]
  cols2 <- cols[grep("\\.1$",cols)]
  tab1 <- tab[,c(colindex,cols1)]
  tab2 <- tab[,c(colindex,cols2)]
  names(tab2) <- sub("\\.1$", "", names(tab2))
  rbind(tab1, tab2)
}
# tabla de frecuencias alélicas
calculate_frAl <- function(colindex, tab, col_pop, pop_frec, cols_sel, size){
  len_ind <- length(colindex)
  tabfrec <- data.frame(a = 1:size)
  ## listar indices de columnas que se renderizaran
  list_2show <- c()
  for(i in cols_sel){
    rgx <- paste("^", i, "$",sep = "")
    list_2show <- append(list_2show, grep(rgx, colnames(tab)))
  }
  # homologar strs para calcular frecuencias, i.e. TH01 y TH01.1 es el mismo
  for (i in 1:length(cols_sel)){
    col1 <- tab[tab[[col_pop]] %in% pop_frec, list_2show[i]]
    col2 <- tab[tab[[col_pop]] %in% pop_frec, list_2show[i]+1]
    tabfrec[colnames(tab)[list_2show[i]]] <- append(col1, col2)
  }
  if(length(tabfrec)>2){
    tabfrec <- tabfrec[,2:length(tabfrec)]
  }else{
    onecol <- grep(colnames(tabfrec)[2],colnames(tabfrec))
    tabfrec <- tabfrec[onecol]
  }
  freqs <- data.frame(matrix(nrow = 0, ncol = 1))
  colnames(freqs) <- "Alelo"
  
  # calcular frecuencias
  for(col in colnames(tabfrec)){
    temp <- data.frame(unclass(table(tabfrec[,col])))
    colnames(temp) <- col
    temp <- cbind(Alelo = rownames(temp), temp)
    rownames(temp) <- 1:nrow(temp)
    colnames(freqs)[1] <- "Alelo"
    freqs <- merge(x = freqs, y = temp, by = "Alelo",all = TRUE)
    freqs$alelonum <- sapply(freqs$Alelo, as.numeric)
    freqs <- freqs[order(freqs$alelonum),]
    freqs <- subset(freqs, select = -alelonum)
    freqs[col] <- freqs[col]/size
  }
  freqs
}
# tabla de genotipos
calculate_geno <- function(tab, col_pop, pop_frec, cols_sel){
  tabpop <- tab[tab[[col_pop]] %in% pop_frec, ]
  tabgen <- data.frame(matrix(nrow = nrow(tabpop), ncol = length(cols_sel)))
  colnames(tabgen) <- cols_sel
  for (col in cols_sel) {
    col1_name <- paste0(col, ".1")
    tabgen[[col]] <- paste(tabpop[[col]], tabpop[[col1_name]], sep = "/")
  }
  tabgen
}
# tabla de frecuencias genotípicas
calculate_frGe <- function(tab, col_pop, pop_frec, cols_sel, size){
  tabgen <- calculate_geno(tab, col_pop, pop_frec, cols_sel)
  frGe <- data.frame(matrix(nrow = 0, ncol = 1))
  colnames(frGe) <- "Genotipo"
  freqhet <- data.frame(matrix(nrow = 1, ncol = length(cols_sel)))
  colnames(freqhet) <- cols_sel
  for(col in cols_sel){
    temp <- data.frame(unclass(table(tabgen[,col])))
    colnames(temp) <- col
    temp <- cbind(Genotipo = rownames(temp), temp)
    rownames(temp) <- 1:nrow(temp)
    colnames(frGe)[1] <- "Genotipo"
    frGe <- merge(x = frGe, y = temp, by = "Genotipo",all = TRUE)
    frGe[col] <- frGe[col]/(size/2)
  }
  frGe
}
# tabla de PEIF (SPFI) en inglés
calculate_stats <- function(cols_sel, freqs, frGe, tabgen, size,cols_str,hwe) {
  # aqui se calcula la frecuencia
  peif <- data.frame(matrix(ncol = length(cols_sel) + 1, nrow = 0))
  colnames(peif) <- append("PEIF", cols_sel)
  list_est <- c("HW","Ho","He", "PIC", "PD", "PE", "PM","TPI")
  
  for (i in 1:length(list_est)) {
    row <- c(list_est[i])
    if (row != "HW"){
      for (col in cols_sel){
        if(list_est[i] %in% c("PD","PM")){
          values <- as.numeric(unlist(frGe[col]))
        }else if(list_est[i] %in% c("Ho","PE","TPI")){
          values <- tabgen[col]
        }else if(list_est[i] %in% c("He","PIC","TPI")){
          values <- as.numeric(unlist(freqs[col]))
        }else{
          values <- cols_str
        }
        values <- values[!is.na(values)]
        row <- c(row, call_function(list_est[i], values, size))
      }
    }
    else{
      if (hwe[1] == "")
        values <- rep("-",length(cols_sel))
      else
        values <- hwe[cols_sel,4]
    }
    peif[i, ] <- c(row,values)
  }
  peif
}
# tabla de valores combinados
calculate_combinedvalues <- function(peif){
  # recibe la tabla de peif de todos los marcadores de la población
  a <- peif
  peif_num <- peif[, -1]
  peif_num <- as.data.frame(sapply(peif_num, as.numeric))
  if(dim(peif_num)[2]>1){
    valor_combinado0 <- apply(peif_num[c(1,4),], 1, function(x) "-") #(HW,PIC)
    valor_combinado1 <- apply(peif_num[2:3,], 1, function(x) (sum(x)/length(peif_num))) #(Ho, He)
    valor_combinado2 <- apply(peif_num[2:3,], 1, function(x) (prod(x)) )#(mean Ho, mean He)
    valor_combinado3 <- apply(peif_num[5:6,], 1, function(x) 1-(1-prod(x))) #(PD, PE)
    valor_combinado4 <- apply(peif_num[7:8,], 1, function(x) prod(x)) #(TPI,PM)
    names <-  c(head(a$PEIF, 3), "Π Ho", "Π He", tail(a$PEIF, -3))
    df_combinado <- data.frame(PEIF = names,comb = c(valor_combinado0[1],valor_combinado1,valor_combinado2,valor_combinado0[2],valor_combinado3,valor_combinado4))
    colnames(df_combinado)[colnames(df_combinado) == "comb"] <- "Valor combinado"
    df_combinado <- df_combinado[df_combinado$`Valor combinado` != "-", ]
    df_combinado$`Valor combinado` <- as.numeric(df_combinado$`Valor combinado`)
    df_combinado
  }
  else{
    peif
  }
}
# heatmap de población-alelo por str
heatmapstr <- function(tabtotal,colindex, str,idioma){
  T <- table(tabtotal[[colindex[2]]],tabtotal[,str])
  T_SUM <- rowSums(T)
  TP <- T/T_SUM
  coul <- colorRampPalette(c("#d4e6f1","#5499c7","#a9cce3", "#cd6155","#d98880"))(30)
  data <- as.matrix(TP)
  
  x <- ifelse(idioma == 'es', "Alelo", "Allele")
  y <- ifelse(idioma == 'es', "Población", "Population")
  
  heatmaply(main = str,xlab = x,ylab = y,
            na.value = "mintcream",as.data.frame.matrix(data),
            dendrogram = 'both',seriate = 'GW',,colors = coul
  ) %>% layout(height = 430)
}

# estadígrafos ---------------------------------------------------------------------------
HW <- function(tabgen, size){
  # recibe la tabla de genotipos
  tabgen <- tabgen %>% mutate_all(~gsub("\\.", "_", .))
  obj <- df2genind(tabgen,sep="/",ploidy=2, ncode=1,)
  set.seed(3744)
  return(hw.test(obj, B = 1000))
}
Ho <- function(col, size){
  # recibe la lista de genotipos de un str y la cantidad de alelos
  homo <- table(col)
  size <- size / 2
  # homocigotos (alelo1/alelo1, alelo2/alelo2, ..)
  homo <- rownames(homo)[grepl("^(\\S+)\\/\\1$", rownames(homo))]
  sum_homo <- 0
  for(h in homo)
    sum_homo <- sum_homo + sum(col==h)
  return(1-sum_homo/size)
}
He <- function(col, size){
  # recibe las frecuencias alélicas de un str y la cantidad de alelos
  # (∑(1 − ∑pij^2))/n
  return ((1-sum(col^2, na.rm=TRUE))*(size/(size-1)))
}
PIC <- function(col, size){
  # recibe las frecuencias alélicas de un str y la cantidad
  s <- 0
  for (i in seq_along(col)[-length(col)]){
    wrp_sum <- 0
    li <- i+1
    for(j in li:length(col)){
      wrp_sum <- wrp_sum + 2*col[i]^2*col[j]^2
    }
    s <- s + wrp_sum
  }
  # 1 −∑ pi^2 - (∑∑ 2pi^2pj^2)
  return (1-sum(col^2, na.rm=TRUE)-s)
}
PD <- function(col, size){
  # recibe la tabla de frecuencias genotípicas
  # 1 - PM
  return(1-sum(col^2))
}
PE <- function(col, size){
  # recibe la lista de genotipos de un str y la cantidad de alelos
  freqhet <- Ho(col, size)
  freqhomo <- 1 - freqhet
  # h^2(1-2hH^2)
  return(freqhet^2*(1-2*freqhet*freqhomo^2))
}
PM <- function(col, size){
  # recibe la tabla de frecuencias genotípicas
  return(sum(col^2))
}
TPI <- function(col, size){
  # recibe la lista de genotipos de un str y la cantidad de ellos
  # aquí size no se divide entre 2 por que ya lo divide Ho
  H <- 1-Ho(col,size)
  return(1/(2*H))
}
combine_pops <- function(populations){
  # recibe la lista de  poblaciones que hay en el dataset
  # calcula las posibles combinaciones entre poblaciones (se usa para las distancias)
  indexes <- list()
  for (i in 1:(length(populations)-1))
    for (j in (i+1):length(populations))
      indexes[[length(indexes)+1]] <- c(i,j)
  indexes
}
Da <- function(allfreqs,populations) {
  # recibe las frecuencias alélicas de cada marcador - población y la lista de poblaciones
  if (length(populations) > 1){
    distances <- data.frame(matrix(ncol = length(populations), nrow = length(populations)))
    rownames(distances) <- populations
    colnames(distances) <- populations
    markers <- names(data.frame(allfreqs[1]))[-1]
    # obtención de todas las combinaciones posibles de pares de poblaciones
    indexes <- combine_pops(populations)
    # iteración sobre los pares de poblaciones
    for (i in seq_along(indexes)) {
      sum_r <- 0
      for (marker in markers) {
        pop_name1 <- populations[indexes[[i]][1]]
        pop_name2 <- populations[indexes[[i]][2]]
        # extracción de frecuencias alélicas para las dos poblaciones
        df1 <- data.frame(allfreqs[indexes[[i]][1]])[,c("Alelo",marker)]
        df2 <- data.frame(allfreqs[indexes[[i]][2]])[,c("Alelo",marker)]
        merged_df <- merge(df1, df2, by = "Alelo", all = TRUE)
        names(merged_df) <- c(names(merged_df)[1],paste0(marker,"_",pop_name1),
                              paste0(marker,"_",pop_name2))
        # cálculo de la distancia de Nei para el marcador actual
        # ∑∑ sqrt(xij*yij)
        sum_mj <- sum(sqrt(merged_df[,ncol(merged_df)-1] * merged_df[,ncol(merged_df)]), na.rm=T)
        sum_r <- sum_r+sum_mj
      }
      # cálculo de la distancia genética para el par de poblaciones actual
      # 1-1/r*∑∑ sqrt(xij*yij)
      value <- 1-(1/length(markers)*sum_r)
      distances[indexes[[i]][2],indexes[[i]][1]] <- value
      distances[indexes[[i]][1],indexes[[i]][2]] <- value
    }
    for (i in seq_along(populations))
      distances[i,i]<-"0 "
    distances <- data.frame(names = row.names(distances), distances)
    colnames(distances) <- c(" ",colnames(distances)[-1])
    distances
  } else{
    NULL
  }
}
Fst <- function(allfreqs,populations){
  # recibe las frecuencias alélicas de cada marcador - población y la lista de poblaciones
  if (length(populations) > 1){
    distances <- data.frame(matrix(ncol = length(populations), nrow = length(populations)))
    rownames(distances) <- populations
    colnames(distances) <- populations
    # tuplas de combinaciones entre poblaciones
    indexes <- combine_pops(populations)
    for(i in seq_along(indexes)){
      # obtención de los dataframes de frecuencias alélicas para las dos poblaciones
      df1 <- allfreqs[[indexes[[i]][1]]]
      df2 <- allfreqs[[indexes[[i]][2]]]
      # homocigosidades esperadas de cada población
      Jx <- (sum(df1[-1]^2, na.rm = TRUE))/(length(colnames(df1[-1])))
      Jy <- (sum(df2[-1]^2, na.rm = TRUE))/(length(colnames(df1[-1])))
      # mezcla de poblaciones x*y
      xy  <- merge(df1, df2, by = "Alelo", all = TRUE)
      for (col in colnames(df1[-1])) {
        xy[[col]] <- xy[[paste0(col, ".x")]] * xy[[paste0(col, ".y")]]
      }
      xy[is.na(xy)] <- 0
      xy <- xy[!grepl("\\.x$|\\.y$", names(xy))]
      # homocigosidad esperada de todas as poblaciones
      Jxy <- (sum(xy[-1], na.rm = TRUE))/(length(colnames(xy[-1])))
      dum <- ((Jx +  Jy)/2)
      value <- (dum - Jxy)/(1 - Jxy)
      distances[indexes[[i]][2],indexes[[i]][1]] <- value
      distances[indexes[[i]][1],indexes[[i]][2]] <- value
    }
    for (i in seq_along(populations))
      distances[i,i]<-"0 "
    distances <- data.frame(names = row.names(distances), distances)
    colnames(distances) <- c(" ",colnames(distances)[-1])
    distances
  } else{
    NULL
  }
}

# interfaz -------------------------- -----------------------------------------------------
ui <- fluidPage(
  theme = bs_theme(version = 4, bootswatch = "cosmo"),
  navbarPage(a(href = "https://www.enacif.unam.mx", target = "_blank", 
               img(src = 'logo unam bw.png', width = "90px", height = "auto")
              ),id="inTabset", 
             theme = bs_theme(version = 4, bootswatch = "cosmo"),
             
             # cambios en estilos y código javascript ------------------------------------------------
             includeCSS("www/styles_shiny.css"),
             tags$head(tags$script(src = "www/script_shiny.js"),
                       tags$script(src = "message-handler.js")
             ),
             tags$div(class = "idioma-select",
                      tags$label(textOutput("tag_language"), `for` = "cambiar_idioma", style = "color: white; margin-bottom: 3px; display: block;"),
                      selectInput("cambiar_idioma", label = NULL,
                                  choices = c("Español" = "es", "English" = "en"),
                                  selected = "es")
             ),
             
             # paneles --------------------------------------------------------------------------
               tabPanel("Foren-STR",
                      sidebarLayout(
                        sidebarPanel(
                          # panel lateral: lector de archivos y seleccionador de columnas
                          tags$div(id = "file_container", fileInput("file", label = NULL)), 
                          h4(textOutput("file_txt")),
                          div(
                            class = "file-wrapper",
                            tags$label(
                              `for` = "file", class = "custom-btn", 
                              textOutput("updatefile_txt")
                            ),
                            div(
                              class = "file-name-container",
                              textOutput("file_name", inline = TRUE)
                            ),
                            actionButton("reloadbtn", "⟳️", class = "reload-btn", title = 'Recargar archivo')
                          ),
                          
                          
                          tags$head(tags$script(src = "message-handler.js")),
                          p(h5(textOutput("var_txt"))),
                          fluidRow(
                            column(6,
                                   pickerInput("ind", textOutput("ind_txt"), choices = c())),
                            column(6,
                                   pickerInput("pop", textOutput("pop_txt"), choices = c()))
                          ),
                          # panel condicional en función de la subpestaña seleccionada
                          # en ciertas pestañas, se pueden elegir varios strs
                          conditionalPanel(condition="input.tabselected==1 & 
                                    input.subtabselected1<4 |
                                    input.tabselected==2 & 
                                    input.subtabselected2==1 ",
                                           p(h5(textOutput("var_out_txt"))),
                                           pickerInput("cols_sel", textOutput("cols_sel_txt"), choices = c(),multiple=TRUE,
                                                       options = options_picker)
                          ),
                          # solo se puede seleccionar 1 str en la pestaña de heatmaps
                          conditionalPanel(condition="input.tabselected==1 & input.subtabselected1==4 ",
                                           p(textOutput("var_out_0_txt")),
                                           pickerInput("col_sel", textOutput("cols_sel_0_txt"), choices = c(),
                                                       options = options_picker)
                          ),
                          
                          # se pueden seleccionar más poblaciones en otras pestañas
                          conditionalPanel(condition="input.tabselected==1 & 
                                    input.subtabselected1>1 & input.subtabselected1<4 |
                                    input.tabselected==2 & 
                                    input.subtabselected2<=2 ",
                                           uiOutput("moreControls")
                          ),
                          width = 3
                        ),
                        mainPanel(
                          # pestañas y subpestañas
                          tabsetPanel(id = "tabselected",
                                      # pestaña 1
                                      tabPanel(title = textOutput("tab_description_txt"), value=1,
                                               mainPanel(
                                                 tabsetPanel(id = "subtabselected1",
                                                             tabPanel(value = 1, title = textOutput("tab_dataset_txt"), DT::dataTableOutput("tabdescrip")),
                                                             tabPanel(value = 2, title = textOutput("tab_allelic_txt"), uiOutput("allelicfreqs")),
                                                             tabPanel(value = 3, title = textOutput("tab_genotypic_txt"), uiOutput("genofreqs")),
                                                             tabPanel(value = 4,
                                                                      title = textOutput("tab_heatmap_txt"),
                                                                      h5(textOutput("tab_heatmap_0_txt")),
                                                                      plotlyOutput("heatmaps"),
                                                                      hr(),
                                                                      actionButton("prevBtn", icon = icon("arrow-left"), ""),
                                                                      actionButton("nextBtn", icon = icon("arrow-right"), "")
                                                             )
                                                 ), width = 11
                                               )
                                      ),
                                      # pestaña 2
                                      tabPanel(title = textOutput("tab_statistics_txt"), value=2,
                                               mainPanel(tabsetPanel(id = "subtabselected2",
                                                                     tabPanel(title=textOutput("stat_txt"), value=1,
                                                                              uiOutput("statistics")
                                                                     ),
                                                                     tabPanel(title=textOutput("com_val_txt"), value=2,
                                                                              uiOutput("statisticsall")
                                                                     ),
                                                                     tabPanel(title=textOutput("distances_txt"), value=3,
                                                                              h5(textOutput("da_distance_txt")),
                                                                              DT::dataTableOutput("tableda"),
                                                                              hr(),
                                                                              h5(textOutput("fst_distance_txt")),
                                                                              DT::dataTableOutput("tablefst"),
                                                                              hr()
                                                                     )),width = 12)
                                      )
                          ), width = 9
                        )
                      )
      
               ),
               tabPanel(title=textOutput("example_txt"),value = "panel1",
                        shinyjs::useShinyjs(),
                        fluidRow(
                          column(1),
                          column(10,
                                 h5(textOutput("structure_txt")),
                                 DT::dataTableOutput("samplefile")
                                 #downloadButton("downloadFile", textOutput("download_file_txt"))
                                 ),
                          column(1,align="center"))),
               tabPanel(title=textOutput("about_txt"),value = "panel1",
                        shinyjs::useShinyjs(),
                        fluidRow(
                          column(1),
                          column(10,
                            tags$ul(
                              h3(textOutput("about_0_txt")),
                              h5(textOutput("about_1_txt")),
                              tags$li(h4(textOutput("about_2_txt"))),
                              h5(textOutput("about_3_txt")),
                              tags$li(h4(textOutput("download_file_txt"))),
                              h5(textOutput("about_download_txt")),
                              downloadButton("downloadFileEsp", textOutput("download_file_esp_txt")),
                              downloadButton("downloadFileEng", textOutput("download_file_eng_txt")),
                              tags$li(h4(textOutput("about_4_txt"))),
                              h5(textOutput("about_5_txt")),
                              tags$li(h4(textOutput("about_6_txt"))),
                              h5("Castillo-Ortiz J., Salinas-Pineda L., Huerta-Pacheco N. S., Guardado-Estrada M."),
                              tags$li(h4(textOutput("about_7_txt"))),
                              tags$ol(
                                  tags$li(h5("Butler JM (2015) Chapter 1 - Data Interpretation Overview. In: Butler JM (ed) Advanced Topics in Forensic DNA Typing: Interpretation. Academic Press, San Diego, pp 3–24")),
                                  tags$li(h5("Kayser M, De Knijff P (2011) Improving human forensics through advances in genetics, genomics and molecular biology. Nat Rev Genet 12:179–192")),
                                  tags$li(h5("Ruitberg CM, Reeder DJ, Butler JM (2001) STRBase: a short tandem repeat DNA database for the human identity testing community")),
                                  tags$li(h5("Aguilar-Velázquez JA, Locia-Aguilar G, López-Saucedo B, et al (2018) Forensic parameters and admixture in seven geographical regions of the Guerrero state (South, Mexico) based on STRs of the Globalfiler ® kit. Ann Hum Biol 45:524–530. https://doi.org/10.1080/03014460.2019.1568548")),
                                  tags$li(h5("Sun G, McGarvey ST, Bayoumi R, et al (2003) Global genetic variation at nine short tandem repeat loci and implications on forensic genetics. European Journal of Human Genetics 11:39–49. https://doi.org/10.1038/sj.ejhg.5200902")),
                                  tags$li(h5("Butler JM (2015) STR Population Data Analysis. In: Advanced Topics in Forensic DNA Typing: Interpretation. Elsevier, pp 239–279")),
                                  tags$li(h5("Verma S, Pal R, Kandpal J, et al (2024) Developmental validation of NeoTyper autosomal STR kit. Human Gene 42:. https://doi.org/10.1016/j.humgen.2024.201348")),
                                  tags$li(h5("Excoffier L, Lischer HEL (2010) Arlequin suite ver 3.5: A new series of programs to perform population genetics analyses under Linux and Windows. Mol Ecol Resour 10:564–567. https://doi.org/10.1111/j.1755-0998.2010.02847.x")),
                                  tags$li(h5("Gouy A, Zieger M (2017) STRAF—A convenient online tool for STR data evaluation in forensic genetics. Forensic Sci Int Genet 30:148–151. https://doi.org/10.1016/j.fsigen.2017.07.007")),
                                  tags$li(h5("Gouy A, Zieger M (2025) STRAF 2: New features and improvements of the STR population data analysis software. Forensic Sci Int Genet 76:. https://doi.org/10.1016/j.fsigen.2024.103207")),
                                  tags$li(h5("Takezaki N, Nei M (1996) Genetic Distances and Reconstruction of Phylogenetic Trees From Microsatellite DNA")),
                                  tags$li(h5("R Core Team (2024) R: A Language and Environment for Statistical Computing")),
                                  tags$li(h5("Chang W, Cheng J, Allaire JJ, et al (2023) shiny: Web Application Framework for R")),
                                  tags$li(h5("Latter BDH (1972) SELECTION IN FINITE POPULATIONS WITH MULTIPLE ALLELES. 111. GENETIC DIVERGENCE WITH CENTRIPETAL SELECTION AND MUTATION"))
                              ),
                            )
                          ),
                          column(1,align="center")
                          )
                      ),
               nav_menu(
                 title = textOutput("links_txt"),
                 align = "right",
                 nav_item(tags$a(textOutput("links_1_txt"), href = "https://www.enacif.unam.mx")),
                 nav_item(tags$a("UNAM", href = "https://www.unam.mx"))
               )
       )
)