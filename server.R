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
sample <- read.csv("www/sample.txt",sep = "\t")

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

# servidor -------------------------------------------------------------------------------
server <- function(input, output, session) {
  # definición de variables --------------------------------------------------------------
  # variables reactivas (esperan un cambio en otra variable) 
  ind <- reactive({
    input$ind 
  }) # columna del individuo
  pop <- reactive({
    input$pop 
  }) # columna de la población
  col_index <- reactive({
    c(input$ind, input$pop) 
  }) # indice (individuo, población)
  cols_sel <- reactive({
    input$cols_sel 
  }) # strs seleccionados
  col_sel <- reactive({
    input$col_sel
  })  # str seleccionado (en caso de que solo se pueda elegir 1)
  num_pops <- reactive({
    input$num_pops
  })  # poblaciones a desplegar
  tab_added <- reactiveVal(FALSE) # evalúa si se agregó a tabla de combinaciones
  hwe_done <- reactiveVal(FALSE)
  stats_of_txt_ <- reactiveVal("")
  freqs_of_txt_ <- reactiveVal("")
  freqs_of_0_txt_ <- reactiveVal("")
  stats_data_txt_ <- reactiveVal("")
  
  # variables globales
  cols_str <- c() # columnas str disponibles
  cols_str1 <- c() # columnas str disponibles (1a copia)
  all_cols <- c() #todas las columnas que se renderizarán
  all_pops <- c() # poblaciones disponibles
  pop_frec <- list() # población seleccionada
  size <- c() # cantidad de individuos
  xlsx <- FALSE # determina si el archivo es .xlsx
  
  tab <- data.frame() # tabla leída
  freqs <- list() # frecuencias alélicas por marcador
  freqstot <- list() # frecuencias alélicas totales
  tabgen <- list() # listado de genotipos por marcador
  frGe <- list() # frecuencias genotipicas por marcador
  peif <- list() # tabla con los parámetros estadísticos por marcador
  peifall <- list() # tabla con los valores combinados
  hwe <- list() # tabla con los p-value de cada marcador
  
  PCO <- data.frame() # tabla de probabilidad de correcta ocurrencia
  distances <- list() # tablas con las distancias Da y FST
  

  # renderización inicial ----------------------------------------------------------------
  # lectura de archivo 
  observeEvent(input$reloadbtn, {
    session$reload()
  })
  
  observeEvent(input$file_clicked,{
    req(!is.null(input$file_clicked))
    if(input$file_clicked == 1)
      if (nrow(tab)>0){
        session$reload()
      }
  })
  
  observeEvent(input$file,{
    req(input$file)
    
    #regresar a la primera pagina pa q se renderice todo bien
    updateTabsetPanel(session, "tabselected", selected = "1")
    updateTabsetPanel(session, "subtabselected1", selected = "1")
    inFile <- input$file
    name <- inFile$datapath
    
    # cuando se lee un xlsx, las columnas diploides aparecen con un "...#num columna"
    xlsx <<- FALSE
    # leer archivo con base en la extensión del archivo 
    if (is.null(inFile))
      return(NULL)
    if(substr(name,nchar(name)-4,nchar(name))==".xlsx"){
      xlsx <<- TRUE
      tab<<-data.frame(read_excel(name))
    }
    else if(substr(name,nchar(name)-3,nchar(name))==".csv"){
      xlsx <<- FALSE
      tab<<-read.csv(name)
    }
    else if(substr(name,nchar(name)-3,nchar(name))==".txt"){
      xlsx <<- FALSE
      tab<<-read.table(name, header = TRUE)
    }
    
  })
  observeEvent(input$file, {
    req(nrow(tab)>0)
    # se "infiere" la columna ind y pob como la primera y la segunda
    updatePickerInput(session, "ind", choices = colnames(tab),
                      selected = colnames(tab)[1])
    updatePickerInput(session, "pop", choices = colnames(tab),
                      selected = colnames(tab)[2])
    if (idioma()=="es")
      session$sendCustomMessage(type = 'testmessage',
                                message = "Por favor verifique o valide sus datos cargados en la herramienta.")
    else
      session$sendCustomMessage(type = 'testmessage',
                                message = "Please check or validate your uploaded data in the tool.")
  })
  observeEvent(col_index(),{
    # definir columnas de indices (individuo y poblacion) 
    req(length(input$ind)>0 & length(pop())>0)
    new_cols <- c()
    
    # redefinir columnas para archivos .xlsx 
    if (xlsx){
      for (c in colnames(tab))
        if (!(c %in% col_index())){
          if (strtoi(sub("((.+)[.]{3})","",c)) %% 2 == 0)
            apnd = ".1"
          else
            apnd = ""
          new_cols <- append(new_cols, 
                             paste0(gsub("([.]{3})(\\d+)", "", c), apnd)) 
        }
      colnames(tab) <<- append(col_index(),new_cols)
    }
    
    # quitar anomalías en nombres de columnas
    newcols <- c()
    for(col in colnames(tab)) {
      # si la columna termina con ".1", se eliminan los puntos antes de ".1"
      if (grepl("\\.1$", col)) {
        newcols <- c(newcols, gsub("\\.(?!1$)", "", col, perl = TRUE))
      } else {
        newcols <- c(newcols, gsub("\\.", "", col))
      }
    }
    colnames(tab) <<- newcols
    
    # sustraer nombres de marcadores sin duplicados
    tab_names <- colnames(tab)
    tab_names <- gsub("([.]{1})(\\d{1})", "", tab_names)
    
    # seleccionar columnas en el input
    cols_str <<- setdiff(tab_names, col_index()) #columnas str disponibles
    cols_str1 <<- setdiff(colnames(tab), cols_str) 
    cols_str1 <<- setdiff(cols_str1, col_index()) #columnas str (copia dipploide)  disponibles
    all_cols <<- col_index() # columnas a renderizar (indice + str)
    tabtotal <<- calculate_tabtotal(tab,col_index()) # concatenado de copias  disponibles
    all_pops <<- unique(tab[,pop()]) # poblaciones disponibles
    updatePickerInput(session, "cols_sel", choices = cols_str[order(cols_str)],
                      selected =  cols_str) # STRs seleccionados 
  })
  
  # actualizar valores en listas desplegabes
  observeEvent(col_index(), {
    req(input$file)
    updatePickerInput(session, "ind", choices = col_index()[1],
                      selected=ind())
    updatePickerInput(session, "pop", choices = col_index()[2],
                      selected=pop())
    updatePickerInput(session, "col_sel", choices = cols_str[order(cols_str)],
                      selected = col_sel())
    updatePickerInput(session, "num_pops", choices = unique(tab[,pop()]),
                      selected = num_pops())
  })
  # evitar que se deseleccionen todas
  observeEvent(input$deselect, {
    if (input$deselect$value == 1) {
      input_id <- input$deselect$input_id
      if (input_id == "cols_sel") {
        updatePickerInput(session, "cols_sel", choices = cols_str[order(cols_str)],
                          selected = cols_str[order(cols_str)][1])
      } else if (input_id == "num_pops") {
        updatePickerInput(session, "num_pops", choices = unique(tab[,pop()]),
                          selected = unique(tab[,pop()])[1])
      }
      updateTextInput(session, "deselect", value = 0)
    }
  })
  disable("ind")
  disable("pop") 
  
  # outputs ------------------------------------------------------------------------------
  
  # archivo sample
  output$samplefile <- DT::renderDataTable({
    return_table(sample, colnames(sample), idioma())
  })
  # tabla principal
  output$tabdescrip <- DT::renderDataTable({
    req(length(input$cols_sel)>0)
    # mostrar columnas seleccionadas
    cols_dup <- c()
    for (col in colnames(tab))
      for (col_rend  in cols_sel())
        if(startsWith(col, col_rend))
          if (!(col %in% cols_dup)){
            cols_dup <- append(cols_dup,col)
          }
    all_cols <<- c(col_index(),cols_dup)
    # renderizar tablas con una sola columna
    updatePickerInput(session, "num_pops", choices = all_pops,
                      selected =  all_pops[1])
    if(length(all_cols)==1){
      onecol <- grep(all_cols, colnames(tab))
      return_table(tab[onecol], all_cols,idioma())
    }else{
      return_table(tab[,all_cols], all_cols,idioma())   
    }
  })
  # controles adicionales (poblaciones a seleccionar)
  output$moreControls <- renderUI({
    tagList(
      pickerInput("num_pops", textOutput("select_pop_txt"), choices = c(),multiple=TRUE,
                  options = pickerOptions(
                    actionsBox = TRUE,
                    hideDisabled = TRUE,
                    maxOptions = 100,
                    liveSearch = TRUE))
    )
  })

  # almaenamiento de tablas necesarias para hacer cálculos
  # por marcador
  observe({
    req(length(input$num_pops)>0)
    lapply(0:length(all_pops)+1, function(i){
      if(length(all_pops)>0){
        if(i == length(all_pops)+1)
          pop_frec[[i]] <<- all_pops
        else
          pop_frec[[i]] <<- all_pops[i]
        # tablas por str
        size[i] <<- 2 * length(tab[tab[[pop()]] %in% c(pop_frec[[i]]),pop()])
        freqs[[i]] <<- calculate_frAl(col_index(), tab, pop(), c(pop_frec[[i]]), cols_sel(), size[i])
        freqstot[[i]] <<- calculate_frAl(col_index(), tab, pop(), c(pop_frec[[i]]), cols_str, size[i])
        frGe[[i]] <<- calculate_frGe(tab, pop(), c(pop_frec[[i]]), cols_sel(), size[i])
        tabgen[[i]] <<- calculate_geno(tab, pop(), c(pop_frec[[i]]), cols_sel())
      }
    })
  })

  # renderizado de tablas (output es el renderizado y observe es el cálculo)
  # tablas de frecuencias
  # frecuencias alélicas 
  output$allelicfreqs <- renderUI({
    req(length(input$num_pops)>0)
    tagList(
      lapply(1:length(num_pops()), function(i) {
        uiOutput(paste0("containertable_pop", i))
      }),
      uiOutput("containertable_joint")  
    )
  })
  
  observe({
    req(length(input$num_pops)>0)
    freqs_of_txt_(traducir("freqs_of_txt"))
    txt <- freqs_of_txt_()
    lapply(1:length(num_pops()), function(i){
      output[[paste0("containertable_pop", i)]] <- renderUI({
        tagList(
          tags$div(
            style = "display: inline-block; float: left; ",
            h5(
              tags$span(txt), 
              tags$span(" "),  
              tags$span(num_pops()[i]) 
            ),
          ),
          DT::dataTableOutput(paste0("table_pop",i)),
          hr()
        )
      }) 
    })
    lapply(1:length(num_pops()), function(i){
      output[[paste0("table_pop", i)]] <- DT::renderDataTable({
        if (length(cols_sel())>0){
          freqsrend <- return_table(freqs[[i]], cols_sel(),idioma())
          freqsrend %>%
            formatRound(columns=(1:length(cols_sel())+1), digits = 4)
        }else{
          data.frame()
        }
      })
    })
    # tabla con la frecuencia conjunta
    output[["containertable_joint"]] <- renderUI({
      tagList(
        tags$div(
          style = "display: inline-block; float: left; ",
          h5(textOutput("freq_all_conj"))
        ),
        DT::dataTableOutput("table_joint"),
        hr()
      )
    })
    output[["table_joint"]] <- DT::renderDataTable({
      if (length(cols_sel())>0){
        joint_freq <- return_table(freqs[[length(all_pops)+1]], cols_sel(), idioma())
        joint_freq %>% formatRound(columns = (1:length(cols_sel())+1), digits = 4)
      }else{
        data.frame()
      }
    })
  })
  # frecuencias genotípicas
  output$genofreqs <- renderUI({
    req(length(input$num_pops)>0)
    if (length(num_pops())>0)
      tagList(
        lapply(1:length(num_pops()), function(i) {
          uiOutput(paste0("containertable_popg", i))
        }),
        uiOutput("containertable_jointg")  
      )
  })
  
  observe({
    req(length(input$num_pops)>0)
    freqs_of_0_txt_(traducir("freqs_of_0_txt"))
    txt <- freqs_of_0_txt_()
    if (length(num_pops())>0){
      lapply(1:length(num_pops()), function(i){
        output[[paste0("containertable_popg", i)]] <- renderUI({
          tagList(
            tags$div(
              style = "display: inline-block; float: left; ",
              h5(
                tags$span(txt), 
                tags$span(" "),  
                tags$span(num_pops()[i]) 
              ),
            ),
            DT::dataTableOutput(paste0("table_popg",i)),
            hr()
          )
        }) 
      })
      lapply(1:length(num_pops()), function(i){
        output[[paste0("table_popg", i)]] <- DT::renderDataTable({
          if (length(cols_sel())>0){
            freqsrend <- return_table(frGe[[i]], cols_sel(),idioma())
            freqsrend %>%
              formatRound(columns=(1:length(cols_sel())+1), digits = 4)
          }else{
            data.frame()
          }
        })
      })
      # tabla con la frecuencia genotípica conjunta
      output[["containertable_jointg"]] <- renderUI({
        tagList(
          tags$div(
            style = "display: inline-block; float: left; ",
            h5(textOutput("freq_gen_conj"))
          ),
          DT::dataTableOutput("table_jointg"),
          hr()
        )
      })
      output[["table_jointg"]] <- DT::renderDataTable({
        if (length(cols_sel())>0){
          joint_genofreq <- return_table(frGe[[length(all_pops)+1]], cols_sel(), idioma())
          joint_genofreq %>% formatRound(columns = (1:length(cols_sel())+1), digits = 4)
        }else{
          data.frame()
        }
      })
      
    }
  })
  
  # tablas de estadígrafos -----------------------------------------------------
  # por marcadores seleccionados
  output$statistics <- renderUI({
    req(length(input$num_pops)>0)
    if (length(num_pops())>0)
      tagList(
        lapply(1:length(num_pops()), function(i) {
          uiOutput(paste0("containertable_est", i))
        }),
        uiOutput("containertable_est_jointg")  
      )
  })
  
  observe({
    req(length(input$num_pops)>0)
    stats_of_txt_(traducir("stats_of_txt"))
    txt <- stats_of_txt_()
    if (!hwe_done())
      if(input$subtabselected2 == "1" && input$tabselected==2 ){
        if (idioma() == "es") {
          showModal(modalDialog(
            title = "Esperando",
            "Por favor espere mientras se ejecuta el proceso..."
          ))
        } else {
          showModal(modalDialog(
            title = "Waiting",
            "Please wait while the process is running..."
          ))
        }
        lapply(0:length(all_pops)+1, function(i){
          hwe[[i]] <<- HW(calculate_geno(tab, pop(), c(pop_frec[[i]]), cols_str))
          peif[[i]] <<- calculate_stats(cols_sel(), freqs[[i]], frGe[[i]],tabgen[[i]], size[i],cols_str,hwe[[i]])
          peifall[[i]] <<- calculate_combinedvalues(peif[[i]])
        })
        hwe_done(TRUE)
        removeModal()
      }

    lapply(1:(length(num_pops())), function(i) {
      p <- if (i == length(num_pops()) + 1) textOutput("all_pops_txt",inline=TRUE) else num_pops()[i]
      output[[paste0("containertable_est", i)]] <- renderUI({
        tagList(
          tags$div(
            style = "display: inline-block; float: left; ",
            h5(
              tags$span(txt), 
              tags$span(" "),  
              tags$span(p) 
            ),
          ),
          DT::dataTableOutput(paste0("table_est", i)),
          hr()
        )
      }) 
    })
    
    lapply(1:(length(num_pops())), function(i) {
      output[[paste0("table_est", i)]] <- DT::renderDataTable({
        if (length(cols_sel()) > 0) {
          peifrend <- return_table(peif[[i]], colnames(peif[[i]]), idioma())
          peifrend %>% 
            formatRound(columns = (1:length(cols_sel())+1), digits = 4)
        } else {
          data.frame()
        }
      })
    })
    # tabla con la los estadísticos para todas las poblaciones
    output[["containertable_est_jointg"]] <- renderUI({
      tagList(
        tags$div(
          style = "display: inline-block; float: left; ",
          h5(textOutput("stats_conj"))
        ),
        DT::dataTableOutput("table_stats_jointg"),
        hr()
      )
    })
    output[["table_stats_jointg"]] <- DT::renderDataTable({
      if (length(cols_sel())>0){
        joint_stats <- return_table(peif[[length(all_pops)+1]], cols_sel(), idioma())
        joint_stats %>% formatRound(columns = (1:length(cols_sel())+1), digits = 4)
      }else{
        data.frame()
      }
    })
  })
  
  # Valor combinado
  output$statisticsall <- renderUI({
    req(length(input$num_pops)>0)
    if (length(num_pops())>0)
      tagList(
        lapply(1:length(num_pops()), function(i) {
          uiOutput(paste0("containertable_estall", i))
        }),
        uiOutput("containertable_estall_jointg")  
      )
  })
  observe({
    req(length(input$num_pops)>0)
    stats_data_txt_(traducir("stats_data_txt"))
    txt <- stats_data_txt_()
    lapply(1:length(num_pops()), function(i){
      p <- num_pops()[i]
      output[[paste0("containertable_estall", i)]] <- renderUI({
        tagList(
          tags$div(
            style = "display: inline-block; float: left; ",
            h5(
              tags$span(txt), 
              tags$span(" "),  
              tags$span(p) 
            ),
          ),
          DT::dataTableOutput(paste0("table_estall",i)),
          hr()
        )
      })
    })
    lapply(1:length(num_pops()), function(i){
      output[[paste0("table_estall", i)]] <- DT::renderDataTable({
        if (length(cols_sel()) > 0) {
          peifrend_df <- peifall[[i]]
          table_data <- return_table(peifrend_df, colnames(peifrend_df), idioma(), formato_valor_combinado = TRUE)
          return(table_data)
        } else {
          return(data.frame())
        }
      })
    })
    # tabla con los valores combinados conjuntos
    output[["containertable_estall_jointg"]] <- renderUI({
      tagList(
        tags$div(
          style = "display: inline-block; float: left; ",
          h5(textOutput("comb_conj"))
        ),
        DT::dataTableOutput("table_estallcombined_jointg"),
        hr()
      )
    })
    output[["table_estallcombined_jointg"]] <- DT::renderDataTable({
      if (length(cols_sel()) > 0) {
        combpeifrend_df <- peifall[[length(all_pops)+1]]
        combtable_data <- return_table(combpeifrend_df, colnames(combpeifrend_df), idioma(), formato_valor_combinado = TRUE)
        return(combtable_data)
      } else {
        return(data.frame())
      }
    })
  })
  
  # distancias
  observe({
    req(length(input$num_pops)>0)
    req(length(freqstot)>=length(all_pops))
    if(input$subtabselected2=="3" & length(all_pops)<2){
      if (idioma()=="es")
        shinyalert("Distancias", 
                   "El archivo debe contener al menos 2 poblaciones para ver esta sección.")
      else
        shinyalert("Distances", 
                   "The file must contain at least 2 populations to view this section.")
    }
    if(length(all_pops)>1){
      distances[[1]] <<- Da(freqstot,all_pops) 
      distances[[2]] <<- Fst(freqstot,all_pops)
    }
  })
  # distancia de Nei
  output$tableda <- DT::renderDataTable({
    req(length(input$num_pops)>0)
    req(length(freqstot)>=length(all_pops))
    if(length(all_pops)>1){
      distancesrend <- return_table(distances[[1]],colnames(distances[[1]]),idioma())
      distancesrend %>%
        formatRound(columns = colnames(distances[[1]])[-1],digits = 4)
    }else{
      data.frame()
    }
  })
  # distancia Fst
  output$tablefst <- DT::renderDataTable({
    req(length(input$num_pops)>0)
    req(length(freqstot)>=length(all_pops))
    if(length(all_pops)>1){
      distancesrend <- return_table(distances[[2]],colnames(distances[[2]]),idioma())
      distancesrend %>%
        formatRound(columns = colnames(distances[[2]])[-1],digits = 4)
    }else{
      data.frame()
    }
  })
  
  # heatmaps
  output$heatmaps <- renderPlotly({
    req(length(input$num_pops)>0)
    if(input$subtabselected1=="4" & length(all_pops)<2){
      if (idioma()=="es")
        shinyalert("Mapa de calor", 
                   "El archivo debe contener al menos 2 poblaciones para ver esta sección.")
      else
        shinyalert("Heatmap", 
                   "The file must contain at least 2 populations to view this section.")

    }
    if(input$subtabselected1=="4" & length(all_pops)>=2){
      if (idioma()=="es")
        showModal(modalDialog(title = "Esperando",
                              "Por favor espere mientras se ejecuta el proceso..."))
      else
        showModal(modalDialog(title = "Waiting",
                              "Please wait while the process is running..."))
      Sys.sleep(2)
      removeModal()
    }
      if(length(all_pops)>1)
        heatmapstr(tabtotal, col_index(),col_sel(),idioma())
  })


  # botones (pasador de imágenes)
  observeEvent(input$prevBtn,{
    req(length(input$num_pops)>0)
    if(col_sel() != cols_str[1]){
      index <- match(col_sel(),cols_str)
      updatePickerInput(session, "col_sel", choices = cols_str,
                        selected =  cols_str[index-1])
    }
  })
  observeEvent(input$nextBtn,{
    req(length(input$num_pops)>0)
    if(col_sel() != cols_str[length(cols_str)]){
      index <- match(col_sel(),cols_str)
      updatePickerInput(session, "col_sel", choices = cols_str,
                        selected =  cols_str[index+1])
    }
  })
  
  # traducciones
  idioma <- reactiveVal("es")  
  
  traducir <- function(id) {
    traducciones[traducciones$id == id, idioma()]
  }
  observeEvent(input$cambiar_idioma,{
    idioma(input$cambiar_idioma)
    session$sendCustomMessage("idioma_cambiado", idioma())
  })
  observeEvent(input$num_pops,{
    idioma(input$cambiar_idioma)
    session$sendCustomMessage("idioma_cambiado", idioma())
  })
  observeEvent(input$col_sel,{
    idioma(input$cambiar_idioma)
    session$sendCustomMessage("idioma_cambiado", idioma())
  })
  observeEvent(input$HWE,{
    idioma(input$cambiar_idioma)
    session$sendCustomMessage("idioma_cambiado", idioma())
    session$sendCustomMessage("showRadioButtonHandler", list())
  })
  observeEvent(input$subtabselected1, {
    idioma(input$cambiar_idioma)
    session$sendCustomMessage("idioma_cambiado", idioma())
  })
  observeEvent(input$subtabselected2, {
    idioma(input$cambiar_idioma)
    session$sendCustomMessage("idioma_cambiado", idioma())
  })
  observeEvent(input$tabselected, {
    idioma(input$cambiar_idioma)
    session$sendCustomMessage("idioma_cambiado", idioma())
  })
  
  # Asignar traducciones dinámicas a los títulos
  output$esp_txt <- renderText(traducir("esp_txt"))
  output$eng_txt <- renderText(traducir("eng_txt"))
  output$tab_description_txt <- renderText(traducir("tab_description_txt"))
  output$tab_dataset_txt <- renderText(traducir("tab_dataset_txt"))
  output$tab_allelic_txt <- renderText(traducir("tab_allelic_txt"))
  output$tab_genotypic_txt <- renderText(traducir("tab_genotypic_txt"))
  output$tab_heatmap_txt <- renderText(traducir("tab_heatmap_txt"))
  output$tab_heatmap_0_txt <- renderText(traducir("tab_heatmap_0_txt"))
  output$file_txt <- renderText(traducir("file_txt"))
  output$var_txt <- renderText(traducir("var_txt"))
  output$ind_txt <- renderText(traducir("ind_txt"))
  output$pop_txt <- renderText(traducir("pop_txt"))
  output$cols_sel_txt <- renderText(traducir("cols_sel_0_txt"))
  output$cols_sel_0_txt <- renderText(traducir("cols_sel_txt"))
  output$var_out_txt <- renderText(traducir("var_out_txt"))
  output$calc_hw_txt <- renderText(traducir("calc_hw_txt"))
  output$stat_txt <- renderText(traducir("stat_txt"))
  output$com_val_txt <- renderText(traducir("com_val_txt"))
  output$distances_txt <- renderText(traducir("distances_txt"))
  output$da_distance_txt <- renderText(traducir("da_distance_txt"))
  output$fst_distance_txt <- renderText(traducir("fst_distance_txt"))
  output$tab_statistics_txt <- renderText(traducir("tab_statistics_txt"))
  output$example_txt <- renderText(traducir("example_txt"))
  output$structure_txt <- renderText(traducir("structure_txt"))
  output$about_txt <- renderText(traducir("about_txt"))
  output$about_0_txt <- renderText(traducir("about_0_txt"))
  output$about_1_txt <- renderText(traducir("about_1_txt"))
  output$about_2_txt <- renderText(traducir("about_2_txt"))
  output$about_3_txt <- renderText(traducir("about_3_txt"))
  output$about_4_txt <- renderText(traducir("about_4_txt"))
  output$about_5_txt <- renderText(traducir("about_5_txt"))
  output$about_6_txt <- renderText(traducir("about_6_txt"))
  output$about_7_txt <- renderText(traducir("about_7_txt"))
  output$links_txt <- renderText(traducir("links_txt"))
  output$links_1_txt <- renderText(traducir("links_1_txt"))
  output$download_file_txt <- renderText(traducir("download_file_txt"))
  output$about_download_txt <- renderText(traducir("about_download_txt"))
  output$download_file_esp_txt <- renderText(traducir("download_file_esp_txt"))
  output$download_file_eng_txt <- renderText(traducir("download_file_eng_txt"))
  output$select_pop_txt <- renderText(traducir("select_pop_txt"))
  output$freqs_of_txt <- renderText(traducir("freqs_of_txt"))
  output$freqs_of_0_txt <- renderText(traducir("freqs_of_0_txt"))
  output$all_pops_txt <- renderText(traducir("all_pops_txt"))
  output$all_pops_0_txt <- renderText(traducir("all_pops_0_txt"))
  output$stats_of_txt <- renderText(traducir("stats_of_txt"))
  output$stats_data_txt <- renderText(traducir("stats_data_txt"))
  output$var_out_0_txt <- renderText(traducir("var_out_0_txt"))
  output$freq_all_conj <- renderText(traducir("freq_all_conj"))
  output$freq_gen_conj <- renderText(traducir("freq_gen_conj"))
  output$stats_conj <- renderText(traducir("stats_conj"))
  output$comb_conj <- renderText(traducir("comb_conj"))
  output$tag_language <- renderText(traducir("tag_language"))
  
  output$updatefile_txt <- renderText(paste0(traducir("updatefile_txt"),"..."))
  output$file_name <- renderText({
    req(input$file)
    input$file$name
  })
  
  #boton recargar
  observe({
    if (idioma()=='es')
      runjs(sprintf('$("#reloadbtn").attr("title", "%s")', "Reiniciar memoria para cargar datos nuevos."))
    else
      runjs(sprintf('$("#reloadbtn").attr("title", "%s")', "Reset memory to load new data."))
  })
  
   # manuales
  manual_esp <- "www/Foren-STR Manual Español.pdf"
  output$downloadFileEsp <- downloadHandler(
    filename = function() {
      "Foren-STR Manual Español.pdf"
    },
    content = function(file) {
      file.copy(manual_esp, file)
    }
  )
  manual_eng <- "www/Foren-STR Manual English.pdf"
  output$downloadFileEng <- downloadHandler(
    filename = function() {
      "Foren-STR Manual English.pdf"
    },
    content = function(file) {
      file.copy(manual_eng, file)
    }
  )
}