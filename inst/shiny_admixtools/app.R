
library(admixtools)
library(igraph)
library(htmlwidgets)
library(shiny)
library(shinythemes)
library(plotly)
library(tidyverse)
library(magrittr)
library(RColorBrewer)
library(shinyjs)
library(shinyBS)
library(shinyFiles)
library(shinyalert)
library(shinyWidgets)
library(shinydashboard)


#options(shiny.sanitize.errors = FALSE)

edgemul = 1000
dtstyle = 'font-size:70%'
pl = 15
do = list(pageLength = pl, lengthMenu = list(c(pl, -1), c(pl, 'All')),
          buttons = c('copy', 'csv', 'excel'), scrollX = TRUE, dom = 'Blfrtip')
actionButton = function(...) actionBttn(..., style = 'unite', size = 'sm')
cols = gg_color_hue(5, 97, 15)
mutfuns = c(`SPR leaves`='spr_leaves', `SPR all`='spr_all', `Swap leaves`='swap_leaves', `Move admixture edge`='move_admixedge_once', `Flip admix`='flipadmix_random')

box = function (..., title = NULL, footer = NULL, status = NULL, solidHeader = FALSE,
                background = NULL, width = 6, height = NULL, collapsible = FALSE,
                collapsed = FALSE)
{
  boxClass <- "box"
  if (solidHeader || !is.null(background)) {
    boxClass <- paste(boxClass, "box-solid")
  }
  if (!is.null(status)) {
    validateStatus(status)
    boxClass <- paste0(boxClass, " box-", status)
  }
  if (collapsible && collapsed) {
    boxClass <- paste(boxClass, "collapsed-box")
  }
  # if (!is.null(background)) {
  #   validateColor(background)
  #   boxClass <- paste0(boxClass, " bg-", background)
  # }
  style <- NULL
  if (!is.null(height)) {
    style <- paste0("background:", background, "; height: ", validateCssUnit(height))
  }
  titleTag <- NULL
  if (!is.null(title)) {
    titleTag <- h3(class = "box-title", title)
  }
  collapseTag <- NULL
  if (collapsible) {
    buttonStatus <- status %OR% "default"
    collapseIcon <- if (collapsed)
      "plus"
    else "minus"
    collapseTag <- div(class = "box-tools pull-right", tags$button(class = paste0("btn btn-box-tool"),
                                                                   `data-widget` = "collapse", shiny::icon(collapseIcon)))
  }
  headerTag <- NULL
  if (!is.null(titleTag) || !is.null(collapseTag)) {
    headerTag <- div(class = "box-header", titleTag, collapseTag)
  }
  div(class = if (!is.null(width))
    paste0("col-sm-", width), div(class = boxClass, style = if (!is.null(style))
      style, headerTag, div(class = "box-body", style = paste('background:', background), ...), if (!is.null(footer))
        div(class = "box-footer", footer)))
}


tt = c('data' = 'Load data here',
       'dir' = 'Directory where precursors for f-statistics will be stored or have been stored.',
       'popfilediv' = 'Select a three column text file with individual labels in column 1 and population labels in column 3.',
       'graphfilediv' = 'A graph in AdmixTools format, or as a two column edge list.',
       'genofile1' = 'Select a Packedancestrymap .geno file to extract data from. This step should only be necessary once. After that, the extracted data can be used for any further analyses.',
       'adjpopdiv' = 'Assign individuals to populations. Population names can be changed, but should match population names in the graph file, if one is provided. Data will be read from disk and f2-statistics will be computed when switching to another sidebar item.',
       'bookmarkdiv' = 'Click this to save the current session. This will write files to disk and generate a URL which can be used to restore the state of the progam.',
       'qpgraph_fit' = 'Click this to see qpgraph fit statistics',
       'qpgraph_similar_minus1' = 'Fit all graphs that result from removing one admixture edge from the current graph',
       'qpgraph_similar_minusplus' = 'Fit all graphs that result from removing one admixture edge from the current graph and adding back one admixture edge',
       'qpgraph_similar_plus1' = 'Fit all graphs that result from adding one admixture edge to the current graph',
       'qpgraph_similar_decomposed' = 'Fit all trees that result from removing one edge from each admixture node',
       'qpgraph_similar_treeneighbors' = 'Fit all trees that result from removing one edge from each admixture node and then reattaching each edge to each other edge',
       'qpgraph_similar_flipadmix' = 'Fit all graphs that result from flipping the direction of each admixture edge (unless this would introduce a cycle)',
       'qpgraph_similar_update' = 'Set the current graph to the displayed similar graph',
       'qpgraph_similar_revert' = 'Display the original graph',
       'addleaf' = 'Choose a population which should be added to the current graph',
       'qpgraph_add' = 'Fit all graphs that result from adding this population to each edge',
       'optim_run' = 'Run optimization. The current graph is modified in a number of different ways, and each modified graph is evaluated. Can be run iteratively across several generations.',
       'optim_ngraphs' = 'Number of graphs to generate and evaluate per generation',
       'optim_ngen' = 'Number of generations. In each generation the best fitting graphs are selected for the next generation.',
       'mutfuns' = 'SPR leaves: Subree prune and regraft leaf edges\n
                    SPR all: Subree prune and regraft any edges\\n
                    Swap leaves: Swap two random leaves\\\n
                    Move admixture edge: Shift one admixture edge\\\\n
                    Flip admix: Flip the direction of one admixture edge',
       'econstraints_update_div' = 'Set limits to minimum and maximum edge lengths',
       'aconstraints_update_div' = 'Set limits to admixture weights',
       'qpgraphbin' = 'Choose the qpGraph executable on your computer',
       'qpadmbin' = 'Choose the qpAdm executable on your computer',
       'qpgraph_genofile' = 'Choose the .geno file for qpGraph to use',
       'qpadm_genofile' = 'Choose the .geno file for qpAdm to use',
       'usef3precomp' = 'Differences between the original and R version of qpGraph can be due to the estimation of f-statistics, or due to the fitting step. Click this to use the f-statistics generated in the original qpGraph to use in qpGraph in R.',
       'init_resample' = 'Re-evaluate graph with different random intializations of starting weights. This is always done by default, but only the parameters for the starting weights that minimize the score are shown. Here, the table will show all estimated graph parameters (score, edge lengths, admixture weights) for all sets of random starting parameters.',
       'ind_resample' = 'Re-evaluate the graph leaving out each sample at a time. Samples are only left out if they are not the only sample in their respective population. This can take long, as data has to be read and f2 computed separately for each sample that is being left out.',
       'snp_resample' = 'Re-evaluate the graph leaving out SNP blocks. With jackknife, each SNP block is left out at a time. With bootstrap, SNP blocks are resampled with replacement to match the original number of SNP blocks (usually around 700). The number of random resamplings can be specified.',
       'clear_edge' = 'Clear the current selection',
       'delete_edge' = 'Delete the selected edge or node. Has to be an admixture edge or leaf node.',
       'add_edge' = 'Add an admixture edge by connecting two selected non-admixture edges. Fails if it introduces cycles or attempts to connect admixture edges.',
       'randgraph' = 'Generate a random admixture graph. Uses the current populations, a specified number of admixture nodes, and a specified outgroup.',
       'lsqmode' = 'Least squares mode. By default, the computation of the graph likelihood score uses the inverse of the f3-statistic covariance matrix. When fitting many populations, this inverse can be unstable. Least squares mode sets the offdiagonal elements of the f3-statistic covariance matrix to zero, which makes the matrix inversion stable, but can introduce bias.',
       'f2_denom' = 'Scale the estimated f2-statistics by a constant factor. A value of around 0.278 usually converts f2 to FST scale.',
       'seed' = 'Set a random seed to always pick the same initial weights and to make the results reproducible.',
       'multiprocess' = 'Enable parallel evaluation of graphs in multiple R processes. Faster, but may not work on all operating systems. See "A Future for R: A Comprehensive Overview" for more details.',
       'collapse_edges' = 'Collapse all nodes which are separated by less then a certain threshold of estimated drift.',
       'f2corr' = 'Untick this to compute f-statistics without bias correction factor.',
       'addpop' = 'Add a population',
       'removepop' = 'Remove a population',
       'showseledge' = 'Selected edges or nodes are shown here',
       'downloadgraph' = 'Download the fit of the current graph as a .tsv file.',
       'numstart' = 'How often optimization should be initialized with random weights',
       'qpgraph_diag' = 'Factor to add to diagonal elements of f3-statistic covariance matrix before inversion. Increase this value if optimization doesn\'t converge',
       'qpadm_fit' = 'Return to the qpadm fit statistics',
       'qpadm_diag' = 'Factor to add to diagonal elements of f4-statistic covariance matrix before inversion. Increase this value if optimization doesn\'t converge',
       'qpadm_constrained' = 'Click this to constrain the admixture weights to be between 0 and 1',
       'data_load' = 'Only a placeholder at the moment. Data is loaded when switching from the "Data" tab to another tab',
       'b' = '',
       'c' = '',
       'f2' = 'These are the estimated f2 statistics',
       'f3' = 'These are the estimated and fitted f3 statistics',
       'opt' = 'This table contains information about the admixture weight optimization step for each iteration, such as the randomly drawn initial weights, as well as the final score ("value", the quantity that is being minimized) and convergence parameters from the R "optim" function. Run "Resample / Initial weights" to see how other estimated graph parameters depend on the initial weights.')


format_table = function(dat) {

  if(is.null(dat)) return(tibble())
  print('ft1')
  print(names(dat))
  popcols = which(str_sub(names(dat), 1, 3) == 'pop')
  names = unique(c(as.matrix(dat[,popcols])))
  print('ft2')
  nam = shortest_unique_prefixes(names, 3)
  pcols = which(names(dat) %in% c('p', 'p.value'))
  print('ft3')
  numcols = setdiff(which(sapply(dat, class) == 'numeric'), pcols)
  print('ft4')

  if(length(popcols) > 0) dat %<>% mutate_at(popcols, ~recode(., !!!set_names(nam, names)))
  if('graph' %in% names(dat)) dat %<>% select(-graph)
  if('edges' %in% names(dat)) dat %<>% select(-edges)
  print('ft5')
  out = dat %>%
    mutate_at(numcols, ~round(., max(2, 2 - ceiling(log10(max(abs(na.omit(.)))))))) %>%
    mutate_at(pcols, format.pval) %>%
    rename_at(vars(pcols), ~'p.value')
  out
}

aa = '.shiny-notification { height: 100px; width: 400px; position: fixed; top: 20%; left: calc(50% - 200px);;}'

ui = function(request) {
  tagList(useShinyjs(),
  useShinyalert(),
  chooseSliderSkin('Modern', color = 'black'),
  bsTooltip('tooltips', 'Click this to turn tooltips on or off. Needs to be reactivated after switching tabs.'),
  dashboardPage(
    dashboardHeader(title = 'AdmixTools', titleWidth = 300,
                    tags$li(class = 'dropdown', actionLink('tooltips', label = '', icon = icon('question'))),
                    tags$li(class = 'dropdown', id = 'bookmarkdiv', actionLink('._bookmark_', label = '', icon = icon('save')))),
  dashboardSidebar(tags$head(tags$style(HTML(aa))), width = 300,
    sidebarMenu(id = 'navbar',
    menuItem('Data', tabName = 'data', expandedName = 'data', id = 'datax', icon = icon('database'), startExpanded = TRUE,
             menuItem('Extract settings', tabName = 'extract_settings',
                      checkboxInput('fix_populations', 'Fix populations'),
                      numericInput('max_miss', 'max missing', value = 0.1, step = 0.01, min = 0, max = 1),
                      splitLayout(
                        numericInput('min_maf', 'min MAF', value = 0, step = 0.01, min = 0, max = 0.5),
                        numericInput('max_maf', 'max MAF', value = 0.5, step = 0.01, min = 0, max = 0.5)),
                      radioButtons('trans_extract', 'Mutations', choices = c('both', 'only transitions', 'only transversions')),
                      fileInput('keepsnps', NULL, buttonLabel = 'SNP list'),
                      numericInput('maxmem', 'max RAM in GB', value = 15, step = 1, min = 0)
                      )),
    menuItem('f2', tabName = 'f2', expandedName = 'f2', icon = icon('dice-two'),
             menuItem('Options', tabName = 'f2_options')),
    menuItem('f3', tabName = 'f3', expandedName = 'f3', icon = icon('dice-three'),
             menuItem('Options', tabName = 'f3_options')),
    menuItem('f4', tabName = 'f4', expandedName = 'f4', icon = icon('dice-four'),
             menuItem('Options', tabName = 'f4_options')),
    menuItem('qpAdm', tabName = 'qpadm', expandedName = 'qpadm', icon = icon('balance-scale'),
             actionLink('qpadm_fit', 'Fit'),
             menuItem('Compare', href = 'qpadm_comparison', tabName = 'qpadm_comparison',
                      splitLayout(
                        shinyFilesButton('qpadmbin', 'qpAdm bin', 'Select qpAdm bin', FALSE),
                        shinyFilesButton('qpadm_genofile', 'Geno file', 'Select Packedancestrymap geno file', FALSE)),
                      p(HTML('<center>or</center>')),
                      fileInput('qpadm_out', NULL, buttonLabel = 'qpAdm output'),
                      hr(),
                      actionButton('run_qpadm', 'Run')),
             menuItem('Resample', tabName = 'qpadm_resample'),
             menuItem('Options', tabName = 'qpadm_options',
                      numericInput('qpadm_diag', 'diag', value = 0.001, step = 0.001),
                      checkboxInput('qpadm_constrained', 'Constrain weights'))),
    menuItem('qpGraph', tabName = 'qpgraph', expandedName = 'qpgraph', id = 'qpgraph', icon = icon('project-diagram'),
             actionLink('qpgraph_fit', 'Fit'),
             menuItem('Load graph', fileInput('graphfile', NULL, placeholder = '', buttonLabel = 'Graph file')),
             menuItem('Modify graph', expandedName = 'qpgraph_modify', id = 'qpgraph_modify', tabName = 'qpgraph_modify',
                      hr(),
                      p('Selected'),
                      div(verbatimTextOutput('showseledge', placeholder = TRUE), tags$head(tags$style('#showseledge{max-height: 300px; background: white}'))),
                      splitLayout(
                        actionButton('clear_edge', 'Clear'),
                        actionButton('delete_edge', 'Delete'),
                        actionButton('add_edge', 'Add')),
                      hr(),
                      p('Randomize graph'),
                      uiOutput('nadmix'),
                      #splitLayout(checkboxInput('fixoutgroup', 'Fix outgroup', value = TRUE),
                      actionButton('randgraph', 'Randomize'),
                      hr()),
             menuItem('Similar graphs', tabName = 'qpgraph_similar', expandedName = 'qpgraph_similar',
                      actionLink('qpgraph_similar_minus1', 'Less admixture'),
                      actionLink('qpgraph_similar_minusplus', 'Same admixture'),
                      actionLink('qpgraph_similar_plus1', 'More admixture'),
                      actionLink('qpgraph_similar_decomposed', 'Component trees'),
                      actionLink('qpgraph_similar_treeneighbors', 'Tree neighbors'),
                      actionLink('qpgraph_similar_flipadmix', 'Flip admixture'),
                      menuItem('Add population', tabName = 'qpgraph_add',
                               uiOutput('qpgraph_add'),
                               actionBttn('qpgraph_add_run', 'Run')),
                      hr(),
                      actionLink('qpgraph_similar_update', 'Update graph'),
                      actionLink('qpgraph_similar_revert', 'Revert graph')),
             menuItem('Optimize', tabName = 'qpgraph_optim',
                      actionButton('optim_run', 'Run'),
                      sliderInput('optim_ngraphs', '# graphs', 1, 100, 10),
                      sliderInput('optim_ngen', '# generations', 1, 20, 1),
                      conditionalPanel('input.optim_ngen > 1', sliderInput('optim_nsel', '# selected', 1, 10, 3)),
                      sliderInput('optim_nrep', '# repeats', 1, 50, 1),
                      checkboxInput('optim_initrand', 'Start randomly'),
                      checkboxGroupInput('mutfuns', 'Graph modifications', choices = mutfuns, selected = mutfuns)),
             menuItem('Constrain', tabName = 'qpgraph_constraints',
                      menuItem('Drift', tabName = 'qpgraph_econstraints',
                        uiOutput('econstraints')),
                      menuItem('Admixture', tabName = 'qpgraph_aconstraints',
                        uiOutput('aconstraints'))),
             menuItem('Compare', tabName = 'qpgraph_comparison',
                      splitLayout(
                        shinyFilesButton('qpgraphbin', 'qpGraph bin', 'Select qpGraph bin', FALSE),
                        shinyFilesButton('qpgraph_genofile', 'Geno file', 'Select Packedancestrymap geno file', FALSE)),
                      checkboxInput('usef3precomp', 'Use qpGraph fstats in R'),
                      p(HTML('<center>or</center>')),
                      fileInput('qpgraph_out', NULL, buttonLabel = 'qpGraph output'),
                      hr(),
                      actionButton('run_qpgraph', 'Run')),
             menuItem('Resample', tabName = 'qpgraph_resample',
                      menuItem('Initial weights', tabName = 'qpgraph_resample_init',
                               sliderInput('numstart2', '# init', 1, 100, 10),
                               actionButton('init_resample', 'Run')),
                      menuItem('Individuals', tabName = 'qpgraph_resample_individuals',
                               actionButton('ind_resample', 'Run')),
                      menuItem('SNPs', tabName = 'qpgraph_resample_snps',
                               checkboxInput('boot', 'Bootstrap'),
                               conditionalPanel('input.boot',
                                 sliderInput('bootnum', '# boot', 0, 1000, 100)),
                               actionButton('snp_resample', 'Run'))),
             menuItem('Options', tabName = 'qpgraph_options',
                      menuItem('Fit', tabName = 'qpgraph_options_fit',
                        actionButton('options_update', 'Update'),
                        numericInput('qpgraph_diag', 'diag', value = 0.001, step = 0.001),
                        numericInput('numstart', '# init', value = 10, min = 1),
                        numericInput('f2_denom', 'f2 scale', value = 0.278, step = 0.001),
                        numericInput('seed', 'Random seed', value = NULL),
                        checkboxInput('lsqmode', 'lsqmode')),
                      menuItem('Plot', tabName = 'qpgraph_options_plot',
                        checkboxGroupInput('plotopt', '',
                          choices = c(`Reorder edges`='reorder_edges',
                                      `Shift edges down`='shift_down',
                                      `Collapse edges`='collapse_edges'),
                          selected = c('shift_down')),
                        sliderInput('collapse_threshold', 'Log10 collapse threshold', -6, 2, -3, step = 0.1))),
             menuItem('Download', tabName = 'qpgraph_download', expandedName = 'qpgraph_download',
                      div(radioButtons('downloadgraph_format', 'Graph format',
                                  choices = c(`ADMIXTOOLS`='admixtools', `Edge list`='edgelist'),
                                  selected = 'admixtools'),
                      downloadButton('downloadgraph', 'Save graph')))),
    menuItem('Options', tabName = 'options', expandedName = 'options', id = 'options', icon = icon('cogs'),
             checkboxInput('multiprocess', 'multiprocess', value = TRUE),
             checkboxInput('f2corr', 'f2 bias correction', value = TRUE)))),
  dashboardBody(uiOutput('dashboardbody'), tags$head(tags$style(HTML('.skin-black .sidebar .shiny-download-link a { color: #444; } .content {background-color: white } .main-header .logo {font-family: "GillSans", "Skia", "Avenir-Medium", "Herculanum", "Hiragino Maru Gothic Pro", "Maru Gothic Pro", "Hiragino", "Comic Sans MS"; font-weight: normal;};')))),
  skin = 'black'
  ))
  }


server = function(input, output, session) {

  showLog()
  dp = '~/mnt/o2/projects/admixprograms/indpairs_v42.1/'
  dp = '~/Downloads/testinds/'
  dp = '/Users/robert/Downloads/countdata_eas5/'
  dp = '~/Downloads/countdata_iosifgraph/'
  dp = ''
  #dp = '~/Documents/countdata_afafam/'
  #dp = '~/Documents/countdata_pavel2/'
  global = reactiveValues(countdir = dp,
                          allinds = list.dirs(paste0(dp, '/pairs'),F,F),
                          poplist = list())
  global$graph = NULL
  home = normalizePath('~')
  volumes = getVolumes()
  volumes = expr(c(home = normalizePath('~'), getVolumes()() %>% set_names))

  shinyDirChoose(input, 'dir', roots = volumes)
  shinyFileChoose(input, 'genofile1', roots = volumes)
  shinyFileChoose(input, 'qpgraph_genofile', roots = volumes)
  shinyFileChoose(input, 'qpadm_genofile', roots = volumes)
  shinyFileChoose(input, 'qpgraphbin', roots = volumes)
  shinyFileChoose(input, 'qpadmbin', roots = volumes)

  dir = reactive(input$dir)
  output$dirout = renderText({global$countdir})
  #shinyjs::hide('show_extract')
  #shinyjs::hide('show_popadjust')
  #shinyjs::hide('show_indselect')
  output$show_popadjust = renderText({as.numeric(show_popadjust())})
  output$show_extract = renderText({as.numeric(show_extract())})
  output$show_indselect = renderText({as.numeric(show_indselect())})
  output$genofile1out = renderText({parseFilePaths(eval(volumes), input$genofile1)$datapath %>% str_replace('\\.geno$', '')})

  observeEvent(input$dir, {
    if (!"path" %in% names(dir())) return()
    volumes = eval(volumes)
    print('volumes')
    print(volumes)
    print(dir()$path)
    print('input$dir')
    print(c(input$dir, class(input$dir)))
    print('xxx')
    print(parseDirPath(volumes, input$dir))
    print(str(parseDirPath(volumes, input$dir)))
    print(as.character(parseDirPath(volumes, input$dir)))
    global$countdir = parseDirPath(volumes, input$dir)
    global$iscountdata = 'indivs' %in% list.dirs(global$countdir,F,F)
    global$isf2data = !'indivs' %in% list.dirs(global$countdir,F,F) &
      'block_lengths.rds' %in% list.files(global$countdir)
    print('global$iscountdata')
    print(global$iscountdata)
    if(global$iscountdata) {
      global$allinds = list.dirs(paste0(global$countdir, '/pairs'),F,F)
    } else if(global$isf2data) {
      global$allinds = list.dirs(global$countdir,F,F)
      global$poplist = global$allinds %>% set_names(global$allinds)
      print('global$allinds')
      print(global$allinds)
    } else {

    }
    })

  observeEvent(input$navbar, {

    print('navbar detected')
    print(input$navbar)

  })

  # observeEvent(input$navbar == 'data', {
  #
  #   print('navbar == data!')
  #   global$databod = get_loaddata()
  #   global$bod = global$databod
  #
  # })

  observeEvent(input$sidebarItemExpanded, {
    print('expanded!')
    exp = input$sidebarItemExpanded
    print(exp)

    if(exp == 'options') {

      return()

    } else if(exp == 'data') {

      print('navbar == data!')
      global$databod = get_loaddata()
      global$bod = global$databod

    } else {

      nam = map(paste0('pop', seq_len(length(global$poplist))), ~input[[.]])
      if(is.null(nam[[1]])) {
        shinyalert('Error', global$poplist)
        return()
      }
      global$poplist = map(nam, ~input[[.]]) %>% set_names(nam)
      print('switch away from data')
      if(length(unlist(global$poplist)) == 0) {
        shinyalert('No samples selected!', '')
        return()
      }

      if(exp == 'f2') {

      print('switch to f2')

      global$bod = fluidRow(
        tabsetPanel(tabPanel('Pop pairs', plotlyOutput(outputId = 'f2heatmap', height = '800', width='auto')),
                    tabPanel('Ind pairs', plotlyOutput(outputId = 'f2heatmap_indiv', height = '800', width='auto')),
                    tabPanel('Pop table', div(DT::DTOutput('f2stats'), style = dtstyle))))

    } else if(exp == 'f4') {

      print('switch to f4')
      choices = names(global$poplist)
      global$f4_left = div(map(1:4, ~selectizeInput(paste0('f4pops', .),
                                 paste0('Population ', .),
                                 choices = choices,
                                 multiple = TRUE,
                                 selected = choices[seq(.*(.-1)/2+1, .*(.-1)/2+.)])),
                           checkboxInput('f4_allowdup', 'Allow duplicates'),
                           checkboxInput('f4_se', 'Show standard errors'))
      popsinuse = unique(c(input$f4pops1, input$f4pops2, input$f4pops3, input$f4pops4))

      # update
      map(1:4, ~{
        fp = paste0('f4pops', .)
        sel = input[[fp]]
        updateSelectizeInput(session, fp, selected = sel,
                             choices = union(sel, setdiff(names(global$poplist), popsinuse)))})

      global$bod = fluidRow(
        column(3, global$f4_left),
        column(9, tabsetPanel(tabPanel('Plot', plotlyOutput(outputId = 'f4plot', height = '800', width='auto')),
                              tabPanel('Table', div(DT::DTOutput('f4stats'), style = dtstyle)))))

    } else if(exp == 'qpadm') {

      print('switch to qpAdm')
      choices = names(global$poplist)
      if(is.null(input$qpadmpops1) || is.null(input$qpadmpops2) || is.null(input$qpadmpops3)) {
        global$qpadmpops1 = choices[1]
        global$qpadmpops2 = choices[2:3]
        global$qpadmpops3 = choices[4:6]
      } else {
        global$qpadmpops1 = input$qpadmpops1
        global$qpadmpops2 = input$qpadmpops2
        global$qpadmpops3 = input$qpadmpops3
      }

      #global$qpadm_leftpanel = qpadm_leftpanel()
      #global$qpadm_rightpanel = qpadm_rightpanel_fit()

      global$bod = fluidRow(
        column(3, qpadm_leftpanel()),
        column(9, qpadm_rightpanel_fit()))

    } else if(exp == 'qpgraph') {

      print('switch to qpGraph')

      if(is.null(global$graph)) {
        #shinyalert('Generating random graph', '')
        global$graph = random_admixturegraph(names(global$poplist), 2)
      } else if(!all(get_leafnames(global$graph) %in% names(global$poplist))) {
        shinyalert('Generating random graph because populations don\'t match',
                   setdiff(get_leafnames(global$graph), names(global$poplist)))
        global$graph = random_admixturegraph(names(global$poplist), numadmix(isolate(global$graph)))
      }

      global$qpg_right = qpg_right_fit()
      global$bod = global$qpgraphbod
    } else {
      shinyalert('Error', 'implement tab')
    }

    }

    if(input$tooltips %% 2 == 1) {
      print('add tt')
      imap(tt, ~addTooltip(session, .y, .x))
    }
  })


  observeEvent(global$graph, {
    print('observe global$graph')
    updateSelectInput(session, 'addleaf',
      choices = setdiff(names(global$poplist), get_leafnames(global$graph)))
  })

  get_graph = reactive({
    print('get graph')
    global$graph
  })


  observeEvent(global$qpadm_rightpanel, {
    print('global$qpadm_rightpanel change detected!')
    choices = names(global$poplist)
    global$qpadmbod = fluidRow(
      column(3, div(
          selectizeInput('qpadmpops1', 'Target', choices = choices, multiple = FALSE, selected = global$qpadmpops1),
          selectizeInput('qpadmpops2', 'Sources', choices = choices, multiple = TRUE, selected = global$qpadmpops2),
          selectizeInput('qpadmpops3', 'Outgroups', choices = choices, multiple = TRUE, selected = global$qpadmpops3),
               actionButton('qpadm_randomize', 'Randomize'))),
      column(9, global$qpadm_rightpanel))
  })
  observeEvent(global$qpg_right, {
    print('global$qpg_right change detected!')
    global$qpgraphbod = fluidRow(
      column(8, plotlyOutput(outputId = 'graphplot', height = '800', width='auto')),
      column(4, global$qpg_right))
  })

  observeEvent(global$databod, {global$bod = global$databod})
  observeEvent(global$qpadmbod, {global$bod = global$qpadmbod})
  observeEvent(global$qpgraphbod, {global$bod = global$qpgraphbod})

  observeEvent(input$navbar, {
    print('navbar selected')
    print(input$navbar)
  })

  observeEvent(input$addpop, {
    pops = names(global$poplist)
    i = 0
    repeat({i=i+1; newpop = paste0('pop', i); if(!newpop %in% pops) break})
    global$poplist[newpop] = NA
  })
  observeEvent(input$removepop, {
    global$poplist = head(global$poplist, -1)
  })

  observeEvent(input$popfile, {
    nam = c('ind', 'sex', 'pop')
    global$poplist = read_table2(input$popfile$datapath, col_names = nam, col_types = cols()) %>%
      select(-sex) %$% split(ind, factor(pop, levels = unique(pop)))

    if(length(global$allinds) == 0) { # f2data dir
      global$allinds = unlist(global$poplist)
    }

    print('observe')
    imap(names(global$poplist), ~{
      bttn = paste0('group', .y)
      print(bttn)
      observeEvent(input[[bttn]], {
        isolate({
          if(input[[bttn]] == 0) return()
          poplist = global$poplist
          #poplistbak = global$poplistbak
          countdir = global$countdir
          inds = poplist[[.x]]
          #if(is.null(global$poplistbak[[.x]])) {
          isgrouped = length(inds) == 1 && is_group(countdir, inds)
          if(!isgrouped) {
            if(.x %in% global$allinds) {
              shinyalert('Error', 'Group name cannot be identical to an existing sample or group name!')
              return()
            }
            print('add group')
            withProgress(message = paste0('grouping ', length(inds), ' samples...'), {
              group_samples(global$countdir, inds, .x, overwrite = TRUE)
            })
            global$allinds = union(global$allinds, .x)
            #global$poplistbak[[.x]] = inds
            global$poplist[[.x]] = .x
            updateSelectizeInput(session, .x, selected = .x, choices = global$allinds)
            #updateButton(session, bttn, 'ungroup')
          } else {
            print('delete group')
            withProgress(message = paste0('ungrouping ', length(inds), ' samples...'), {
              global$poplist[[.x]] = delete_groups(global$countdir, .x)
            })
            #global$poplist[[.x]] = global$poplistbak[[.x]]
            #global$poplistbak[[.x]] = NULL
            global$allinds = setdiff(global$allinds, .x)
            updateSelectizeInput(session, .x, selected = global$poplist[[.x]], choices = global$allinds)
            #updateButton(session, bttn, 'group')
          }
        })
      })
    })
  })


  observeEvent(input$graphfile, {
    gf = input$graphfile$datapath
    print('read graphfile')
    print(gf)
    colnums = read_lines(gf) %>% str_split('\\s+') %>% map_dbl(length)
    if(length(unique(colnums)) > 1) {
      print('admixtools format')
      global$graph = parse_qpgraph_graphfile(gf) %>% graph_from_edgelist()
    } else {
      print('two column format')
      global$graph = read_table2(gf, col_types = cols()) %>% select(1:2) %>%
        as.matrix %>% graph_from_edgelist()
    }
    print('read graphfile done')
  })

  observeEvent(input$qpgraph_out, {
    global$qpgraphout = parse_qpgraph_output(input$qpgraph_out$datapath)
    global$qpg_right = plotlyOutput('graphcomparison')
  })

  observeEvent(input$randgraph, {
    print('randg')
    #if(input$fixoutgroup) outpop = as_edgelist(global$graph)[1,2]
    global$graph = random_admixturegraph(get_leafnames(global$graph), input$nadmix, outpop = input$outpop)
  })

  observeEvent(input$qpgraph_similar_update, {
    print('qpgraph_similar_update')
    global$graph = global$selgraph
  })
  observeEvent(input$qpgraph_similar_revert, {
    print('qpgraph_similar_revert')
    global$selgraph = global$graph
    fit = get_fit()
    global$edges = fit$edges
    global$score = fit$score
    print(global$score)
  })

  observeEvent(input$navbar, {
    print('navbar:')
    print(input$navbar)
  })

  observeEvent(input$optim_run, {

    numgraphs = input$optim_ngraphs
    numgen = input$optim_ngen
    numsel = input$optim_nsel
    numrep = input$optim_nrep
    fudge = input$qpgraph_diag
    f2_denom = as.numeric(input$f2_denom)
    f2blocks = get_f2blocks()
    pops = dimnames(f2blocks)[[1]]
    f2blocks = f2blocks[pops, pops,]
    outpop = get_outpop(global$graph)
    nadmix = numadmix(global$graph)
    g = NULL
    if(!input$optim_initrand) g = global$graph

    print('running opt...')
    print(pops)
    print(dim(f2blocks))
    print(paste(numgraphs, numgen, numsel, nadmix))
    mutfuns = input$mutfuns
    print(mutfuns)

    withProgress(message = paste('Evaluating ', numsel+(numgraphs-numsel)*numgen*numrep,' graphs...'), {
    opt_results = find_graphs(f2blocks, outpop = outpop, numrep = numrep,
                              numgraphs = numgraphs, numgen = numgen, numsel = numsel,
                              numadmix = nadmix, initgraph = g, mutfuns = mutfuns, debug = FALSE,
                              fudge = fudge, f2_denom = f2_denom)
    })
    print(opt_results)
    opt_results %<>% filter(!is.na(score))
    print('opt done')
    if(nrow(opt_results) == 0) {
      shinyalert('No valid graphs in output!', 'Try increasing the "diag" parameter.')
    } else {
      winner = opt_results %>% top_n(1, -jitter(score)) %$% igraph[[1]]
      oldgraph = global$graph
      oldscore = global$score
      oldedges = global$edges
      withProgress(message = 'Fitting best graph...', {
        global$graph = winner
        fit = get_fit()
      })
      print('fit')
      print(str(fit))

      if(fit$score < oldscore) {
        shinyalert('Better graph found!', paste('Old score: ', round(oldscore, 2),
                                                '\nNew score: ', round(fit$score, 2)))
      } else {
        global$graph = oldgraph
        global$score = oldscore
        global$edges = oldedges
        fit = get_fit()
        shinyalert('No better graph found!')
      }
      print('graph')
      print(global$graph)
    }
    global$seed = NULL
  })


  qpgraphfun = reactive({

    input$options_update
    print('qpgraphfun')

    seed = isolate(as.numeric(input$seed))
    if(is.na(seed)) seed = global$seed
    if(input$usef3precomp && !is.null(global$precomp)) f3precomp = global$precomp
    else f3precomp = NULL

    isolate({
    numstart = as.numeric(input$numstart)
    f2_denom = as.numeric(input$f2_denom)
    if(is.null(f2_denom) || length(f2_denom) == 0) f2_denom = 1
    fudge = input$qpgraph_diag
    lsqmode = input$lsqmode
    print(paste('qpgraphfun:', numstart, seed, fudge, lsqmode, f2_denom))
    function(x, y, ...) {
      args = list(...)
      if(!'numstart' %in% names(args)) args[['numstart']] = numstart
      args = c(list(x, y), args, fudge = fudge, f2_denom = f2_denom, lsqmode = lsqmode,
               seed = seed, f3precomp = list(f3precomp))
      do.call(quietly(qpgraph), args)$result
    }
    })
  })

  get_fit = reactive({
    print('get_fit')
    f2blocks = get_f2blocks()
    g = global$graph
    input$qpgraph_similar_revert
    input$qpgraph_similar_update
    input$aconstraints_update
    input$econstraints_update
    global$qpgraph_ranges = NULL
    global$qpgraph_scorerange = NULL
    req(g, f2blocks)

    replace_null = function(x, y) ifelse(is.null(x), y, x)

    edges = g %>% igraph::as_edgelist() %>% as_tibble(.name_repair = ~c('from', 'to')) %>%
      group_by(to) %>% mutate(type = ifelse(n() > 1, 'admix', 'normal')) %>% ungroup %>%
      mutate(eid = paste(from, to, sep = ' -> '), eid2 = str_remove_all(eid, '\\.|>| '))
    isolate({
      if(isTRUE(global$useconstraints)) {
        edges %<>% mutate(lower = map_dbl(eid2, ~replace_null(input[[.]][1], 0)),
                          upper = map_dbl(eid2, ~replace_null(input[[.]][2], 1e9))) %>%
          mutate(upper = ifelse(type == 'normal', upper/edgemul, upper),
                 lower = ifelse(type == 'normal', lower/edgemul, lower)) %>% select(-type, -eid2)
      }
    })

    print('get_fit2')
    print(str(edges))
    withProgress(message = 'Fitting graph...', {
      fit = qpgraphfun()(f2blocks, edges)
    })
    print('get_fit3')
    global$edges = fit$edges
    global$score = fit$score
    fit
  })

  get_pdat = reactive({
    print('get_pdat')
    if(is.null(global$edges)) return()
    print('get_pdat2')
    init_plotdata()
  })

  get_f2blocks = reactive({
    print('get f2')
    poplist = global$poplist
    req(poplist, global$countdir)
    if(length(poplist) == 0) return()
    pops = rep(names(poplist), sapply(poplist, length))
    inds = NULL
    if(global$iscountdata) inds = unlist(poplist)
    withProgress(message = 'Reading data and computing f2...', {
      f2blocks = f2_from_precomp(global$countdir, inds = inds, pops = pops,
                                 apply_corr = input$f2corr)
    })
    print('get f2 done')
    f2blocks
  })

  get_f2dat_indiv = reactive({
    print('get f2 indiv')
    poplist = global$poplist
    req(poplist, global$countdir)
    if(length(poplist) == 0) return()
    withProgress(message = 'Reading data and computing f2...', {
      f2dat = f2_from_precomp(global$countdir, inds = unlist(poplist), pops = unlist(poplist), return_array = FALSE)
    })
    print('get f2 indiv done')
    f2dat
  })

  init_plotdata = reactive({

    print('init')
    input$qpgraph_similar_revert
    input$options_update
    graph = global$edges
    if('igraph' %in% class(graph)) {
      edges = graph %>% as_edgelist %>% as_tibble
    } else {
      edges = graph %>% as_tibble
      graph %<>% select(1:2) %>% as.matrix %>% graph_from_edgelist()
    }
    names(edges)[1:2] = c('from', 'to')
    if('collapse_edges' %in% input$plotopt) {
      withProgress(message = 'collapsing edges...', {
        edges %<>% collapse_edges(10^input$collapse_threshold)
      })
    }

    admixnodes = unique(edges[[2]][edges[[2]] %in% names(which(table(edges[[2]]) > 1))])
    graph = edges %>% select(1:2) %>% as.matrix %>% graph_from_edgelist()

    pos = data.frame(names(V(graph)), igraph::layout_as_tree(graph), stringsAsFactors = F) %>%
      set_colnames(c('node', 'x', 'y'))

    eg = edges %>% left_join(pos, by=c('from'='node')) %>%
      left_join(pos %>% transmute(to=node, xend=x, yend=y), by='to') %>%
      mutate(type = ifelse(to %in% admixnodes, 'admix', 'normal'))

    if('reorder_edges' %in% input$plotopt) {
      withProgress(message = 'reordering edges...', {
        eg = fix_layout(eg %>% rename(name = from), graph) %>% rename(from = name)
      })
    }
    if('shift_down' %in% input$plotopt) {
      withProgress(message = 'reordering edges...', {
        eg = admixtools:::fix_shiftdown(eg %>% rename(name = from), graph) %>% rename(from = name)
      })
    }

    if('weight' %in% names(edges)) {
      fa = function(x) paste0(round(x*100), '%')
      fe = function(x) round(x*edgemul)
      edges$label = ifelse(edges$type == 'admix', fa(edges$weight), fe(edges$weight))
      if(!is.null(global$qpgraph_ranges)) {
        edges %<>% left_join(global$qpgraph_ranges, by = c('from', 'to')) %>%
          mutate(label = ifelse(type == 'admix',
                                paste(fa(mid), paste0('[', fa(lo), '-', fa(hi), ']'), sep = '\n'),
                                paste(fe(mid), paste0('[', fe(lo), '-', fe(hi), ']'), sep = '\n')))
      }
    }
    if(!'label' %in% names(edges)) edges %<>% mutate(label='')
    e2 = edges %>% transmute(from, to, label) %>% left_join(count(., to), by = c('from'='to'))
    eg %<>% left_join(e2, by=c('from', 'to')) %>% mutate(indegree = replace_na(n, 0))
    nodes = eg %>% filter(to %in% get_leafnames(graph)) %>% rename(x=xend, y=yend, xend=x, yend=y) %>%
      transmute(name = to, x, y, xend, yend, to=NA, rownum = 1:n())

    namedList(nodes, eg)
  })

  show_popadjust = reactive(length(global$poplist) > 0)
  show_extract = reactive(isFALSE(global$iscountdata) && isFALSE(global$isf2data))
  show_indselect = reactive(show_extract() || isTRUE(global$iscountdata))

  observeEvent(input$minus1_cell_clicked, {
    row = input$minus1_rows_selected
    req(row)
    sel = get_minus1() %>% slice(row)
    global$selgraph = sel$graph[[1]]
    global$edges = sel$edges[[1]]
    global$score = sel$score[[1]]
    print('newgraph set...')
    print(numadmix(global$graph))
  })
  observeEvent(input$plus1_cell_clicked, {
    row = input$plus1_rows_selected
    req(row)
    sel = get_plus1() %>% slice(row)
    global$selgraph = sel$graph[[1]]
    global$edges = sel$edges[[1]]
    global$score = sel$score[[1]]
  })
  observeEvent(input$minusplus_cell_clicked, {
    row = input$minusplus_rows_selected
    req(row)
    sel = get_minusplus() %>% slice(row)
    global$selgraph = sel$graph[[1]]
    global$edges = sel$edges[[1]]
    global$score = sel$score[[1]]
  })
  observeEvent(input$decomposed_cell_clicked, {
    row = input$decomposed_rows_selected
    req(row)
    sel = get_decomposed() %>% slice(row)
    global$selgraph = sel$graph[[1]]
    global$edges = sel$edges[[1]]
    global$score = sel$score[[1]]
  })
  observeEvent(input$treeneighbors_cell_clicked, {
    row = input$treeneighbors_rows_selected
    req(row)
    sel = get_treeneighbors() %>% slice(row)
    global$selgraph = sel$graph[[1]]
    global$edges = sel$edges[[1]]
    global$score = sel$score[[1]]
  })
  observeEvent(input$flipadmix_cell_clicked, {
    row = input$flipadmix_rows_selected
    req(row)
    sel = get_flipadmix() %>% slice(row)
    global$selgraph = sel$graph[[1]]
    global$edges = sel$edges[[1]]
    global$score = sel$score[[1]]
  })
  observeEvent(input$addleaf_cell_clicked, {
    row = input$addleaf_rows_selected
    req(row)
    sel = get_addleaf() %>% slice(row)
    global$selgraph = sel$graph[[1]]
    global$edges = sel$edges[[1]]
    global$score = sel$score[[1]]
  })

  observeEvent(input$initresample_cell_clicked, {
    row = input$initresample_rows_selected
    req(row)
    sel = get_initresample0() %>% slice(row)
    global$edges = sel$edges[[1]]
    global$score = sel$score[[1]]
  })
  observeEvent(input$indresample_cell_clicked, {
    row = input$indresample_rows_selected
    req(row)
    sel = get_indresample0() %>% slice(row)
    global$edges = sel$edges[[1]]
    global$score = sel$score[[1]]
  })
  observeEvent(input$snpresample_cell_clicked, {
    row = input$snpresample_rows_selected
    req(row)
    sel = get_snpresample0() %>% slice(row)
    global$edges = sel$edges[[1]]
    global$score = sel$score[[1]]
  })


  observeEvent(input$multiprocess, {
    print('multiprocess')
    if(input$multiprocess) future::plan('multiprocess')
    else future::plan('sequential')
  })


  make_reactive_graphtabs = function(graphfun, inp) {
    eventReactive(inp, {
      print('graph neighbors')
      f2blocks = get_f2blocks()
      print('graph neighbors2')
      #g = global$graph
      g = get_graph()
      print(g)
      print(graphfun)
      withProgress(message = 'Finding graphs...', {
        out = g %>% graphfun
      })
      print(paste('found', nrow(out), 'graphs'))
      withProgress(message = paste('Evaluating', nrow(out), 'graphs...'), {
        qpgfun = qpgraphfun()
        out %<>% mutate(res = furrr::future_map(graph, ~qpgfun(f2blocks, .))) %>%
        unnest_wider(res) %>%
        select(starts_with('from'), starts_with('to'), starts_with('score'), starts_with('graph'), starts_with('edges'))
      })
      if(nrow(out) == 0) {
        shinyalert('No graphs could be evaluated!')
        return()
      }
      out %>% arrange(score)
    })
  }

  get_minus1 = make_reactive_graphtabs(graph_minusone, input$qpgraph_similar_minus1)
  get_plus1 = make_reactive_graphtabs(graph_plusone, input$qpgraph_similar_plus1)
  get_minusplus = make_reactive_graphtabs(graph_minusplus, input$qpgraph_similar_minusplus)
  get_decomposed = make_reactive_graphtabs(admixtools::decompose_graph, input$qpgraph_similar_decomposed)
  get_treeneighbors = make_reactive_graphtabs(decomposed_tree_neighbors, input$qpgraph_similar_treeneighbors)
  get_flipadmix = make_reactive_graphtabs(admixtools::graph_flipadmix, input$qpgraph_similar_flipadmix)

  get_addleaf = eventReactive(input$qpgraph_add_run, {
    f2blocks = get_f2blocks()
    g = global$graph
    withProgress(message = 'Finding graphs...', {
      out = g %>% graph_addleaf(input$addleaf)
    })
    withProgress(message = paste('Evaluating', nrow(out), 'graphs...'), {
      qpgfun = qpgraphfun()
      out %<>% mutate(res = furrr::future_map(graph, ~qpgfun(f2blocks, .))) %>%
        unnest_wider(res) %>%
        select(starts_with('from'), starts_with('to'), starts_with('score'), starts_with('graph'), starts_with('edges'))
    })
    if(nrow(out) > 0) out %<>% arrange(score)
    out
  })

  qpgraph_scorerange = function(out) {
    out %>%
      select(score) %>%
      filter(!is.na(score)) %>%
      arrange(score) %>%
      slice(1, round(n()/2), n()) %>%
      add_column(type = c('lo', 'mid', 'hi'), .before = 1) %>%
      deframe
  }
  qpgraph_ranges = function(out) {
    out %>%
      filter(!is.na(weight)) %>%
      group_by(from, to) %>%
      arrange(weight) %>%
      slice(1, round(n()/2), n()) %>%
      mutate(type = c('lo', 'mid', 'hi')) %>%
      ungroup %>%
      select(from, to, type, weight) %>%
      pivot_wider(names_from = type, values_from = weight)
  }

  get_initresample0 = reactive({
    print('get initresample')
    input$init_resample
    f2blocks = get_f2blocks()
    g = global$graph
    numstart = as.numeric(isolate(input$numstart2))

    withProgress(message = paste('Evaluating graph ', numstart, ' times...'), {
      qpgfun = qpgraphfun()
      out = tibble(id = seq_len(numstart)) %>%
        mutate(res = furrr::future_map(id, ~qpgfun(f2blocks, g, numstart = 1))) %>%
        unnest_wider(res) %>%
        select(id, starts_with('score'), starts_with('opt'), starts_with('graph'), starts_with('edges'))
    })
  })
  get_initresample = reactive({
    out = get_initresample0()
    if(nrow(out) > 0 && 'edges' %in% names(out)) {
      global$qpgraph_scorerange = qpgraph_scorerange(out)
      out %<>% unnest(edges)
      global$qpgraph_ranges = qpgraph_ranges(out)
      out %<>% mutate(edge = paste(from, to, sep = ' -> ')) %>%
        select(id, opt, score, edge, weight) %>%
        pivot_wider(names_from = edge, values_from = weight) %>%
        unnest(opt)
    }
    print('initresample done')
    out
  })

  get_indresample0 = reactive({
    print('get indresample')
    input$ind_resample
    poplist = global$poplist
    countdir = global$countdir
    inds = unlist(poplist)
    numsingle = sum(sapply(poplist, length) == 1)
    pops = rep(names(poplist), sapply(poplist, length))
    f2blocks = get_f2blocks()
    g = global$graph
    qpgfun = qpgraphfun()
    fun = make_resample_inds_fun(~qpgfun(..., y = g))
    withProgress(message = paste('Evaluating graph ', length(inds)-numsingle, ' times...'), {
      out = fun(countdir, inds = inds, pops = pops, numstart = 1)
    })
    out
  })
  get_indresample = reactive({
    out = get_indresample0()
    if(nrow(out) > 0) {
      print('head(out)')
      print(head(out))
      print(head(out$error))
      global$qpgraph_scorerange = qpgraph_scorerange(out)
      out %<>% unnest(edges)
      global$qpgraph_ranges = qpgraph_ranges(out)
      out %<>%
        mutate(edge = paste(from, to, sep = ' -> ')) %>%
        select(ind, pop, edge, score, weight) %>%
        pivot_wider(names_from = edge, values_from = weight)
    }
    print('indresample done')
    out
  })

  get_snpresample0 = reactive({
    f2blocks = get_f2blocks()
    input$snp_resample
    g = global$graph
    bootnum = ifelse(input$boot, input$bootnum, FALSE)
    withProgress(message = paste('Evaluating graph ', ifelse(bootnum, bootnum, dim(f2blocks)[[3]]),' times...'), {
      out = make_resample_snps_fun(~qpgraphfun()(..., y = g))(f2blocks, boot = bootnum, numstart = 1)
    })
  })
  get_snpresample = reactive({
    out = get_snpresample0()
    if(nrow(out) > 0 && 'edges' %in% names(out)) {
      print('head(out)')
      print(head(out))
      global$qpgraph_scorerange = qpgraph_scorerange(out)
      out %<>% unnest(edges)
      global$qpgraph_ranges = qpgraph_ranges(out)
      out %<>%
        mutate(edge = paste(from, to, sep = ' -> ')) %>%
        select(id, edge, score, weight) %>%
        pivot_wider(names_from = edge, values_from = weight)
    }
    print('snpres done')
    out
  })

  output$dashboardbody = renderUI({

    print('dashboardbody')
    global$bod

  })

  output$popselectors = renderUI({

    poplist = global$poplist
    indnames = unique(global$allinds)
    print('popselectors')
    print(indnames)

    sellist = imap(names(poplist), ~{
      #buttlab = ifelse(is.null(global$poplistbak[[.x]]), 'group', 'ungroup')
      isgrouped = length(poplist[[.x]]) == 1 && is_group(global$countdir, poplist[[.x]])
      buttlab = ifelse(isgrouped, 'ungroup', 'group')
      fillCol(tagList(textInput(paste0('pop', .y), NULL, .x, placeholder = 'pop name'),
                      shiny::actionButton(paste0('group', .y), buttlab),
                      selectizeInput(.x, NULL, indnames,
                            poplist[[.x]],
                            multiple = TRUE, width = 'auto'),
                   tags$style(type="text/css", paste0('.selectize-control{white-space:pre-wrap}'))), height = '1000px')})
    print('popselectors2')
    do.call(splitLayout, sellist)
  })

  output$econstraints = renderUI({
    global$useconstraints = FALSE
    print('econstraints')
    edgecols = isolate(get_edgecols())
    edges = isolate(get_fit()$edges) %>%
      mutate(nam = paste(from, to, sep = ' -> '), nam2 = str_remove_all(nam, '\\.|>| ')) %>%
      left_join(edgecols, by = c('nam2'='edge')) %>%
      filter(type == 'edge') %>% select(nam, nam2, weight, col)

    maxlen = ceiling(max(edges$weight))*edgemul

    print('econstraints 2')
    div(div(id = 'econstraints_update_div', actionButton('econstraints_update', 'Update')),
        sliderInput('any_econstraints', 'All edges', 0, maxlen, c(0, maxlen), dragRange = FALSE),
        edges %>%
        pmap(~{div(tags$style(HTML(paste0('[for=',..2,']+span>.irs-bar, [for=',..2,']+span>.irs>.irs-from, [for=',..2,']+span>.irs>.irs-to {background: ', ..4,'}'))),
                   sliderInput(..2, ..1, min = 0, max = maxlen, value = c(0, maxlen),
                           dragRange = FALSE, width = '100%'))}))
  })

  output$aconstraints = renderUI({
    global$useconstraints = FALSE
    print('aconstraints')
    input$add_edge
    input$delete_edge
    edgecols = isolate(get_edgecols())
    edges = isolate(get_fit()$edges) %>%
      mutate(nam = paste(from, to, sep = ' -> '), nam2 = str_remove_all(nam, '\\.|>| ')) %>%
      left_join(edgecols, by = c('nam2'='edge')) %>%
      filter(type == 'admix') %>% select(nam, nam2, weight, col)

    div(div(id = 'aconstraints_update_div', actionButton('aconstraints_update', 'Update')),
        edges %>%
        pmap(~{div(tags$style(HTML(paste0('[for=',..2,']+span>.irs-bar, [for=',..2,']+span>.irs>.irs-from, [for=',..2,']+span>.irs>.irs-to {background: ', ..4,'}'))),
                   sliderInput(..2, ..1, min = 0, max = 1, value = c(0, 1),
                               dragRange = FALSE, width = '100%'))})
    )
  })

  output$nadmix = renderUI({
    nad = 3
    if(!is.null(global$graph)) nad = numadmix(global$graph)
    div(sliderInput('nadmix', '# admix', 0, 10, nad),
    selectizeInput('outpop', 'Outgroup', names(global$poplist), selected = get_outpop(global$graph)))
  })

  output$qpgraph_add = renderUI({
    print('selectInput addleaf')
    choices = setdiff(names(global$poplist), get_leafnames(global$graph))
    print(choices)
    print(names(global$poplist))
    print(get_leafnames(global$graph))
    selectInput('addleaf', '', choices = choices)
  })

  output$showseledge = renderText({
    input$clear_edge
    print('showseledge1')
    se = global$seledge
    print('se')
    print(global$seledge)
    print(se)
    if(!is.null(se) && (is.na(se) || se == 'NA')) {
      print('is na')
      se = NULL
    }
    print('se2')
    print(se)
    se
  })

  observeEvent(input$any_econstraints, {

    print('update triggered by any_econstraints')
    req(global$graph)
    print('any_econstraints 2')
    edges = global$graph %>% igraph::as_edgelist() %>% as_tibble(.name_repair = ~c('from', 'to')) %>%
      group_by(to) %>% mutate(type = ifelse(n() > 1, 'admix', 'normal')) %>% ungroup %>%
      mutate(eid = paste(from, to, sep = '->') %>% str_remove_all('\\.|>| ')) %$% eid
    print('any_econstraints 3')
    if(is.null(isolate(input[[edges[1]]]))) return()
    edges %>% walk(~{updateSliderInput(session, ., value = c(input$any_econstraints[1],
                                                             input$any_econstraints[2]))})
    print('any_econstraints 4')
  })

  observeEvent(input$delete_edge, {

    print('del 1')
    eg = get_seledge()
    if(is.null(eg)) return()
    g = global$graph
    leaves = get_leafnames(g)
    print('del 2')
    print(eg)
    eg %<>% str_split('\n') %>% `[[`(1)
    if(str_detect(eg[1], ' -> ')) {
      eg %<>% str_split(' -> ')
      print(eg)
      for(i in 1:length(eg)) {
        gnew = delete_admix(g, eg[[i]][1], eg[[i]][2])
        g <<- g
        gnew <<- gnew
        if(!is.dag(gnew)) shinyalert('not dag')
        if(!is.connected(gnew)) shinyalert('not connected')
        g = gnew
      }
    } else {
      for(i in 1:length(eg)) {
        g = delete_leaf(g, eg[i])
      }
    }
    newleaves = setdiff(get_leafnames(g), leaves)
    if(length(newleaves) > 0) {
      newleaves <<- newleaves
      oldgraph <<- global$graph
      newgraph <<- g
      shinyalert('New leaves!', newleaves)
    }
    global$graph = g
    global$seledge = NULL
    print('delete_edge done')
    print(global$seledge)
  })

  observeEvent(input$clear_edge, {
    print('clear')
    global$seledge = NULL
  })

  observeEvent(input$add_edge, {
    print('add')
    eg = get_seledge() %>% str_split('\n') %>% `[[`(1) %>% str_split(' -> ')
    if(length(eg) != 2) {
      shinyalert('Error', 'Can only add edge which connects exactly two other edges!')
      return()
    }
    from = eg[[1]][2]
    to = eg[[2]][2]
    print(paste(from, to))
    g = global$graph
    leaves = get_leafnames(g)
    alert = function(x) shinyalert('Could not insert edge!', as.character(x))
    tryCatch({
      gnew = insert_admix_igraph(g, from, to, allow_below_admix = TRUE, desimplify = TRUE)
      }, warning = alert, error = alert)
    if(!exists('gnew') || is.null(gnew)) return()
    if(!igraph::is.dag(gnew)) {
      shinyalert('Could not insert edge!', 'Edge would create a cycle!')
      return()
    }
    global$graph = gnew
    global$seledge = NULL
    newleaves = setdiff(get_leafnames(gnew), leaves)
    if(length(newleaves) > 0) shinyalert('New leaves!', newleaves)
    print('add2')
  })

  observeEvent(input$run_qpadm, {
    print('observe run_qpadm')
    global$qpadmout = get_qpadmout()
    global$qpadm_rightpanel = plotOutput('qpadmcomparison')
  })
  observeEvent(input$run_qpgraph, {
      print('observe run_qpgraph')
      global$qpgraphout = get_qpgraphout()
      global$qpg_right = plotlyOutput('graphcomparison')
  })

  map(1:4,
  ~observeEvent(input[[paste0('f4pops', .)]], {
    if(input$f4_allowdup) return()
    popsinuse = unique(c(input$f4pops1, input$f4pops2, input$f4pops3, input$f4pops4))
    map(setdiff(1:4, .), ~{
      fp = paste0('f4pops', .)
      sel = input[[fp]]
    updateSelectizeInput(session, fp, selected = sel,
      choices = union(sel, setdiff(names(global$poplist), popsinuse)))})
  }))

  observeEvent(input$f4_allowdup, {
    map(1:4, ~{
      fp = paste0('f4pops', .)
      sel = input[[fp]]
      updateSelectizeInput(session, fp, selected = sel,
                           choices = names(global$poplist))})
  })

  map(1:3,
      ~observeEvent(input[[paste0('qpadmpops', .)]], {
        print('qpadmpops change')
        popsinuse = unique(c(input$qpadmpops1, input$qpadmpops2, input$qpadmpops3))
        map(setdiff(1:3, .), ~{
          fp = paste0('qpadmpops', .)
          sel = input[[fp]]
          global[[fp]] = sel
          updateSelectizeInput(session, fp, selected = sel,
                               choices = union(sel, setdiff(names(global$poplist), popsinuse)))})
      }))

  observeEvent(input$qpadm_randomize, {
    nam = sample(names(global$poplist))
    nl = floor(length(nam)/2)
    choices = list(nam[1], nam[2:nl], nam[(nl+1):length(nam)])
    map(1:3, ~{
      global[[paste0('qpadmpops', .)]] = choices[[.]]
      updateSelectizeInput(session, paste0('qpadmpops', .),
                                  selected = choices[[.]], choices = nam)})
  })

  get_qpadmout = reactive({

    print('get_qpadmout')
    get_qpadm()
    volumes %<>% eval
    pref = parseFilePaths(volumes, input$qpadm_genofile)$datapath
    pref %<>% str_replace('\\.geno$', '') %>% normalizePath(mustWork = FALSE)
    bin = parseFilePaths(volumes, input$qpadmbin)$datapath %>% normalizePath(mustWork = FALSE)
    if(length(bin) == 0) bin = '~/Downloads/admixtoools_hub_reldec19/bin/qpAdm'
    #if(length(pref) == 0) pref = '~/Downloads/eas2'
    #if(length(pref) == 0) pref = '~/Downloads/v42.1_small'
    #if(length(pref) == 0) pref = '~/Documents/v42.1_pavel'
    if(length(pref) == 0) pref = '~/Documents/v42.1_afafam'
    g = global$graph
    target = input$qpadmpops1
    left = input$qpadmpops2
    right = input$qpadmpops3

    withProgress(message = 'running original qpAdm...', {
      dir = tempdir()
      out = tryCatch(qpadm_wrapper(target, left, right, bin = bin, pref = pref, outdir = dir),
                     error = function(e) shinyalert('Error!', as.character(e)))
    })
    out
  })
  get_qpgraphout = reactive({

    print('get_qpgraphout')
    get_fit()
    volumes %<>% eval
    pref = parseFilePaths(eval(volumes), input$qpgraph_genofile)$datapath
    pref %<>% str_replace('\\.geno$', '') %>% normalizePath(mustWork = FALSE)
    bin = parseFilePaths(volumes, input$qpgraphbin)$datapath %>% normalizePath(mustWork = FALSE)
    print('bin')
    print(bin)
    if(length(bin) == 0) bin = '~/Downloads/admixtoools_hub_reldec19/bin/qpGraph'
    #if(length(pref) == 0) pref = '~/Downloads/eas2'
    #if(length(pref) == 0) pref = '~/Downloads/v42.1_small'
    #if(length(pref) == 0) pref = '~/Documents/v42.1_pavel'
    if(length(pref) == 0) pref = '~/Documents/v42.1_afafam'
    g = global$graph

    withProgress(message = 'running original qpGraph...', {
      dir = tempdir() %>% normalizePath()
      print('pref')
      print(pref)
      print('dir')
      print(dir)
      out = qpgraph_wrapper(g, bin, pref, outdir = dir, zthresh = 0)
      out$f2 = NULL
      print('parse fstats')
      global$precomp = parse_fstats(paste0(dir, '/fstats.out'))
    })
    out
  })

  observeEvent(input$econstraints_update, {
    global$useconstraints = TRUE
  })
  observeEvent(input$aconstraints_update, {
    global$useconstraints = TRUE
  })

  dto = function(x) div(DT::DTOutput(x), style = dtstyle)

  observeEvent(input$qpgraph_similar_minus1, { global$qpg_right = dto('minus1') })
  observeEvent(input$qpgraph_similar_plus1, { global$qpg_right = dto('plus1') })
  observeEvent(input$qpgraph_similar_minusplus, { global$qpg_right = dto('minusplus') })
  observeEvent(input$qpgraph_similar_decomposed, { global$qpg_right = dto('decomposed') })
  observeEvent(input$qpgraph_similar_treeneighbors, { global$qpg_right = dto('treeneighbors') })
  observeEvent(input$qpgraph_similar_flipadmix, { global$qpg_right = dto('flipadmix') })
  #observeEvent(input$addleaf, { global$qpg_right = dto('addleaf') })

  observeEvent(input$qpgraph_add_run, {global$qpg_right = dto('addleaf') })

  observeEvent(input$init_resample, {
    print('initresample!')
    global$qpg_right = div(DT::DTOutput('initresample'), style = dtstyle)
  })

  observeEvent(input$ind_resample, {
    print('indresample!')
    global$qpg_right = div(DT::DTOutput('indresample'), style = dtstyle)
  })

  observeEvent(input$snp_resample, {
    print('snp_resample!')
    global$qpg_right = div(DT::DTOutput('snpresample'), style = dtstyle)
  })

  observeEvent(input$qpadm_fit, {
    print('fit qpadm!')
    global$qpadm_rightpanel = qpadm_rightpanel_fit()
  })
  observeEvent(input$qpgraph_fit, {
    print('fit!')
    global$qpgraph_ranges = NULL
    global$qpgraph_scorerange = NULL
    global$qpg_right = qpg_right_fit()
  })

  qpadm_leftpanel = reactive({
    print('reactive qpadm_leftpanel')
    choices = names(global$poplist)
    div(
      selectizeInput('qpadmpops1', 'Target', choices = choices, multiple = FALSE, selected = global$qpadmpops1),
      selectizeInput('qpadmpops2', 'Sources', choices = choices, multiple = TRUE, selected = global$qpadmpops2),
      selectizeInput('qpadmpops3', 'Outgroups', choices = choices, multiple = TRUE, selected = global$qpadmpops3),
      actionButton('qpadm_randomize', 'Randomize'))
})

  qpadm_rightpanel_fit = reactive({
    tabsetPanel(id = 'qpadm_tabset',
                tabPanel('Weights', div(DT::DTOutput('qpadm_weights'), style = dtstyle)),
                tabPanel('f4', div(DT::DTOutput('qpadm_f4'), style = dtstyle)),
                tabPanel('Rank drop', div(DT::DTOutput('qpadm_rankdrop'), style = dtstyle)),
                tabPanel('Pop drop', div(DT::DTOutput('qpadm_popdrop'), style = dtstyle)))
  })

  qpg_right_fit = reactive({
    tabsetPanel(tabPanel('f2', id = 'f2', div(DT::DTOutput('f2'), style = dtstyle)),
                tabPanel('f3', div(DT::DTOutput('f3'), style = dtstyle)),
                tabPanel('opt', div(DT::DTOutput('opt'), style = dtstyle)),
                id = 'qpgraph_tabset')
  })


  get_loaddata = reactive({
    print('load_data')
    bh = '140px'
    div(
      fluidRow(
        box(width=4, height=bh, background = cols[1], h4('Select data directory'),
          div(
          shinyDirButton('dir', 'Browse', 'Upload'),
          verbatimTextOutput('dirout', placeholder = TRUE))),
        conditionalPanel('output.show_extract == "1"',
                         fluidRow(box(width=4, height=bh, background = cols[2], h4('Extract data'),
               splitLayout(
               div(shinyFilesButton('genofile1', 'Geno file', 'Select Packedancestrymap geno file', FALSE),
                   verbatimTextOutput('genofile1out', placeholder = TRUE)),
               actionButton('extract_counts', 'Extract counts'))))),
        conditionalPanel('output.show_indselect == "1"', (
          box(width=4, height=bh, background = cols[3], h4('Select .ind file'),
            div(fileInput('popfile', NULL, placeholder = '', buttonLabel = 'Ind file'), id = 'popfilediv')))),
        #column(3, box(width=12, height=bh, background = cols[4], h4('Select graph file'),
        #    div(fileInput('graphfile', NULL, placeholder = '', buttonLabel = 'Graph file'), id = 'graphfilediv'))),
      fluidRow(column(12, conditionalPanel('input.extract_counts > 0', verbatimTextOutput('console')))),
      #box(width=6, textOutput('console')),
      div(style = 'visibility: hidden', verbatimTextOutput('show_popadjust')),
      div(style = 'visibility: hidden', verbatimTextOutput('show_extract')),
      div(style = 'visibility: hidden', verbatimTextOutput('show_indselect')),
      conditionalPanel('output.show_popadjust == "1"', fluidRow(column(12,
        box(width=12, background = cols[5],
            fluidRow(column(2, h4('Adjust populations', id = 'adjpopdiv')),
                    column(1, div(actionButton('removepop', '–'),
                        actionButton('addpop', '+')))),
            fluidRow(column(11, uiOutput('popselectors')))))))))
  })

  observeEvent(input$extract_counts, {
    if(is.null(input$genofile1)) {
      shinyalert('Error!', 'Need to select .geno file!')
      return()
    }
    volumes %<>% eval
    pref = parseFilePaths(volumes, input$genofile1)$datapath %>% str_replace('\\.geno$', '')

    oldnam = names(global$poplist)
    nam = map(paste0('pop', seq_len(length(global$poplist))), ~input[[.]])
    if(is.null(nam[[1]])) return()
    global$poplist = map(oldnam, ~input[[.]]) %>% set_names(nam)

    poplist = global$poplist
    inds = unlist(poplist)
    pops = rep(names(poplist), sapply(poplist, length))
    transitions = input$trans_extract %in% c('both', 'only transitions')
    transversions = input$trans_extract %in% c('both', 'only transversions')

    if(input$fix_populations) extract_data = function(...) extract_f2(pops = pops, ...)
    else extract_data = extract_counts

    print('inds')
    print(inds)

    print('pref')
    print(pref)

    print('global$countdir')
    print(global$countdir)

    withCallingHandlers({
      shinyjs::html('console', '')
      print(list(pref, global$countdir, inds = inds,
                 maxmiss = input$max_miss, minmaf = input$min_maf, maxmaf = input$max_maf,
                 transitions = transitions, transversions = transversions,
                 keepsnps = input$keepsnps, maxmem = input$maxmem*1e3))
      extract_data(pref, global$countdir, inds = inds,
                   maxmiss = input$max_miss, minmaf = input$min_maf, maxmaf = input$max_maf,
                   transitions = transitions, transversions = transversions,
                   keepsnps = input$keepsnps, maxmem = input$maxmem*1e3)
    },
    message = function(m) {
      shinyjs::html(id = 'console', html = m$message, add = TRUE)
    })

  })


  observeEvent(event_data('plotly_click', source = 'src'), {

    print('observe seledge')
    ed = event_data('plotly_click', source = 'src')
    isolate({
      plt = plotly_graph()
      newe = ed_to_name(ed, plt)
      oldedge = global$seledge
      if(!is.null(oldedge) && is.na(oldedge)) {
        print('seledge is na. why?')
        oldedge = NULL
      }
      global$seledge = newe
    })
    if(!is.null(newe) && !is.null(oldedge) && newe != oldedge) global$seledge = paste(c(oldedge, newe), collapse = '\n')
  })

    get_seledge = reactive({
      print('get_seledge')
      print(global$seledge)
      global$seledge
    })

  ed_to_name = function(ed, plt) {
    dat = plt$x$data[[ed$curveNumber + 1]]
    type = ifelse(length(dat) == 14, 'text', ifelse('hovertext' %in% names(dat), 'edge', 'node'))
    pt = ed$pointNumber + 1

    if (type == 'edge') {
      name = dat$hovertext[pt]
    } else if (type == 'node') {
      col = ifelse(duplicated(dat$text)[pt], 5, 3)
      name = dat$text[pt]
    } else {
      name = dat$text[pt]
    }
    print('ed_to_name22')
    name
  }

  get_edgecols = reactive({
    plt = plotly_graph()
    edgecols = plt$x$data %>% map(~pluck(., 'line', 'color')) %>% discard(is.null)
    txt = plt$x$data[seq_len(length(edgecols))] %>% map(~pluck(., 'text'))
    tibble(edge = str_remove_all(unlist(txt), '\\.|>| '),
           col = rep(as.vector(edgecols), map_dbl(txt, length))) %>%
      unnest(col) %>% filter(!duplicated(.), !is.na(edge))
  })


  plotly_graph = reactive({

    print('plotly_graph')
    poplist = global$poplist
    pdat = get_pdat()
    req(pdat, poplist)
    nadmix = numadmix(graph_from_edgelist(as.matrix(global$edges)[,1:2]))
    textsize = 2.5
    print('plotly_graph 2')

    eg = pdat$eg
    nodelabs = imap(poplist, ~paste(.y, paste(.x, collapse = '\n'), sep = '\n')) %>%
      unlist %>% enframe(name = 'name', value = 'text')
    #nodelabs = poplist %>% group_by(pop) %>%
    #  summarize(text = paste(pop[1], paste(ind, collapse='\n'), sep = '\n'))
    nodes = pdat$nodes %>% left_join(nodelabs, by = c('name'))
    allnodes = eg %>% transmute(from, x, y, xend=0, yend=0, to=0) %>% filter(!duplicated(.))
    segtext = "ifelse(indegree == 1, to, from)"
    segtext = "paste(from, to, sep = ' -> ')"
    print('plotly_graph 3')
    sr = global$qpgraph_scorerange
    scoretext = ifelse(is.null(sr), round(global$score, 2),
                       paste0(round(sr['mid'], 2), '\n[', round(sr['lo'], 0), ' - ', round(sr['hi'], 0), ']'))

    gg = eg %>% mutate(rownum = 1:n()) %>%
      ggplot(aes(x=x, xend=xend, y=y, yend=yend, from=from, to=to)) +
      geom_segment(aes_string(linetype = 'type', col = 'as.factor(y)', text = segtext),
                   arrow=arrow(type = 'closed', angle = 10, length=unit(0.15, 'inches'))) +
      geom_text(aes(x = (x+xend)/2, y = (y+yend)/2, label = label, text = paste(from, to, sep = ' -> ')),
                size = textsize) +
      geom_text(data = nodes, aes_string(label = 'name', col = 'as.factor(yend)',
                                         from = NA, text = 'text'), size = textsize) +
      geom_point(data = allnodes, aes(x, y, text = from), col = 'black', alpha = 0) +
      annotate('text', x = min(nodes$x), y = max(nodes$yend), hjust = 0,
                label = paste0('score: ', scoretext, '\nadmix: ', nadmix)) +
      theme(panel.background = element_blank(),
            axis.line = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            legend.position = 'none') +
      xlab('') + ylab('') +
      scale_linetype_manual(values = c(admix=3, normal=1)) +
      ggtitle('') +
      scale_x_continuous(expand = c(0.1, 0.1))

    print('plotly_graph 4')
    plt = plotly::ggplotly(gg, source = 'src', tooltip=c('text')) #%>% event_register('plotly_click')
    print('plotly_graph 5')
    plt
  })

  output$graphplot <- renderPlotly({

    print('graphplot')
    suppressWarnings({withProgress(message = 'Processing...', {plt = plotly_graph()})})
    #event_register(plt, 'plotly_hover')
    print('graphplot2')
    plt
  })

  get_f4 = reactive({

    print('get f4')
    input$f4_run
    f2blocks = get_f2blocks()
    p1 = input$f4pops1
    p2 = input$f4pops2
    p3 = input$f4pops3
    p4 = input$f4pops4
    req(f2blocks, p1, p2, p3, p4)
    f4(f2blocks, p1, p2, p3, p4)
  })

  get_qpadm = reactive({

    print('get qpadm')
    input$qpadm_run
    f2blocks = get_f2blocks()
    p1 = input$qpadmpops1
    p2 = input$qpadmpops2
    p3 = input$qpadmpops3
    req(f2blocks, p1, p2, p3)
    qpadm(f2blocks, p1, p2, p3)
  })

  plotly_f2 = reactive({

    print('get f2')
    f2blocks = get_f2blocks()

    pdat = f2blocks %>%
      apply(1:2, mean, na.rm = T) %>%
      as_tibble(rownames = 'pop1') %>%
      pivot_longer(-pop1, 'pop2', values_to = 'est')

    gg = pdat %>%
      ggplot(aes(pop1, pop2, text = paste(pop1, pop2, round(est, 3), sep='\n'))) +
      geom_tile(aes(fill = est)) +
      xlab('') +
      ylab('') +
      theme(axis.text.x = element_text(angle = 90, hjust = 1), panel.background = element_blank()) +
      scale_fill_distiller(palette = 'Purples')

    plt = plotly::ggplotly(gg, source = 'src_f2', tooltip=c('text'))
    plt
  })

  plotly_f2_indiv = reactive({

    f2dat = get_f2dat_indiv()
    poplist = global$poplist
    ind = enframe(poplist, 'pop', 'ind') %>% unnest(ind)

    dat = f2dat %>%
      transmute(ind1 = pop1, ind2 = pop2, f2corr = f2, f2uncorr) %>%
      left_join(ind %>% transmute(ind1 = ind, pop1 = pop), by = 'ind1') %>%
      left_join(ind %>% transmute(ind2 = ind, pop2 = pop), by = 'ind2') %>%
      arrange(pop1, pop2) %>%
      mutate(i1 = as.numeric(factor(ind1, levels = unique(ind1[order(pop1)]))),
             i2 = as.numeric(factor(ind2, levels = unique(ind2[order(pop2)]))),
             f2 = ifelse(i1 > i2, f2corr, f2uncorr))
    d2 = enframe(sapply(poplist, length)) %>%
      arrange(name) %>%
      mutate(end = cumsum(value), start = lag(end, default = 0), mn = start + value/2 + 0.5)

    exp = 0
    gg = dat %>%
      ggplot(aes(i1, i2)) +
      geom_tile(aes(fill = f2, text = paste(ind1, ind2, round(f2, 3), sep = '\n'))) +
      geom_vline(aes(xintercept = end+0.5), data = d2 %>% head(-1) %>% rename(ind2 = name)) +
      geom_hline(aes(yintercept = end+0.5), data = d2 %>% head(-1) %>% rename(ind2 = name)) +
      scale_x_continuous(labels = d2$name, breaks = d2$mn, expand = c(exp, exp)) +
      scale_y_continuous(labels = d2$name, breaks = d2$mn, expand = c(exp, exp)) +
      xlab('') + ylab('') +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
            panel.background = element_blank()) +
      scale_fill_distiller(palette = 'Purples')
    plt = plotly::ggplotly(gg, source = 'src_f2', tooltip = 'text')
    plt
  })

  plotly_f4 = reactive({

    print('get f4')
    input$f4_run
    #isolate({
      f2blocks = get_f2blocks()
      p1 = input$f4pops1[1]
      p2 = input$f4pops2[1:2]
      f4dat = get_f4() %>% filter(pop1 == p1, pop2 %in% p2)
      req(f4dat)
    #})
    f4dat %<>%
        mutate(pop2 = recode(pop2, !!sym(p2[1]):='x', !!sym(p2[2]):='y')) %>%
        pivot_wider(names_from = pop2, values_from = c(est, se, z, p))
    xmin = min(f4dat$est_x - f4dat$se_x)
    xmax = max(f4dat$est_x + f4dat$se_x)
    ymin = min(f4dat$est_y - f4dat$se_y)
    ymax = max(f4dat$est_y + f4dat$se_y)
    gg = f4dat %>%
      ggplot(aes(est_x, est_y, col = pop3, text = paste0('f4(', p1, ', X; ', pop3, ', ', pop4, ')'))) +
      geom_smooth(method = 'lm', se = FALSE, formula=y~x-1) +
      geom_hline(yintercept = 0) +
      geom_vline(xintercept = 0) +
      geom_point() +
      theme(panel.background = element_blank()) +
      xlab(paste0('f4(', p1, ', ', p2[1], '; Pop 3, Pop 4)' )) +
      ylab(paste0('f4(', p1, ', ', p2[2], '; Pop 3, Pop 4)' )) +
      scale_x_continuous(limits = c(xmin, xmax)) +
      scale_y_continuous(limits = c(ymin, ymax))
    if(input$f4_se) {
      gg = gg +
        geom_errorbar(aes(ymin = est_y - se_y, ymax = est_y + se_y), width = 0) +
        geom_errorbarh(aes(xmin = est_x - se_x, xmax = est_x + se_x), height = 0)
    }
    plt = plotly::ggplotly(gg, source = 'src_f4', tooltip=c('text'))
  })

  output$f2heatmap = renderPlotly({

    print('f2plot')
    plt = plotly_f2()
    plt
  })

  output$f2heatmap_indiv = renderPlotly({

    print('f2plot indiv')
    plt = plotly_f2_indiv()
    plt
  })

  output$f4plot = renderPlotly({

    print('f4plot')
    plt = plotly_f4()
    plt
  })

  output$qpadmcomparison = renderPlot({

    print('qpadmcomparison')

    withProgress(message = 'running qpadm in R...', {
      qpadm_R = get_qpadm()
    })

    original = global$qpadmout
    gg = plot_comparison(original, qpadm_R)
    #plotly::ggplotly(gg, tooltip = 'label')
    gg
  })
  output$graphcomparison = renderPlotly({

    print('graphcomparison')

    withProgress(message = 'running qpgraph in R...', {
      qpgraph_R = get_fit()
    })

    original = global$qpgraphout
    gg = plot_comparison(original, qpgraph_R)
    plotly::ggplotly(gg, tooltip = 'label')
  })



  output$downloadgraph = downloadHandler('graph.tsv', function(file) {
    edges = get_fit()$edges
    if(input$downloadgraph_format == 'admixtools') edges %>% fit_to_qpgraph_format %>% write(file)
    else if(input$downloadgraph_format == 'edgelist') edges %>% write_tsv(file)
    else stop()
  })

  #dtfun = function(x) DT::renderDataTable({format_table(x)}, extensions = 'Buttons', options = do)
  #dtfun2 = function(x) reactive({DT::renderDataTable({format_table(x)}, extensions = 'Buttons', options = do)})
  ex = c('Buttons', 'Scroller')
  dtfun = function(...) DT::renderDataTable(..., extensions = ex, options = do, selection = 'single')

  output$f2 = dtfun({format_table(get_fit()$f2)})
  output$f3 = dtfun({format_table(get_fit()$f3)})
  output$opt = dtfun({format_table(get_fit()$opt)})
  output$minus1 = dtfun({format_table(get_minus1())})
  output$minusplus = dtfun({format_table(get_minusplus())})
  output$plus1 = dtfun(format_table(get_plus1()))
  output$decomposed = dtfun({format_table(get_decomposed())})
  output$treeneighbors = dtfun({format_table(get_treeneighbors())})
  output$flipadmix = dtfun({format_table(get_flipadmix())})
  output$addleaf = dtfun({format_table(get_addleaf())})

  output$initresample = dtfun({format_table(get_initresample())})
  output$indresample = dtfun({format_table(get_indresample())})
  output$snpresample = dtfun({format_table(get_snpresample())})

  output$f4stats = dtfun({format_table(get_f4())})

  output$qpadm_weights = dtfun({format_table(get_qpadm()$weights %>% rename(source = left))})
  output$qpadm_f4 = dtfun({format_table(get_qpadm()$f4)})
  output$qpadm_rankdrop = dtfun({format_table(get_qpadm()$rankdrop)})
  output$qpadm_popdrop = dtfun({format_table(get_qpadm()$popdrop)})


  onBookmark(function(state) {
    names(global) %>% map(~{state$values[[.]] = global[[.]]})
  })
  onRestore(function(state) {
    names(state$values) %>% map(~{global[[.]] = state$values[[.]]})
  })

  #setBookmarkExclude('sidebarItemExpanded')
  #setBookmarkExclude('navbar')

  observeEvent(input$tooltips, {
    print('refresh tooltips')
    if(input$tooltips %% 2 == 1) {
      updateActionButton(session, 'tooltips', icon = icon('exclamation'))
      imap(tt, ~removeTooltip(session, .y))
      imap(tt, ~addTooltip(session, .y, .x))
    } else {
      updateActionButton(session, 'tooltips', icon = icon('question'))
      imap(tt, ~removeTooltip(session, .y))
    }
  })
}

shinyApp(ui = ui, server = server, enableBookmarking = 'server')

#runApp(list(ui = ui, server = server), display.mode = 'showcase')


