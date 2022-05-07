#' Compare two models
#'
#' Plot a comparison of two runs of qpgraph or two runs of qpadm
#' @export
#' @param out1 First model
#' @param out2 Second model
#' @param name1 Optional name for first model
#' @param name2 Optional name for second model
#' @return A ggplot object.
#' @examples
#' fit1 = qpgraph(example_f2_blocks, example_graph, lsqmode = FALSE)
#' fit2 = qpgraph(example_f2_blocks, example_graph, lsqmode = TRUE)
#' plot_comparison(fit1, fit2)
plot_comparison = function(out1, out2, name1 = NULL, name2 = NULL) {
  # plots a comparison of two qpgraph or qpadm outputs
  if(is.null(name1)) name1 = enexpr(out1)
  if(is.null(name2)) name2 = enexpr(out2)

  f = ifelse('score' %in% names(out1), plot_comparison_qpgraph, plot_comparison_qpadm)
  f(out1, out2, name1, name2)
}

#' Compare two models
#'
#' Plot a comparison of two runs of qpgraph or two runs of qpadm
#' @export
#' @param out1 First model
#' @param out2 Second model
#' @param name1 Optional name for first model
#' @param name2 Optional name for second model
#' @return A plotly object.
#' @examples
#' fit1 = qpgraph(example_f2_blocks, example_graph, lsqmode = FALSE)
#' fit2 = qpgraph(example_f2_blocks, example_graph, lsqmode = TRUE)
#' plotly_comparison(fit1, fit2)
plotly_comparison = function(out1, out2, name1 = NULL, name2 = NULL) {

  if(is.null(name1)) name1 = enexpr(out1)
  if(is.null(name2)) name2 = enexpr(out2)
  p = plot_comparison(out1, out2, name1 = name1, name2 = name2)
  plotly::ggplotly(p)
}


plot_comparison_qpgraph = function(out1, out2, name1 = NULL, name2 = NULL) {

  f2comp = f3comp = NULL
  if(!is.null(out1$f2) && !is.null(out2$f2)) {

    show = c('pop1', 'pop2', 'est', 'fit', 'diff', 'se', 'z')
    common = intersect(names(out1$f2), names(out2$f2))
    f2comp = out1$f2 %>% select(all_of(intersect(show, common))) %>% mutate(i = 'x') %>%
      bind_rows(out2$f2 %>% select(all_of(intersect(show, common))) %>% mutate(i = 'y')) %>%
      gather('type', 'v', -c(pop1, pop2, i)) %>% spread(i, v) %>%
      rename(from = pop1, to = pop2) %>%
      mutate(type = paste0('f2_', type))
  }

  if(!is.null(out1$f3) && !is.null(out2$f3)) {

    show = c('pop2', 'pop3', 'est', 'fit', 'diff', 'se', 'z')
    common = intersect(names(out1$f3), names(out2$f3))
    f3comp = out1$f3 %>% select(all_of(intersect(show, common))) %>% mutate(i = 'x') %>%
      bind_rows(out2$f3 %>% select(all_of(intersect(show, common))) %>% mutate(i = 'y')) %>%
      gather('type', 'v', -c(pop2, pop3, i)) %>% spread(i, v) %>%
      rename(from = pop2, to = pop3) %>% filter(!is.na(x), !is.na(y)) %>%
      mutate(type = paste0('f3_', type))
  }

  if(!'low' %in% names(out1$edges)) out1$edges %<>% mutate(low = NA)
  if(!'low' %in% names(out2$edges)) out2$edges %<>% mutate(low = NA)
  if(!'high' %in% names(out1$edges)) out1$edges %<>% mutate(high = NA)
  if(!'high' %in% names(out2$edges)) out2$edges %<>% mutate(high = NA)

  out1$edges %>% rename(x = weight, xmin = low, xmax = high) %>%
    left_join(out2$edges %>% select(-'type') %>%
                rename(y = weight, ymin = low, ymax = high), by=c('from', 'to')) %>%
    bind_rows(f2comp, f3comp) %>%
    mutate_at(vars(c(xmin, xmax)), ~ifelse(is.na(.), x, .)) %>%
    mutate_at(vars(c(ymin, ymax)), ~ifelse(is.na(.), y, .)) %>%
    mutate(e = paste(from, to, sep = ' -> ')) %>%
    ggplot(aes(x=x, y=y, label = e)) +
    geom_errorbarh(aes(xmin = replace_na(xmin, 0), xmax = replace_na(xmax, 0)), height = 0, col = 'grey') +
    geom_errorbar(aes(ymin = replace_na(ymin, 0), ymax = replace_na(ymax, 0)), width = 0, col = 'grey') +
    geom_point(aes(col = from == to)) +
    facet_wrap('type', scales='free', dir = 'v', ncol = 3) +
    geom_abline() +
    xlab(paste0(name1, ' (score: ', round(out1$score, 2),')')) +
    ylab(paste0(name2, ' (score: ', round(out2$score, 2),')')) +
    theme(panel.background = element_blank(), axis.line = element_line(), legend.position = 'none') +
    geom_smooth(method='lm', se = FALSE, formula=y~x-1)
}


#' Plot an admixture graph
#' @export
#' @param graph An admixture graph. If it's an edge list with a \code{label} column,
#' those values will displayed on the edges
#' @param fix If \code{TRUE}, there will be an attempt to rearrange the nodes to minimize
#' the number of intersecting edges. This can take very long for large graphs.
#' By default this is only done for graphs with fewer than 10 leaves.
#' @param title Plot title
#' @param color Plot it in color or greyscale
#' @param textsize Size of edge and node labels
#' @param highlight_unidentifiable Highlight unidentifiable edges in red. Can be slow for large graphs. See \code{\link{unidentifiable_edges}}.
#' @param pos Optional data frame with node coordinates (columns `node`, `x`, `y`)
#' @param dates Optional named vector with dates for each node to plot dates on the y-axis
#' @param neff Optional named vector with effective population sizes for each node to show population sizes
#' @param scale_y If TRUE, scale the y-axis according to \code{dates} vector. The default is `FALSE`.
#' @param hide_weights Boolean value specifying if the drift weights that are equal to 1 should be hidden. The default is `FALSE`
#' @return A ggplot object
#' @examples
#' plot_graph(example_graph)
#' 
#' # Plot a random simulation output. Show dates and population sizes on the plot
#' out = random_sim(nleaf=5, nadmix=1)
#' plot_graph(out$edges, dates=out$dates, neff=out$neff)
plot_graph = function(graph, fix = NULL, title = '', color = TRUE, textsize = 2.5, highlight_unidentifiable = FALSE,
                      pos = NULL, dates = NULL, neff = NULL, scale_y = FALSE, hide_weights = FALSE) {
  
  pdat = graph_to_plotdat(graph, fix = fix, fix_down = TRUE, pos = pos)
  
  # ## Correct y-axis values if leaf nodes are fixed to 0 (i.e., if fix_leaf=TRUE)
  # if (sum(dates[pdat$nodes$name]) == 0 && !is.null(dates)){
  #   if ("igraph" %in% class(graph)) g = graph
  #   else g = edges_to_igraph(graph)
  #   pdat$nodes %<>%
  #     mutate(diff = y-min(y),
  #            xdirection = ifelse(xend-x<0, 'pos', 'neg'),
  #            x = ifelse(xdirection=='pos', x+diff/2, x-diff/2),
  #            y = min(y)) %>%
  #     select(-diff, -xdirection)
  #   
  #   nodes = rename(pdat$nodes, x=xend, y=yend, xend=x, yend=y, to=name)
  #   
  #   for (i in nodes$to){
  #     pdat$eg[pdat$eg$to == i, c("x", "y", "xend", "yend", "to")] = nodes[nodes$to == i, c("x", "y", "xend", "yend", "to")]
  #   }
  # }
  
  pdat[['comb']] = bind_rows(select(pdat$nodes, name, x, y), 
                             select(pdat$eg, name, x, y)) %>% distinct()
  
  layers = vector('list')
  if (!is.null(dates)){
    if (!all(sort(unique(pdat$comb$name)) == sort(unique(names(dates))))) stop("dates must be a named vector with a number for each node!")
    
    date_df = pdat$comb %>% left_join(tibble(name=names(dates), date=dates), by="name") %>% 
      group_by(y) %>% mutate(date = max(date)) %>%  # To fix the mismatch when fix_leaf=TRUE
      select(y, date) %>% distinct()
    
    if (isTRUE(scale_y)){
      scaled_y = setNames(round(date_df$date/mean(date_df$date), 4), as.character(date_df$y))
      date_df %<>%
        mutate(y = scaled_y[as.character(y)])
      pdat$eg %<>%
        mutate(y = scaled_y[as.character(y)],
               yend = scaled_y[as.character(yend)])
      pdat$nodes %<>%
        mutate(y = scaled_y[as.character(y)],
               yend = scaled_y[as.character(yend)])
      pdat$comb %<>%
        mutate(y = scaled_y[as.character(y)])
    }
    
    layers$date = list(geom_hline(yintercept = date_df$y, linetype=3, alpha=0.6),
                       scale_y_continuous(breaks=date_df$y, labels=date_df$date),
                       ylab('Time (in generation ago)'),
                       theme(panel.background = element_blank(),
                             axis.line = element_blank(),
                             axis.text.x = element_blank(),
                             axis.ticks = element_blank(),
                             legend.position = 'none'))
    
  } else{
    layers$date = list(ylab(''),
                       theme(panel.background = element_blank(),
                             axis.line = element_blank(),
                             axis.text = element_blank(),
                             axis.ticks = element_blank(),
                             legend.position = 'none'))
  }
  if (!is.null(neff)){
    if (!all(sort(unique(pdat$comb$name)) == sort(unique(names(neff))))) stop("neff must be a named vector with a number for each node!")
    neff_df = pdat$eg %>%
      mutate(x = (x+xend)/2,
             y = (y+yend)/2) %>%
      select(-c(name, label)) %>%
      rename(name = to) %>%
      left_join(tibble(name=names(neff), label=neff), by="name") %>%
      select(-c(name, xend, yend, type))
    
    if ("R" %in% pdat$eg$name ){
      rootlab = pdat$eg %>%
        filter(name == "R") %>%
        select(name, x, y) %>%
        distinct() %>%
        mutate(label = neff['R']) %>%
        select(-name)
      
      neff_df %<>% bind_rows(rootlab)
    }
    layers$neff = geom_text(data=neff_df, aes_string(x='x', y='y', label='label'), col='black', inherit.aes = FALSE, alpha=0.8, size = textsize - 0.5, nudge_x = .0, nudge_y = .15)
  }
  
  if (isTRUE(hide_weights)){
    pdat$eg %<>%
      mutate(label = ifelse(label==1000, " ", as.character(label)))
  }
  
  if(color) {
    layers$gs = geom_segment(aes_string(linetype = 'type', col = 'as.factor(y)'),
                             arrow=arrow(type = 'closed', angle = 10, length=unit(0.15, 'inches')))
    layers$gl = geom_label(data=pdat$nodes, aes_string(label = 'name', col='as.factor(yend)'), size = textsize)
  } else {
    # gs = geom_segment(aes_string(linetype = 'type', col = "as.factor((to %in% c('X1', 'X2') & name %in% c('AB1', 'AB2')) + (name == 'AB1' & to == 'X1'))"),
    #                   arrow=arrow(type = 'closed', angle = 10, length=unit(0.15, 'inches')))
    layers$gs = geom_segment(aes_string(linetype = 'type'),
                             arrow=arrow(type = 'closed', angle = 10, length=unit(0.15, 'inches')))
    layers$gl = geom_label(data=pdat$nodes, aes_string(label = 'name'), col = 'black', size = textsize)
  }
  layers$gl2 = if(is.null(pdat$internal)) NULL else geom_text(data=pdat$internal, aes_string(label = 'name'), size = textsize, nudge_x = -.1, nudge_y = .1)
  
  layers$aes = list(geom_text(aes(x = (x+xend)/2, y = (y+yend)/2, label = label), size = textsize, nudge_y = -.15),
                    xlab(''),
                    scale_linetype_manual(values = c(admix=3, normal=1)),
                    ggtitle(title),
                    scale_x_continuous(expand = c(0.15, 0.15)))
  
  if(highlight_unidentifiable) {
    if ("igraph" %in% class(graph)) g = graph
    else g = edges_to_igraph(graph)
    unid = unidentifiable_edges(g)
    unid2 = pdat$eg %>% rename(from = name) %>% right_join(unid %>% select(-type), by = c('from', 'to'))
    layers$highlight = geom_segment(aes_string(linetype = 'type'), col = 'red', size = 1, data = unid2,
                                    arrow=arrow(type = 'closed', angle = 10, length=unit(0.15, 'inches')))
  }
 
  ggplot(pdat$eg, aes(x=x, xend=xend, y=y, yend=yend)) + 
    layers 
}


graph_to_plotdat = function(graph, fix = NULL, fix_down = TRUE, pos = NULL) {

  if(class(graph)[1] == 'igraph') {
    graph = graph
    edges = graph %>% igraph::as_edgelist() %>% as_tibble(.name_repair = ~c('V1', 'V2'))
  } else {
    if(is.null(colnames(graph))) colnames(graph) = paste0('V', 1:ncol(graph))
    edges = graph %>% as_tibble
    names(edges)[1:2] = c('V1', 'V2')
    graph = igraph::graph_from_edgelist(as.matrix(graph)[,1:2])
  }
  leaves = get_leafnames(graph)
  #edges = as_tibble(edges, .name_repair = ~c('V1', 'V2'))
  admixnodes = unique(edges[[2]][edges[[2]] %in% names(which(table(edges[[2]]) > 1))])

  if(is.null(pos)) {
    pos = data.frame(names(V(graph)), igraph::layout_as_tree(graph), stringsAsFactors = F) %>%
      set_colnames(c('node', 'x', 'y'))
  }
  eg = graph %>% igraph::as_edgelist() %>%
    as_tibble(.name_repair = ~c('V1', 'V2')) %>%
    left_join(pos, by=c('V1'='node')) %>%
    left_join(pos %>% transmute(V2=node, xend=x, yend=y), by='V2') %>%
    mutate(type = ifelse(V2 %in% admixnodes, 'admix', 'normal')) %>% rename(name=V1, to=V2)

  if(isTRUE(fix) || is.null(fix) && length(leaves) < 10) eg = fix_layout(eg, graph)
  if(fix_down) eg = fix_shiftdown(eg, graph)

  if(!'label' %in% names(edges)) {
    if('weight' %in% names(edges)) {
      lab = ifelse(edges$type == 'admix', paste0(round(edges$weight*100), '%'), round(edges$weight*1000))
    }
    else lab = ''
    edges %<>% mutate(label = lab)
  }
  eg %<>% left_join(edges %>% transmute(name=V1, to=V2, label), by=c('name', 'to'))
  if(any(str_detect(eg$name, ' \\(100\\)'))) {
    internal = eg %>% filter(!to %in% leaves) %>%
      rename(x=xend, y=yend, xend=x, yend=y) %>% transmute(name = to, x, y, xend, yend) %>%
      mutate(name = str_replace_all(name, '.+\\(|\\)', ''))
  } else internal = NULL
  nodes = eg %>% filter(to %in% leaves) %>%
      rename(x=xend, y=yend, xend=x, yend=y) %>% transmute(name = to, x, y, xend, yend)
  namedList(eg, nodes, internal)
}


fix_layout = function(coord, graph) {
  # rearranges nodes in tree layout, so that there are fewer long branches
  dst = igraph::distances(graph, V(graph)[1], mode='out')[1,]
  totdist = function(x, xend) {
    sum(abs(x-xend), na.rm=TRUE)
  }
  swap_subtree = function(graph, coord, node) {
    center = coord$x[match(node, coord$name)[1]]
    offspring = names(subcomponent(graph, node, mode='out'))
    dst2 = igraph::distances(graph, node, mode='out')[1,]
    for(n in offspring) {
      if(dst[node] + dst2[n] <= dst[n])
        coord %<>% mutate(x = ifelse(name == n, 2*center-x, x),
                          xend = ifelse(to == n, 2*center-xend, xend))
    }
    coord
  }
  swap_subtree2 = function(graph, coord, node) {
    center = coord$x[match(node, coord$name)[1]]
    children = names(neighbors(graph, node, mode='out'))
    dst2 = igraph::distances(graph, node, mode='out')[1,]
    if(length(children) == 2) {
      for(child in children) {
        offspring = names(subcomponent(graph, child, mode='out'))
        xx = coord %>% filter(name == node, to == child) %$% xend
        shift = 2 * (xx - center)
        for(n in offspring) {
          if(dst[node] + dst2[n] <= dst[n])
            coord %<>% mutate(x = ifelse(name == n, x-shift, x),
                              xend = ifelse(to == n, xend-shift, xend))
        }
      }
    }
    coord
  }
  inner_bf = setdiff(names(subcomponent(graph, V(graph)[1])), get_leafnames(graph))
  for(node in rep(inner_bf, 1)) {

    d1 = totdist(coord$x, coord$xend)
    n1 = num_intersecting(coord)
    s1 = num_samepos(coord)
    c2 = swap_subtree(graph, coord, node)
    d2 = totdist(c2$x, c2$xend)
    n2 = num_intersecting(c2)
    s2 = num_samepos(c2)
    if((n2 < n1 | n2 == n1 & d2 < d1) & s2 <= s1) coord = c2
  }

  for(node in rep(inner_bf, 1)) {
    #if(node == 'Rrrlx') browser()
    d1 = totdist(coord$x, coord$xend)
    n1 = num_intersecting(coord)
    s1 = num_samepos(coord)
    c2 = swap_subtree2(graph, coord, node)
    d2 = totdist(c2$x, c2$xend)
    n2 = num_intersecting(c2)
    s2 = num_samepos(c2)
    if((n2 < n1 | n2 == n1 & d2 < d1) & s2 <= s1) coord = c2
  }
  #fix_layout2(coord)
  coord
}

fix_layout2 = function(coord) {
  # rearranges trees so that horizontal edges on same level don't hide each other
  edges = coord %>% filter(!is.na(to))
  leaves = coord %>% filter(is.na(to))

  nodeshift = edges %>% filter(y == yend)
  if(nrow(nodeshift) > 0) {
    nodeshift %<>% group_by(g = interaction(y, yend, cumsum(cummax(lag(xend, default = first(xend))) < x))) %>%
      mutate(tot = n(), this = (seq_len(tot[1]))-1, shift = this/tot) %>% filter(shift != 0) %>%
      dplyr::select(g, name, to, shift) %>% ungroup %>% gather(irl, node, name, to) %>% dplyr::select(node, shift)
    if(nrow(nodeshift) > 0) {
      edges %<>% left_join(nodeshift %>% transmute(name=node, s1 = shift), by='name') %>%
        left_join(nodeshift %>% transmute(to=node, s2 = shift), by='to') %>%
        mutate(y = ifelse(is.na(s1), y, y+s1), yend = ifelse(is.na(s2), yend, yend+s2))
      leaves %<>% left_join(nodeshift %>% transmute(name=node, shift), by='name') %>%
        mutate(y=ifelse(is.na(shift), y, y+shift))
    }
  }
  bind_rows(edges, leaves)
}

fix_shiftdown = function(coord, graph) {
  adm = find_admixedges(graph) %>%
    left_join(coord, by = c('from'='name', 'to'='to')) %>%
    filter(yend >= y) %>% mutate(by = yend - y + 1)
  for(i in seq_len(nrow(adm))) {
    by = coord %>% filter(name == adm$from[i], to == adm$to[i]) %$% {yend - y + 1}
    coord = shift_subgraph_down(graph, coord, adm$to[i], by = by)
  }
  coord
}

shift_subgraph_down = function(graph, coord, node, by) {
  # returns new coords with some nodes shifted down
  if(by <= 0) return(coord)
  nodes = graph %>% subcomponent(node, mode = 'out') %>% names
  coord %>% mutate(shift = name %in% nodes,
                   shift1 = to %in% nodes,
                   y = ifelse(shift, y-by, y),
                   yend = ifelse(shift | shift1, yend-by, yend)) %>%
    select(-shift, -shift1)
}

intersecting = function(x1, y1, x2, y2, x3, y3, x4, y4) {
  # returns TRUE iff segment 12 intersects with segment 34
  # extends segments by small amounts, jitters x3, x4 to prevent singularity
  eps = 1e-3
  xd1 = x1-x2
  yd1 = y1-y2
  xd2 = x3-x4
  yd2 = y3-y4

  x1 = x1 + xd1*eps
  x2 = x2 - xd1*eps
  y1 = y1 + yd1*eps
  y2 = y2 - yd1*eps
  x3 = x3 + xd2*eps
  x4 = x4 - xd2*eps
  y3 = y3 + yd2*eps
  y4 = y4 - yd2*eps

  if(x1 == x2 & x3 == x4) return(x1 == x3 & max(y1, y2) > min(y3, y4) & min(y1, y2) < max(y3, y4))

  if(x1 == x2) {
    sl = (y4-y3)/(x4-x3)
    inters = max(y1, y2) > min(y3, y4) & min(y1, y2) < max(y3, y4)
    inters = inters & between(x1, min(x3, x4), max(x3, x4))
    inters = inters & between(sl * x1 + (y3 - sl*x3), min(y1, y2), max(y1, y2))
    return(inters)
  }
  if(x3 == x4) {
    sl = (y2-y1)/(x2-x1)
    inters = max(y3, y4) > min(y1, y2) & min(y3, y4) < max(y1, y2)
    inters = inters & between(x3, min(x1, x2), max(x1, x2))
    inters = inters & between(sl * x3 + (y1 - sl*x1), min(y3, y4), max(y3, y4))
    return(inters)
  }

  A1 = (y1-y2)/(x1-x2)
  A2 = (y3-y4)/(x3-x4)
  b1 = y1 - A1*x1
  b2 = y3 - A2*x3
  denom = A1 - A2
  Xa = (b2 - b1) / denom
  isTRUE(Xa > max(min(x1,x2), min(x3,x4)) & Xa < min(max(x1,x2), max(x3,x4)))
}

num_intersecting = function(dat) {
  # counts the number of intersecting segments
  # dat has x, xend, y, yend
  pairs = combn(seq_len(nrow(dat)), 2)
  fun = function(r1, r2) {dat %$% {intersecting(x[r1], y[r1], xend[r1], yend[r1], x[r2], y[r2], xend[r2], yend[r2]) & (length(unique(c(name[r1], to[r1], name[r2], to[r2]))) == 4)}}
  sum(map2_lgl(pairs[1,], pairs[2,], fun))

}

num_samepos = function(dat) {
  # counts the number of nodes at the same position
  dat %<>% select(to, xend, yend) %>% distinct
  nnodes = dat %>% nrow
  dat %<>% select(-to) %>% distinct
  samepos = nnodes - nrow(dat)
  samepos
}


collapse_edges = function(edges, below=0) {
  # input: graph edgelist
  # requires column 'weight'
  # collapses nodes which are separated by 'below' or less
  names(edges)[1:2] = c('from', 'to')
  admixnodenames = edges %>% filter(type == 'admix') %$% unique(to)
  leaves = setdiff(edges[[2]], edges[[1]])
  excl = c(admixnodenames, leaves)
  adjmat = edges %>% select(1, 2, weight) %>% bind_rows(rename(., from=to, to=from)) %>%
    spread(to, weight, fill = 0) %>% column_to_rownames(var = 'from') %>% as.matrix
  adjmat[adjmat > below] = 0
  adjmat[rownames(adjmat) %in% excl, ] = 0
  adjmat[, colnames(adjmat) %in% excl] = 0
  diag(adjmat) = 1
  comp = igraph::components(graph_from_adjacency_matrix(as.matrix(adjmat), weighted='x'))
  newnam = comp$membership %>% enframe %>% group_by(value) %>%
    mutate(name, groupname = paste(name, collapse='.')) %>% ungroup %>% select(-value) %>% deframe
  edges %>% mutate(from = newnam[from], to = newnam[to]) %>% filter(from != to)
}


plot_graph_pcs = function(graph, pcs) {

  leafcoords = pcs %>% rename(lon=PC1, lat=PC2)
  leafcoords2 = leafcoords %>% group_by(group) %>% summarize(x = mean(lon), y = mean(lat))
  sg = graph %>% simplify_graph
  node_coord = node_coords_3d(sg, leafcoords2)
  edge_coord = sg %>% as_edgelist %>% as.data.frame(stringsAsFactors=F) %>%
    set_colnames(c('from', 'to')) %>% mutate(eid = 1:n()) %>% gather(type, node, -eid) %>%
    left_join(node_coord, by='node') %>% arrange(eid)
  edge_coord %<>% mutate(admix = eid %in% (filter(., type == 'to') %>% group_by(node) %>%
                                             mutate(cnt = n()) %>% filter(cnt == 2) %$% eid),
                         admix=ifelse(admix, 'dot', 'solid'))
  node_coord2 = leafcoords %>% rename(x=lon, y=lat) %>%
    bind_rows(node_coord %>% filter(z==0) %>% transmute(s = 'center', group=node, x, y)) %>%
    mutate(z=0, sz = ((s == 'center')*5+2))

  ax = list(visible=FALSE)

  plotly::plot_ly(type = 'scatter3d', mode = 'lines', hoverinfo = 'text') %>%
    plotly::add_trace(data = edge_coord %>% mutate(const='grey'), x = ~x, y = ~y, z = ~z,
              split = ~eid, line=list(color=~const, dash=~admix), showlegend=F) %>%
    plotly::add_markers(data = node_coord2, x = ~x, y = ~y, z = ~z,
                hovertext = ~group, marker=list(size=~((s=='center')*5+5)), color = ~group, colors='Set1') %>%
    graphics::layout(scene = list(xaxis=list(title='PC1'), yaxis=list(title='PC2'), zaxis=ax,
                                  camera = list(projection = list(type = 'orthographic'),
                                                eye = list(x = 0, y = 0, z = 1),
                                                up = list(x = 0, y = 1, z = 0))))
}


#' Plot an admixture graph on a map
#' @export
#' @param graph a two column matrix specifying the admixture graph. first column is source node,
#' second column is target node. the first edge has to be root -> outgroup.
#' admixture nodes are inferred as those nodes which are the targets of two different sources.
#' @param leafcoords data frame with columns \code{group}, \code{lon}, \code{lat}
#' @param shapedata shapedata
#' @return a plotly object
#' @examples
# #' \dontrun{
#' plot_graph_map(example_igraph, example_anno)
# #' }
plot_graph_map = function(graph, leafcoords, shapedata = NULL) {

  stopifnot(all(c('iid', 'group', 'lat', 'lon') %in% names(leafcoords)))
  leafcoords %<>% filter(group %in% get_leafnames(graph))
  leafcoords2 = leafcoords %>% group_by(group) %>% summarize(x = mean(lon), y = mean(lat))
  sg = graph %>% simplify_graph
  node_coord = node_coords_3d(sg, leafcoords2)
  edge_coord = sg %>% as_edgelist %>% as.data.frame(stringsAsFactors=F) %>%
    set_colnames(c('from', 'to')) %>% mutate(eid = 1:n()) %>% gather(type, node, -eid) %>%
    left_join(node_coord, by='node') %>% arrange(eid)
  edge_coord %<>% mutate(admix = eid %in%
                           (filter(., type == 'to') %>% group_by(node) %>%
                              mutate(cnt = n()) %>% filter(cnt == 2) %$% eid),
                         admix=ifelse(admix, 'dot', 'solid'))
  node_coord2 = leafcoords %>% rename(x=lon, y=lat) %>%
    bind_rows(node_coord %>% filter(z==0) %>% transmute(iid = 'center', group=node, x, y)) %>%
    mutate(z=0, sz = ((iid == 'center')*5+2))

  ax = list(visible=FALSE)
  if(is.null(shapedata)) shapedata = shapedata_basic
  dat = shapedata %>% filter(between(long, min(node_coord$x), max(node_coord$x)),
                             between(lat, min(node_coord$y), max(node_coord$y))) %>%
    mutate(z = 0, group=1) %>% mutate(const='black')

  plotly::plot_ly(dat, type = 'scatter3d', mode = 'lines', hoverinfo = 'text') %>%
    plotly::add_trace(x = ~long, y = ~lat, z = ~z, split = ~id, line=list(color=~const), showlegend=F) %>%
    plotly::add_trace(data = edge_coord %>% mutate(const='grey'), x = ~x, y = ~y, z = ~z,
              split = ~eid, line=list(color=~const, dash=~admix), showlegend=F) %>%
    plotly::add_markers(data = node_coord2, x = ~x, y = ~y, z = ~z,
                hovertext = ~group, marker=list(size=~((iid=='center')*5+5)), color = ~group, colors='Set1') %>%
    plotly::layout(scene = list(xaxis=ax, yaxis=ax, zaxis=ax,
                                  camera = list(projection = list(type = 'orthographic'),
                                                eye = list(x = 0, y = 0, z = 1),
                                                up = list(x = 0, y = 1, z = 0))))
}


# Plot an admixture graph on a map
# @param graph an admixture graph as an igraph object.
# @param leafcoords data frame with columns \code{group}, \code{lon}, \code{lat}
# @param map_layout 1 or 2
# @return a ggplot object.
# @examples
# \dontrun{
# plot_graph_map2(example_igraph, example_anno, 1)
# }
plot_graph_map2 = function(graph, leafcoords, map_layout = 1) {

  stopifnot(all(c('iid', 'group', 'lat', 'lon') %in% names(leafcoords)))
  leafcoords %<>% filter(group %in% get_leafnames(graph))
  leafcoords2 = leafcoords %>% group_by(group) %>% summarize(x = mean(lon), y = mean(lat))
  sg = graph %>% simplify_graph
  node_coord = node_coords_3d(sg, leafcoords2)
  edge_coord = sg %>% as_edgelist %>% as.data.frame(stringsAsFactors=F) %>%
    set_colnames(c('from', 'to')) %>% mutate(eid = 1:n()) %>% gather(type, node, -eid) %>%
    left_join(node_coord, by='node') %>% arrange(eid)
  node_coord2 = leafcoords %>% rename(x=lon, y=lat) %>%
    bind_rows(node_coord %>% filter(z==0) %>% transmute(iid = 'center', group=node, x, y)) %>%
    mutate(z=0, sz = ((iid == 'center')*5+2))

  ax = list(visible=FALSE)

  src1 = 'https://basemap.nationalmap.gov/arcgis/rest/services/USGSImageryOnly/MapServer/tile/{z}/{y}/{x}'
  if(map_layout == 1) lo = function(x) plotly::layout(x, mapbox = list(style = 'stamen-terrain', zoom = 2),
                                                      xaxis=ax, yaxis=ax)
  if(map_layout == 2) lo = function(x) plotly::layout(
    x, mapbox = list(style = 'white-bg', zoom = 2,
                     layers = list(list(below = 'traces', sourcetype = 'raster', source = list(src1)))))

  plotly::plot_ly(type = 'scattermapbox', mode='lines', hoverinfo = 'text') %>%
    plotly::add_trace(data = edge_coord %>% mutate(const='grey', lon=x, lat=y), lon = ~x, lat = ~y,
              split = ~eid, line=list(color=~const, width=1), showlegend=F) %>%
    plotly::add_trace(data = node_coord2, mode='markers', x = ~x, y = ~y,
              hovertext = ~group, marker=list(size=~((iid=='center')*5+5)), color = ~group, colors='Set1') %>%
    lo
}


#' Plot samples on a map
#' @export
#' @param leafcoords data frame with columns \code{iid}, \code{lon}, \code{lat}
#' @param map_layout 1 or 2
#' @param color color
#' @param colorscale colorscale
#' @param collog10 collog10
#' @return a plotly object.
#' @examples
#' \dontrun{
#' plot_map(example_anno, 1)
#' }
plot_map = function(leafcoords, map_layout = 1, color = 'yearsbp', colorscale = 'Portland', collog10 = TRUE) {

  stopifnot(all(c('iid', 'lat', 'lon', color) %in% names(leafcoords)))
  if(!color %in% names(leafcoords)) leafcoords %<>% mutate(color = 1)
  else leafcoords %<>% mutate(color = !!sym(color))
  leafcoords %<>% filter(!is.na(color), between(lat, -90, 90), between(lon, -180, 180))
  if(collog10) leafcoords %<>% filter(color > 0) %>% mutate(color = log10(color))

  ax = list(visible=FALSE)

  src1 = 'https://basemap.nationalmap.gov/arcgis/rest/services/USGSImageryOnly/MapServer/tile/{z}/{y}/{x}'

  if(map_layout == 1) lo = function(x) plotly::layout(
    x, mapbox = list(style = 'stamen-terrain', center = list(lat = 48.2, lon = 16.3), zoom = 2), xaxis=ax, yaxis=ax)

  if(map_layout == 2) lo = function(x) plotly::layout(
    x, mapbox = list(style = 'white-bg', center = list(lat = 48.2, lon = 16.3),
                     zoom = 2, layers = list(list(below = 'traces', sourcetype = 'raster', source = list(src1)))))

  plotly::plot_ly(leafcoords, type = 'scattermapbox', mode='markers',
                  marker=list(size=4, color = ~color, colorscale = colorscale, showscale = TRUE,
                              colorbar = list(thickness=15, outlinewidth = 0,
                                              title = paste0(ifelse(collog10, 'log10\n', ''), color))),
                  hoverinfo = 'text', x = ~lon, y = ~lat, hovertext = ~group) %>% lo
}


node_coords_3d = function(graph, leafcoords, rootlon=NA, rootlat=NA) {
  # given a graph and 2d coordinates of leaves,
  # this function returns 3d coords for all nodes
  root = names(V(graph)[1])
  leaves = get_leafnames(graph)
  #leafcoords %<>% filter(!is.na(x), !is.na(y))
  leafcoords %<>% mutate(x = ifelse(is.na(x), runif(1, min(x,na.rm=T), max(x,na.rm=T)), x),
                        y = ifelse(is.na(y), runif(1, min(y,na.rm=T), max(y,na.rm=T)), y))
  stopifnot(all(leaves %in% leafcoords$group))
  paths = graph %>% igraph::all_simple_paths(root, leaves, mode='out')
  out = leafcoords %>% mutate(z = 0) %>%
    bind_rows(tibble(group = root, x = mean(.$x, na.rm=T), y = mean(.$y, na.rm=T), z = 1))
  rootcoord = out %>% filter(group == root) %$% c(x=x, y=y, z=z)
  if(!is.na(rootlon) & !is.na(rootlat)) {
    rootcoord = c(x=rootlon, y=rootlat, z=1)
  }

  path_coords = function(path, startcoord, endcoord) {
    # given a sequence of nodes where the first and last node have 3d coords, return xyz for all intermediate nodes
    map2(startcoord, endcoord, ~seq(.x, .y, length.out=length(path))) %>%
      bind_cols %>% add_column(node=path, .before=1)
  }

  paths %>% map(names) %>%
    map(~path_coords(., rootcoord, filter(out, group == tail(., 1)) %$% c(x=x, y=y, z=z))) %>%
    bind_rows(.id='pid') %>% group_by(node) %>% summarize(x = mean(x), y = mean(y), z = mean(z))
}


make_favicon = function() {

  g = random_admixturegraph(str_sub('       ', 1, 1:6), 1)

  g %<>% simplify_graph
  graph = g
  edges = igraph::as_edgelist(g)

  edges = as_tibble(edges, .name_repair = ~c('V1', 'V2'))
  admixnodes = unique(edges[[2]][edges[[2]] %in% names(which(table(edges[[2]]) > 1))])

  pos = data.frame(names(V(graph)), igraph::layout_as_tree(graph), stringsAsFactors = F) %>%
    set_colnames(c('node', 'x', 'y'))
  eg = as_tibble(as_edgelist(graph)) %>% left_join(pos, by=c('V1'='node')) %>%
    left_join(pos %>% transmute(V2=node, xend=x, yend=y), by='V2') %>%
    mutate(type = ifelse(V2 %in% admixnodes, 'admix', 'normal')) %>% rename(name=V1, to=V2)

  eg = fix_layout(eg, graph)
  edges %<>% mutate(label='')
  eg %<>% left_join(edges %>% transmute(name=V1, to=V2, label), by=c('name', 'to'))
  nodes = eg %>% filter(to %in% get_leafnames(graph)) %>% rename(x=xend, y=yend, xend=x, yend=y) %>%
    transmute(name = to, x, y, xend, yend)

  eg %>%
    ggplot(aes(x=x, xend=xend, y=y, yend=yend)) +
    geom_segment(aes_string(col = 'as.factor(y)'), size=12) +
    # arrow=arrow(type = 'closed', angle = 20, length=unit(0.15, 'inches'))
    theme(panel.background = element_rect(fill = "transparent"),
          rect = element_rect(fill = "transparent"),
          axis.line = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.position = 'none', plot.margin=grid::unit(c(0,0,0,0), "mm")) +
    xlab('') + ylab('') + scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0,0))

  ggsave('logo.svg', width=4, height=4, bg = "transparent")

  #pkgdown::build_favicons(overwrite = T)
  #system('mv pkgdown/favicon/* docs/')
}


plot_graph_interactive = function(graph, fix = NULL, title = '', color = TRUE) {

  if(class(graph)[1] == 'igraph') {
    graph = graph
    edges = as_edgelist(graph) %>% as_tibble(.name_repair = ~c('V1', 'V2'))
  } else {
    graph = igraph::graph_from_edgelist(as.matrix(graph)[,1:2])
    edges = graph %>% as_tibble
    names(edges)[1:2] = c('V1', 'V2')
  }
  #edges = as_tibble(edges, .name_repair = ~c('V1', 'V2'))
  admixnodes = unique(edges[[2]][edges[[2]] %in% names(which(table(edges[[2]]) > 1))])

  pos = data.frame(names(V(graph)), igraph::layout_as_tree(graph), stringsAsFactors = F) %>%
    set_colnames(c('node', 'x', 'y'))
  eg = as_tibble(as_edgelist(graph)) %>% left_join(pos, by=c('V1'='node')) %>%
    left_join(pos %>% transmute(V2=node, xend=x, yend=y), by='V2') %>%
    mutate(type = ifelse(V2 %in% admixnodes, 'admix', 'normal')) %>% rename(name=V1, to=V2)

  if(isTRUE(fix) || is.null(fix) && length(V(graph)) < 10) eg = admixtools:::fix_layout(eg, graph)

  if(!'label' %in% names(edges)) edges %<>% mutate(label='')
  eg %<>% left_join(edges %>% transmute(name=V1, to=V2, label), by=c('name', 'to'))
  nodes = eg %>% filter(to %in% get_leafnames(graph)) %>% rename(x=xend, y=yend, xend=x, yend=y) %>%
    transmute(name = to, x, y, xend, yend, to=NA, rownum = 1:n())

  plt = eg %>% mutate(rownum = 1:n()) %>%
    ggplot(aes(x=x, xend=xend, y=y, yend=yend, name=name, to=to, rownum=rownum)) +
    geom_segment(aes_string(linetype = 'type', col = 'as.factor(y)'),
                 arrow=arrow(type = 'closed', angle = 10, length=unit(0.15, 'inches'))) +
    geom_text(aes(x = (x+xend)/2, y = (y+yend)/2, label = label)) +
    geom_text(data=nodes, aes_string(label = 'name', col='as.factor(yend)', rownum='rownum'), size=3) +
    theme(panel.background = element_blank(),
          axis.line = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.position = 'none') +
    xlab('') + ylab('') +
    scale_linetype_manual(values = c(admix=3, normal=1)) +
    ggtitle(title) +
    scale_x_continuous(expand = c(0.15, 0.15))

  plt
}


plot_comparison_qpadm = function(out1, out2, name1 = NULL, name2 = NULL) {

  yscale = diff(range(out2$weights$weight))
  p1 = out1$weights %>% left_join(out2$weights, by = c('target', 'left')) %>%
    ggplot(aes(weight.x, weight.y)) +
    geom_rect(xmin = 0, xmax = 1, ymin = 0, ymax = 1, fill = 'blue', alpha = 0.05) +
    geom_point() +
    geom_abline() +
    geom_errorbar(aes(ymin = weight.y - se.y, ymax = weight.y + se.y), width=0) +
    geom_errorbarh(aes(xmin = weight.x - se.x, xmax = weight.x + se.x), height=0) +
    geom_text(aes(label = left), nudge_y = yscale/10) +
    theme(panel.background = element_blank(), axis.line = element_line(), legend.position = 'none') +
    xlab(paste0(name1, ' weights')) +
    ylab(paste0(name2, ' weights'))

  if(nrow(out1$popdrop) < 2 || nrow(out2$popdrop) < 2) return(p1)

  nams = c('pat', 'wt', 'dof', 'chisq', 'p', 'f4rank', 'feasible', 'best', 'dofdiff', 'chisqdiff', 'p_nested')
  pops1 = setdiff(names(out1$popdrop), nams)
  pops2 = setdiff(names(out2$popdrop), nams)
  stopifnot(all(sort(pops1) == sort(pops2)))

  p2 = out1$popdrop %>% select(pat, dof, chisq, pops1) %>%
    type_convert(col_types = cols()) %>% mutate(lab = 'x') %>%
    bind_rows(out2$popdrop %>% select(pat, dof, chisq, pops1) %>%
                type_convert(col_types = cols()) %>% mutate(lab = 'y')) %>%
    pivot_longer(pops1, 'pop', values_to = 'weight') %>%
    pivot_longer(c(chisq, weight), 'key', values_to = 'value') %>%
    pivot_wider(names_from = lab, values_from = value) %>%
    filter(!is.na(x), !is.na(y)) %>%
    ggplot(aes(x, y, col = pat)) + geom_point() + facet_wrap(~ key, scales = 'free') +
    geom_abline() +
    theme(panel.background = element_blank(), axis.line = element_line(), legend.position = 'none') +
    xlab(paste0(name1)) +
    ylab(paste0(name2))

  gridExtra::grid.arrange(p1, p2, layout_matrix = matrix(c(1, 1, 2, 2), 2))
}

#' Plot an admixture graph using plotly
#' @export
#' @param graph An admixture graph
#' @param collapse_threshold Collapse nodes if they are separated by less than this distance (for fitted graphs)
#' @param fix If \code{TRUE}, there will be an attempt to rearrange the nodes to minimize
#' the number of intersecting edges. This can take very long for large graphs.
#' By default this is only done for graphs with fewer than 10 leaves.
#' @param highlight_unidentifiable Highlight unidentifiable edges in red. Can be slow for large graphs. See \code{\link{unidentifiable_edges}}.
#' @param pos Optional data frame with node coordinates (columns `node`, `x`, `y`)
#' @return A plotly object
#' @examples
#' plotly_graph(example_graph)
plotly_graph = function(graph, collapse_threshold = 0, fix = FALSE,
                        print_highlow = FALSE, highlight_unidentifiable = FALSE, pos = NULL,
                        nudge_y = -0.1, annot = '') {


  if(class(graph)[1] == 'igraph') {
    edges = as_edgelist(graph) %>% as_tibble(.name_repair = ~c('from', 'to'))
  } else {
    if(is.null(colnames(graph))) colnames(graph) = paste0('V', seq_len(ncol(graph)))
    edges = graph %>% as_tibble()
    names(edges)[1:2] = c('from', 'to')
    graph = igraph::graph_from_edgelist(as.matrix(graph)[,1:2])
  }

  if(collapse_threshold > 0) {
      edges %<>% collapse_edges(collapse_threshold)
  }

  admixnodes = unique(edges[[2]][edges[[2]] %in% names(which(table(edges[[2]]) > 1))])
  graph = edges %>% select(1:2) %>% as.matrix %>% igraph::graph_from_edgelist()

  if(is.null(pos)) {
    #pos = dot_coords(graph)
    coordmat = igraph::layout_as_tree(graph)
    pos = data.frame(names(V(graph)), coordmat, stringsAsFactors = F) %>%
      set_colnames(c('node', 'x', 'y'))
  }

  eg = edges %>% left_join(pos, by=c('from'='node')) %>%
    left_join(pos %>% transmute(to=node, xend=x, yend=y), by='to') %>%
    mutate(type = ifelse(to %in% admixnodes, 'admix', 'normal'))

  if(fix) {
    eg = fix_layout(eg %>% rename(name = from), graph) %>% rename(from = name)
  }
  if(TRUE) {
      eg = fix_shiftdown(eg %>% rename(name = from), graph) %>% rename(from = name)
  }

  if('weight' %in% names(edges)) {
    edgemul = 1000
    fa = function(x) paste0(round(x*100), '%')
    fe = function(x) round(x*edgemul)
    edges$label = ifelse(edges$type == 'admix', fa(edges$weight), fe(edges$weight))
    if(print_highlow) {
      edges %<>% mutate(label = ifelse(is.na(low), label, ifelse(type == 'admix',
                                paste(label, paste0('[', fa(low), '-', fa(high), ']'), sep = '\n'),
                                paste(label, paste0('[', fe(low), '-', fe(high), ']'), sep = '\n'))))
    }
    # if(!is.null(qpgraph_ranges)) {
    #   edges %<>% left_join(qpgraph_ranges, by = c('from', 'to')) %>%
    #     mutate(label = ifelse(type == 'admix',
    #                           paste(fa(mid), paste0('[', fa(lo), '-', fa(hi), ']'), sep = '\n'),
    #                           paste(fe(mid), paste0('[', fe(lo), '-', fe(hi), ']'), sep = '\n')))
    # }
  }
  if(!'label' %in% names(edges)) edges %<>% mutate(label='')
  e2 = edges %>% transmute(from, to, label) %>% left_join(count(., to), by = c('from'='to'))
  eg %<>% left_join(e2, by=c('from', 'to')) %>% mutate(indegree = replace_na(n, 0))
  nodes = eg %>% filter(to %in% get_leafnames(graph)) %>% rename(x=xend, y=yend, xend=x, yend=y) %>%
    transmute(name = to, x, y, xend, yend, to=NA, rownum = 1:n())

  nadmix = length(admixnodes)
  textsize = 2.5

  allnodes = eg %>% transmute(from, x, y, xend=0, yend=0, to=0) %>% filter(!duplicated(.))
  segtext = "ifelse(indegree == 1, to, from)"
  segtext = "paste(from, to, sep = ' -> ')"

  suppressWarnings({
    gg = eg %>% mutate(rownum = 1:n()) %>%
    ggplot(aes(x = x, xend = xend, y = y, yend = yend, from = from, to = to)) +
    geom_segment(aes_string(linetype = 'type', col = 'as.factor(y)', text = segtext),
                 arrow=arrow(type = 'closed', angle = 10, length = unit(0.15, 'inches'))) +
    geom_point(data = allnodes, aes(x, y, text = from), col = 'black', alpha = 0) +
    theme(panel.background = element_blank(),
          axis.line = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.position = 'none') +
    xlab('') + ylab('') +
    scale_linetype_manual(values = c(admix=3, normal=1)) +
    ggtitle('') +
    scale_x_continuous(expand = c(0.1, 0.1))

    if(highlight_unidentifiable) {
      unid = unidentifiable_edges(graph)
      unid2 = eg %>% right_join(unid %>% select(-type), by = c('from', 'to'))
      gg = gg + geom_segment(aes_string(linetype = 'type'), col = 'red', data = unid2, size = 1)
    }
    gg = gg + geom_text(aes(x = (x+xend)/2, y = (y+yend)/2, label = label, text = paste(from, to, sep = ' -> ')),
                            size = textsize) +
      geom_text(data = nodes, aes_string(label = 'name', col = 'as.factor(yend)', from = NA),
                        size = textsize, nudge_y = nudge_y) + expand_limits(y = c(0,1)) +
      annotate('text', min(eg$x), max(eg$y), label = annot)
  })
  plt = plotly::ggplotly(gg, tooltip=c('text'))
  plt
}



dot_coords = function(graph) {

  z = Rgraphviz::agopen(igraph::igraph.to.graphNEL(graph), 'bla', layoutType = 'dot')
  coordmat = map(seq_along(z@AgNode), ~{coords = z@AgNode[[.]]@center; c(coords@x, coords@y)}) %>%
    do.call(rbind, .) %>% magrittr::divide_by(100)
  data.frame(names(V(graph)), coordmat, stringsAsFactors = F) %>%
    set_colnames(c('node', 'x', 'y'))
}


static_coords = function(graph, pops = NULL) {

  if(is.null(pops)) pops = get_leafnames(graph)
  pcs = map(1:length(pops), ~apply(combn(pops, .), 2, function(x) paste(x, collapse = ' ')))
  internal = setdiff(names(V(graph)), pops)

  coord = graph %>% distances(names(V(graph)), pops, mode = 'out') %>% as_tibble(rownames = 'from') %>%
    pivot_longer(-from, names_to = 'to', values_to = 'order') %>% filter(is.finite(order)) %>%
    select(-order) %>% group_by(from) %>% summarize(id = paste(to, collapse = ' '), n = n()) %>% rowwise %>% mutate(xx = match(id, pcs[[n]])/length(pcs[[n]])-1/(length(pcs[[n]])+1))

  enframe(names(V(graph))) %>% left_join(coord, by = c('value' = 'from')) %>%
    mutate(n = (length(pops)-n)/length(pops), xx=xx*2*pi) %>% transmute(node = value, x = n, y = xx) #%>%
    #mutate(xx = x*cos(y), yy = x*sin(y)) %>% transmute(node, x = xx, y = yy)
}


