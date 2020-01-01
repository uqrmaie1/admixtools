#' Plot a comparison of two lists of qpgraph output
#' @export
#' @param out1 first qpgraph output.
#' @param out2 second qpgraph output.
#' @return a ggplot object.
#' @examples
#' fit1 = qpgraph(graph1, f2_blocks, block_lengths, cpp = TRUE)
#' fit2 = qpgraph(graph1, f2_blocks, block_lengths, cpp = FALSE)
#' plot_comparison(fit1, fit2)
plot_comparison = function(out1, out2) {
  # plots a comparison of two qpgraph outputs

  if(max(nchar(out1$f2$pop1))==3 || max(nchar(out2$f2$pop1))==3) {
    stopifnot(length(unique(c(out1$f2$pop1, out1$f2$pop2))) ==
                length(unique(substr(c(out1$f2$pop1, out1$f2$pop2), 1, 3))) &&
                length(unique(c(out2$f2$pop1, out2$f2$pop2))) ==
                length(unique(substr(c(out2$f2$pop1, out2$f2$pop2), 1, 3))))
    out1$f2$pop1 = substr(out1$f2$pop1, 1, 3)
    out1$f2$pop2 = substr(out1$f2$pop2, 1, 3)
    out2$f2$pop1 = substr(out2$f2$pop1, 1, 3)
    out2$f2$pop2 = substr(out2$f2$pop2, 1, 3)
  }

  f2comp = out1$f2 %>% select('pop1', 'pop2', 'f2est', 'se') %>% mutate(i = 'x') %>%
    bind_rows(out2$f2 %>% select('pop1', 'pop2', 'f2est', 'se') %>% mutate(i = 'y')) %>%
    gather('type', 'v', .data$f2est, .data$se) %>% spread(.data$i, .data$v) %>%
    rename(from=.data$pop1, to=.data$pop2)

  out1$edges %>% rename(x = .data$weight) %>%
    left_join(out2$edges %>% select(-'type') %>% rename(y = .data$weight), by=c('from', 'to')) %>%
    bind_rows(f2comp) %>%
    ggplot(aes(.data$x, .data$y)) + geom_point() + facet_wrap('type', scales='free') + geom_abline() +
    xlab(paste0('stats 1 (score: ', round(out1$score, 2),')')) +
    ylab(paste0('stats 2 (score: ', round(out2$score, 2),')')) +
    theme(panel.background = element_blank(), axis.line = element_line()) +
    geom_smooth(method='lm', se=F, formula=y~x-1)

}


# Plot an admixture graph
# @param edges a two column matrix specifying the admixture graph. first column is source node, second column is target node. the first edge has to be root -> outgroup. admixture nodes are inferred as those nodes which are the targets of two different sources.
# @param layout argument passed to \code{\link[ggdag]{tidy_dagitty}} to modify plot layout.
# @return a ggplot object.
plot_graph_old = function(edges, layout='tree', fix=TRUE, title = '') {

  if (!requireNamespace('ggdag', quietly = TRUE)) {
    stop('Package "ggdag" needed for this function to work. Please install it.', call. = FALSE)
  }
  if(class(edges)[1] == 'igraph') edges = as_edgelist(edges)
  edges = as_tibble(edges, .name_repair='unique')
  names(edges)[1:2] = c('V1', 'V2')
  x = ggdag::dag(paste(edges[[1]], edges[[2]], sep='->', collapse=' '))
  admixnodes = unique(edges[[2]][edges[[2]] %in% names(which(table(edges[[2]]) > 1))])

  dat = ggdag::tidy_dagitty(x, layout=layout)
  if(layout == 'tree' && fix) dat$data = fix_layout(dat$data, igraph::graph_from_edgelist(as.matrix(edges)[,1:2]))
  dat %<>% group_by(.data$to) %>% mutate(cnt = n()) %>% ungroup %>%
    mutate(admix = (.data$cnt > 1)*2+1) %>% ungroup %>%
    #select(-.data$type) %>%
    mutate(type = ifelse(.data$name == 'R', 'root',
                         ifelse(!.data$name %in% edges[[1]], 'leaf',
                                ifelse(.data$name %in% admixnodes, 'admix', 'normal'))))
  if(!'label' %in% names(edges)) edges %<>% mutate(label='')
  dat$data %<>% left_join(y=edges, by=c('name'='V1','to'='V2'))
  plt = dat %>% ggplot(aes(x = .data$x, y = .data$y, xend = .data$xend, yend = .data$yend)) +
    ggdag::geom_dag_edges(aes(edge_linetype=.data$admix, label=.data$label), edge_colour='grey') +
    ggdag::geom_dag_label_repel(aes(fill=.data$type, label=.data$name), force=0, box.padding=0) +
    ggdag::theme_dag(legend.position='none') + ggtitle(label = title) +
    scale_fill_manual(values=c('admix'=gg_color_hue(3)[1],
                               'leaf'=gg_color_hue(3)[2],
                               'normal'=gg_color_hue(3)[3],
                               'root'='#FFFFFF'))
  plt
}

#' Plot an admixture graph
#' @export
#' @param edges an admixture graph. If it's an edge list matrix or data frame with a \code{label} column, those values will displayed on the edges
#' @param title plot title
#' @return a ggplot object
#' @examples
#' plot_graph(graph1)
plot_graph = function(edges, fix = TRUE, title = '', color = TRUE) {

  if(class(edges)[1] == 'igraph') {
    grph = edges
    edges = as_edgelist(edges)
  } else {
    grph = igraph::graph_from_edgelist(as.matrix(edges)[,1:2])
  }
  edges = purrr::quietly(as_tibble)(edges, .name_repair='unique')[[1]]
  names(edges)[1:2] = c('V1', 'V2')
  admixnodes = unique(edges[[2]][edges[[2]] %in% names(which(table(edges[[2]]) > 1))])

  pos = data.frame(names(V(grph)), igraph::layout_as_tree(grph), stringsAsFactors = F) %>% set_colnames(c('node', 'x', 'y'))
  eg = as_tibble(as_edgelist(grph)) %>% left_join(pos, by=c('V1'='node')) %>% left_join(pos %>% transmute(V2=node, xend=x, yend=y), by='V2') %>% mutate(type = ifelse(V2 %in% admixnodes, 'admix', 'normal')) %>% rename(name=V1, to=V2)

  if(fix) eg = fix_layout(eg, grph)
  # need to update fix to work here; fix assumes that eg has node rows in the end; only relevant in fix_layout2?

  if(!'label' %in% names(edges)) edges %<>% mutate(label='')
  eg %<>% left_join(edges %>% transmute(name=V1, to=V2, label), by=c('name', 'to'))
  nodes = eg %>% filter(to %in% get_leafnames(grph)) %>% rename(x=xend, y=yend, xend=x, yend=y) %>% transmute(name = to, x, y, xend, yend)

  if(color) {
    gs = geom_segment(aes_string(linetype = 'type', col = 'as.factor(y)'),
                      arrow=arrow(type = 'closed', angle = 10, length=unit(0.15, 'inches')))
    gl = geom_label(data=nodes, aes_string(label = 'name', col='as.factor(yend)'), size=3)
  } else {
    gs = geom_segment(aes_string(linetype = 'type'), col = 'grey',
                      arrow=arrow(type = 'closed', angle = 10, length=unit(0.15, 'inches')))
    gl = geom_label(data=nodes, aes_string(label = 'name'), col = 'black', size=3)
  }

  plt = eg %>%
    ggplot(aes(x=x, xend=xend, y=y, yend=yend)) +
    gs +
    geom_text(aes(x = (x+xend)/2, y = (y+yend)/2, label = label)) +
    gl +
    theme(panel.background = element_blank(),
          axis.line = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.position = 'none') +
    xlab('') + ylab('') +
    scale_linetype_manual(values = c(3, 1)) +
    ggtitle(title) +
    scale_x_continuous(expand = c(0.15, 0.15))

  plt

}



fix_layout = function(coord, grph) {
  # rearranges nodes in tree layout, so that there are fewer long branches
  dst = igraph::distances(grph, V(grph)[1], mode='out')[1,]
  totdist = function(x, xend) {
    sum(abs(x-xend), na.rm=TRUE)
  }
  swap_subtree = function(grph, coord, node) {
    center = coord$x[match(node, coord$name)[1]]
    offspring = names(subcomponent(grph, node, mode='out'))
    dst2 = igraph::distances(grph, node, mode='out')[1,]
    for(n in offspring) {
      if(dst[node] + dst2[n] <= dst[n])
        coord %<>% mutate(x = ifelse(name == n, 2*center-x, x),
                          xend = ifelse(to == n, 2*center-xend, xend))
    }
    coord
  }
  swap_subtree2 = function(grph, coord, node) {
    center = coord$x[match(node, coord$name)[1]]
    children = names(neighbors(grph, node, mode='out'))
    dst2 = igraph::distances(grph, node, mode='out')[1,]
    if(length(children) == 2) {
      for(child in children) {
        offspring = names(subcomponent(grph, child, mode='out'))
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
  inner_bf = setdiff(names(subcomponent(grph, V(grph)[1])), get_leafnames(grph))
  for(node in rep(inner_bf, 1)) {

    d1 = totdist(coord$x, coord$xend)
    n1 = num_intersecting(coord)
    s1 = num_samepos(coord)
    c2 = swap_subtree(grph, coord, node)
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
    c2 = swap_subtree2(grph, coord, node)
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
    nodeshift %<>% group_by(g = interaction(y, yend, cumsum(cummax(lag(xend, default = first(xend))) < x))) %>% mutate(tot = n(), this = (seq_len(tot[1]))-1, shift = this/tot) %>% filter(shift != 0) %>% dplyr::select(g, name, to, shift) %>% ungroup %>% gather(irl, node, name, to) %>% dplyr::select(node, shift)
    if(nrow(nodeshift) > 0) {
      edges %<>% left_join(nodeshift %>% transmute(name=node, s1 = shift), by='name') %>% left_join(nodeshift %>% transmute(to=node, s2 = shift), by='to') %>% mutate(y = ifelse(is.na(s1), y, y+s1), yend = ifelse(is.na(s2), yend, yend+s2))
      leaves %<>% left_join(nodeshift %>% transmute(name=node, shift), by='name') %>% mutate(y=ifelse(is.na(shift), y, y+shift))
    }
  }
  bind_rows(edges, leaves)
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


#' @export
collapse_edges = function(edges, below=0) {
  # input: graph edgelist
  # requires column 'weight'
  # collapses nodes which are separated by 'below' or less
  names(edges)[1:2] = c('from', 'to')
  admixnodenames = edges %>% filter(type == 'admix') %$% unique(to)
  leaves = setdiff(edges[[2]], edges[[1]])
  excl = c(admixnodenames, leaves)
  adjmat = edges %>% select(1, 2, weight) %>% bind_rows(rename(., from=to, to=from)) %>%  spread(to, weight, fill = 0) %>% column_to_rownames(var = 'from') %>% as.matrix
  adjmat[adjmat > below] = 0
  adjmat[rownames(adjmat) %in% excl, ] = 0
  adjmat[, colnames(adjmat) %in% excl] = 0
  diag(adjmat) = 1
  comp = components(graph_from_adjacency_matrix(as.matrix(adjmat), weighted='x'))
  newnam = comp$membership %>% enframe %>% group_by(value) %>% mutate(name, groupname = paste(name, collapse='.')) %>% ungroup %>% select(-value) %>% deframe
  edges %>% mutate(from = newnam[from], to = newnam[to]) %>% filter(from != to)

}



#' @export
plot_graph_pcs = function(grph, pcs) {

  leafcoords = pcs %>% rename(lon=PC1, lat=PC2)
  leafcoords2 = leafcoords %>% group_by(group) %>% summarize(x = mean(lon), y = mean(lat))
  sg = grph %>% simplify_graph
  node_coord = node_coords_3d(sg, leafcoords2)
  edge_coord = sg %>% as_edgelist %>% as.data.frame(stringsAsFactors=F) %>% set_colnames(c('from', 'to')) %>% mutate(eid = 1:n()) %>% gather(type, node, -eid) %>% left_join(node_coord, by='node') %>% arrange(eid)
  edge_coord %<>% mutate(admix = eid %in% (filter(., type == 'to') %>% group_by(node) %>% mutate(cnt = n()) %>% filter(cnt == 2) %$% eid), admix=ifelse(admix, 'dot', 'solid'))
  node_coord2 = leafcoords %>% rename(x=lon, y=lat) %>% bind_rows(node_coord %>% filter(z==0) %>% transmute(s = 'center', group=node, x, y)) %>% mutate(z=0, sz = ((s == 'center')*5+2))

  ax = list(visible=FALSE)

  plot_ly(type = 'scatter3d', mode = 'lines', hoverinfo = 'text') %>%
    add_trace(data = edge_coord %>% mutate(const='grey'), x = ~x, y = ~y, z = ~z,
              split = ~eid, line=list(color=~const, dash=~admix), showlegend=F) %>%
    add_markers(data = node_coord2, x = ~x, y = ~y, z = ~z,
                hovertext = ~group, marker=list(size=~((s=='center')*5+5)), color = ~group, colors='Set1') %>%
    graphics::layout(scene = list(xaxis=list(title='PC1'), yaxis=list(title='PC2'), zaxis=ax,
                                  camera = list(projection = list(type = 'orthographic'),
                                                eye = list(x = 0, y = 0, z = 1),
                                                up = list(x = 0, y = 1, z = 0))))

}




#' Plot an admixture graph on a map
#' @export
#' @param grph a two column matrix specifying the admixture graph. first column is source node, second column is target node. the first edge has to be root -> outgroup. admixture nodes are inferred as those nodes which are the targets of two different sources.
#' @param leafcoords data frame with columns \code{group}, \code{lon}, \code{lat}
#' @return a plotly object.
#' @examples
# #' \dontrun{
#' plot_graph_map(igraph1, anno)
# #' }
plot_graph_map = function(grph, leafcoords, shapedata=NULL) {

  stopifnot(all(c('iid', 'group', 'lat', 'lon') %in% names(leafcoords)))
  leafcoords %<>% filter(group %in% get_leafnames(grph))
  leafcoords2 = leafcoords %>% group_by(group) %>% summarize(x = mean(lon), y = mean(lat))
  sg = grph %>% simplify_graph
  node_coord = node_coords_3d(sg, leafcoords2)
  edge_coord = sg %>% as_edgelist %>% as.data.frame(stringsAsFactors=F) %>% set_colnames(c('from', 'to')) %>% mutate(eid = 1:n()) %>% gather(type, node, -eid) %>% left_join(node_coord, by='node') %>% arrange(eid)
  edge_coord %<>% mutate(admix = eid %in% (filter(., type == 'to') %>% group_by(node) %>% mutate(cnt = n()) %>% filter(cnt == 2) %$% eid), admix=ifelse(admix, 'dot', 'solid'))
  node_coord2 = leafcoords %>% rename(x=lon, y=lat) %>% bind_rows(node_coord %>% filter(z==0) %>% transmute(iid = 'center', group=node, x, y)) %>% mutate(z=0, sz = ((iid == 'center')*5+2))

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
# @param grph an admixture graph as an igraph object.
# @param leafcoords data frame with columns \code{group}, \code{lon}, \code{lat}
# @param map_layout 1 or 2
# @return a ggplot object.
# @examples
# \dontrun{
# plot_graph_map2(igraph1, anno, 1)
# plot_graph_map2(igraph1, anno, 2)
# }
plot_graph_map2 = function(grph, leafcoords, map_layout = 1) {

  stopifnot(all(c('iid', 'group', 'lat', 'lon') %in% names(leafcoords)))
  leafcoords %<>% filter(group %in% get_leafnames(grph))
  leafcoords2 = leafcoords %>% group_by(group) %>% summarize(x = mean(lon), y = mean(lat))
  sg = grph %>% simplify_graph
  node_coord = node_coords_3d(sg, leafcoords2)
  edge_coord = sg %>% as_edgelist %>% as.data.frame(stringsAsFactors=F) %>% set_colnames(c('from', 'to')) %>% mutate(eid = 1:n()) %>% gather(type, node, -eid) %>% left_join(node_coord, by='node') %>% arrange(eid)
  node_coord2 = leafcoords %>% rename(x=lon, y=lat) %>% bind_rows(node_coord %>% filter(z==0) %>% transmute(iid = 'center', group=node, x, y)) %>% mutate(z=0, sz = ((iid == 'center')*5+2))

  ax = list(visible=FALSE)

  src1 = 'https://basemap.nationalmap.gov/arcgis/rest/services/USGSImageryOnly/MapServer/tile/{z}/{y}/{x}'
  if(map_layout == 1) lo = function(x) plotly::layout(x, mapbox = list(style = 'stamen-terrain'), xaxis=ax, yaxis=ax)
  if(map_layout == 2) lo = function(x) plotly::layout(x, mapbox = list(style = 'white-bg', zoom = 1, layers = list(list(below = 'traces', sourcetype = 'raster', source = list(src1)))))

  plotly::plot_ly(type = 'scattermapbox', mode='lines', hoverinfo = 'text') %>%
    plotly::add_trace(data = edge_coord %>% mutate(const='grey', lon=x, lat=y), lon = ~x, lat = ~y,
              split = ~eid, line=list(color=~const, width=1), showlegend=F) %>%
    plotly::add_trace(data = node_coord2, mode='markers', x = ~x, y = ~y,
              hovertext = ~group, marker=list(size=~((iid=='center')*5+5)), color = ~group, colors='Set1') %>%
    lo

}

#' Plot samples on a map
#' @export
#' @param leafcoords data frame with columns \code{group}, \code{lon}, \code{lat}
#' @param map_layout 1 or 2
#' @return a plotly object.
#' @examples
#' \dontrun{
#' plot_map(anno, 1)
#' p1 = plot_map(anno, 1)
#' }
plot_map = function(leafcoords, map_layout = 1) {

  stopifnot(all(c('iid', 'lat', 'lon', 'yearsbp') %in% names(leafcoords)))
  #leafcoords2 = leafcoords %>% group_by(group) %>% summarize(lon = mean(lon), lat = mean(lat))

  ax = list(visible=FALSE)

  src1 = 'https://basemap.nationalmap.gov/arcgis/rest/services/USGSImageryOnly/MapServer/tile/{z}/{y}/{x}'
  if(map_layout == 1) lo = function(x) plotly::layout(x, mapbox = list(style = 'stamen-terrain', center = list(lat = 48.2, lon = 16.3)), xaxis=ax, yaxis=ax)
  if(map_layout == 2) lo = function(x) plotly::layout(x, mapbox = list(style = 'white-bg', center = list(lat = 48.2, lon = 16.3),  zoom = 1, layers = list(list(below = 'traces', sourcetype = 'raster', source = list(src1)))))

  plotly::plot_ly(leafcoords %>% filter(between(lat, -90, 90), between(lon, -180, 180), yearsbp > 0), type = 'scattermapbox', mode='markers', hoverinfo = 'text', x = ~lon, y = ~lat, color = ~log10(yearsbp), hovertext = ~group) %>% lo

}



#' @export
node_coords_3d = function(grph, leafcoords, rootlon=NA, rootlat=NA) {
  # given a graph and 2d coordinates of leaves,
  # this function returns 3d coords for all nodes
  root = names(V(grph)[1])
  leaves = get_leafnames(grph)
  #leafcoords %<>% filter(!is.na(x), !is.na(y))
  leafcoords %<>% mutate(x = ifelse(is.na(x), runif(1, min(x,na.rm=T), max(x,na.rm=T)), x),
                        y = ifelse(is.na(y), runif(1, min(y,na.rm=T), max(y,na.rm=T)), y))
  stopifnot(all(leaves %in% leafcoords$group))
  paths = grph %>% igraph::all_simple_paths(root, leaves, mode='out')
  out = leafcoords %>% mutate(z = 0) %>%
    bind_rows(tibble(group = root, x = mean(.$x, na.rm=T), y = mean(.$y, na.rm=T), z = 1))
  rootcoord = out %>% filter(group == root) %$% c(x=x, y=y, z=z)
  if(!is.na(rootlon) & !is.na(rootlat)) {
    rootcoord = c(x=rootlon, y=rootlat, z=1)
  }

  path_coords = function(path, startcoord, endcoord) {
    # given a sequence of nodes where the first and last node have 3d coords, return xyz for all intermediate nodes
    map2(startcoord, endcoord, ~seq(.x, .y, length.out=length(path))) %>% bind_cols %>% add_column(node=path, .before=1)
  }

  paths %>% map(names) %>% map(~path_coords(., rootcoord, filter(out, group == tail(., 1)) %$% c(x=x, y=y, z=z))) %>%
    bind_rows(.id='pid') %>% group_by(node) %>% summarize(x = mean(x), y = mean(y), z = mean(z))

}




