
# nodeHeights gives: (1) the time of the ancestor node of the edge, (2) the time of the descendant node of the edge

# GET DESCENDANTS (my way)
# once is a logical that asks if you want the descendants only one depth upwards (not all descendant nodes, but the two forward ones)
descendants <- function(tree, node, once = F) {
  desc <- c()  # to store descendants
  while (T) {
    node <- tree$edge[tree$edge[,1] %in% node,2]
    desc <- c(desc, node)
    
    if (once) {
      break
    }
    
    if (length(node) == 0) {
      break
    }
  }
  
  return(sort(desc))
}

# PLOTTING INTERNAL NODE NUMBERS
plot.internal <- function (phy, add = F, ...) {
  phy$node.label <- (length(phy$tip.label) + 1):max(phy$edge) # the numbers for all internal nodes, added as an attribute of the phy object
  
  if (add) {
    par(new = T)
    plot.phylo(phy, show.node.label=T, edge.color=opacity("white", alpha = 0.0000001), ...)
  } else {
    plot.phylo(phy, show.node.label = T, ...)
  }
}


# ARC PHYLOGENY

phylo_arc <- function(tree,
                      max_th = 180,
                      inner_radius = 0.5,
                      scale_tree = 1, # number between 1 and 0 that reduces the size of the tree
                      axis_lim = 10,  # a single number, what you want the maximum x and y value to be in the tree
                      edge.lwd = 1,
                      edge.col = "black",
                      plot = T) {
  
  descendants <- function(tree, node, once = F) {
    desc <- c()  # to store descendants
    while (T) {
      node <- tree$edge[tree$edge[,1] %in% node,2]
      desc <- c(desc, node)
      
      if (once) {
        break
      }
      
      if (length(node) == 0) {
        break
      }
    }
    
    return(sort(desc))
  }
  
  ntip = length(tree$tip.label)
  ntot <- ntip + tree$Nnode
  radius = axis_lim * scale_tree
  inner <- inner_radius * radius
  
  y_vals = node.height(tree) - 1  # gets y values (theta)  # need to subtract one cause the first node's y should be y = 0 not 1
  x_vals = branching.times(tree) #  gets x values, backwards (1 is at the left) (r)
  max_x <- max(x_vals)
  max_y <- max(y_vals)
  
  y_vals <- y_vals / max_y   # scale y values to the 0 angle
  x_vals <- 1 - (x_vals / max_x)  # scale x values to a radius of 0
  
  th <- y_vals * max_th 
  r <- (x_vals * (radius - inner)) + inner
  r <- c(rep(radius, ntip), r)  # add in the tip nodes
  
  length(r) == ntot
  length(th) == ntot
  
  # now we convert to cartesian and plot
  node_pos <- toCart(r, th, deg = T)
  
  if (plot) {
    # plotting
    plot.new()
    pin <- par("pin")  # taking current plot dimensions
  
    if (max_th <= 180) {
      xlim <- c(-axis_lim, axis_lim)
      ylim <- c(0, axis_lim)
    } else {
      xlim <- c(-axis_lim, axis_lim)
      ylim <- c(toCart(axis_lim, ifelse(max_th < 270, max_th, 270))$y, axis_lim)
    }
  
    if (pin[1] > pin[2]) {
      xlim <- (pin[1]/pin[2]) * xlim
    } else {
      ylim <- (pin[2]/pin[1]) * ylim
    }
    plot.window(xlim, ylim, "", asp = 1)
    #points(node_pos$x, node_pos$y, pch = 16)
    
    # connect the dots
    for (i in ntot:(ntip + 1)) {  # for each internal node, backwards
      desc <- descendants(tree, i, once = T)
      
      th_desc <- th[desc]
      r_desc <- r[desc]
      
      # draw the vertical lines
      for (j in 1:length(desc)) {
        endpoints <- toCart(r = c(r_desc[j], r[i]), th = th_desc[j])
        segments(x0 = endpoints$x[1], y0 = endpoints$y[1], x1 = endpoints$x[2], y1 = endpoints$y[2],
                 lwd = edge.lwd, col = edge.col)
      }
      
      # draw the horizontal lines
      edges <- floor((max(th_desc) - min(th_desc)) * 10)
      endpoints <- toCart(r = r[i], th = seq.int(max(th_desc), min(th_desc), length.out = edges))
      for (j in 1:(edges-1)) {
        segments(x0 = endpoints$x[j], y0 = endpoints$y[j], x1 = endpoints$x[j+1], endpoints$y[j+1],
                 lwd = edge.lwd, col = edge.col)
      }
    }
  }
}

tree <- chronopl(rtree(20), lambda=1)

phylo_arc(tree, edge.lwd = 2, inner_radius = 0.5, scale_tree = 1, max_th = 350)
plotTree(tree, ftype = "off", type = "arc")

# COLLAPSE GENUS

# function will collapse genus names to as few letters as possible to make each collapsed genus unique
# function keeps entire species name intact
# x must be a vector containing genus and species names, connected by any separator
# max_genus_letters constrains the maximum number of letters in the genus name (may not lead to unique genus names)
# min_genus_letters does the same but constrains the minimum number of letters (default minimum is 1)
# sep sets the separating character for the final collapsed genus and species
# genus_punct allows addition of characters after the collapsed genus (e.g. "T." instead of just "T")

collapseGenus <- function(x, max_genus_letters = NULL, min_genus_letters = NULL, sep = " ", genus_punct = "") {
  
  input_sep <- unique(gsub("[[:alnum:]]+(\\s+|[[:punct:]]+)\\w+.*", "\\1", x))  # i want to select the thing that separates the genus and species
  
  genera <- unique(unlist(lapply(strsplit(x, split = input_sep), "[", 1)))
  genera_spl <- strsplit(genera, split = "")
  
  collapse_to_which <- rep(1, length(genera)) |> setNames(nm = genera)
  i = 1
  letters_i <- unlist(lapply(genera_spl, "[", i))  # extract all letters at place i in the names
  
  while (sum(duplicated(letters_i)) > 0) {
    index <- duplicated(letters_i) | duplicated(letters_i, fromLast = T)
    index2 <- collapse_to_which == i
    index3 <- ((1:length(collapse_to_which))[index2])[index]
    collapse_to_which[index3] <- collapse_to_which[index3] + 1
    i = i + 1
    letters_i <- unlist(lapply(genera_spl, "[", i))[index3]
  }
  
  if (!is.null(max_genus_letters)) {   # if want a maximum number of letters that genus can be collapsed to
    collapse_to_which[collapse_to_which > max_genus_letters] <- max_genus_letters
  }
  if (!is.null(min_genus_letters)) {
    collapse_to_which[collapse_to_which < min_genus_letters] <- min_genus_letters
  }
  
  genera_spl_len <- unlist(lapply(genera_spl, length))
  collapsed_nm <- character(length(x))
  for (i in seq_along(x)) {
    name <- unlist(strsplit(x[i], split = input_sep))
    genus <- name[1]
    species <- name[2]
    col_to <- collapse_to_which[names(collapse_to_which) == genus]
    
    flag = ifelse(length(unlist(strsplit(genus, split = ""))) <= col_to, T, F)
    collapsed_nm[i] <- paste(gsub(paste0("([[:alpha:]]{", col_to, "}).*"), 
                                  ifelse(as.logical(flag), "\\1", paste0("\\1", genus_punct)), 
                                  genus), 
                             species, sep = sep)
  }
  
  return(collapsed_nm)
}

x <- c("Leptophis_ahaetulla", "Oxybelis_aenus", "Chironius_fuscus", "Leptophis_mexicanus", 
       "Leptodeira_annulata", "Crotalus_durissus", "Boa_constrictor", "Tantilla_he", "Tantilla_re")
collapseGenus(x, genus_punct = ".")
collapseGenus(x, max_genus_letters = 4, genus_punct = ".")
collapseGenus(x, min_genus_letters = 2, sep = "-")

# PLOTTING RECTANGLES ON BRANCHES IN A TREE

# into EDGE, pass a vector with the nodes in the right order, ancestor node first then descendant node
# into NREC, pass integer of number of rectangles
# into HEIGHT, pass number with height of rectangle that you want (actually half the total height you want)
# into COLS, pass vector of colors same length as NREC or else it'll recycle
# into BORDER, pass the color for the borders of rectangles
# need to have plotted the right phylogeny before running this, or else it won't work

rect_phy <- function(tree, edge, nrec, height, width, cols, border) {
  edges <- cbind(tree$edge, tree$edge.length, nodeHeights(tree))
  index <- which(edges[,1] == edge[1] & edges[,2] == edge[2])
  target <- edges[index,]
  
  slice <- rep(target[3] / (nrec + 1), nrec) |> cumsum()
  stakes <- target[4] + slice
  
  obj <- get("last_plot.phylo", envir=.PlotPhyloEnv)
  mid = obj$yy[edge[2]]
  
  for (i in 1:nrec) {
    polygon(x = c(stakes[i] - width, stakes[i] + width, stakes[i] + width, stakes[i] - width), 
            y = c(mid - height, mid - height, mid + height, mid + height),
            col = cols[i], border = border)
  }
}

# PLOTTING A SQUIGGLY TO SHOW PASSAGE OF TIME ON A BRANCH

# into EDGE, pass a vector with the nodes in the right order, ancestor node first then descendant node
# into WIDTH, pass a decimal from 0 to 1 that represents the width of the squiggle that you want, it is the percentage of the branch you want to modify that will be taken up by the squiggle
# into PEAKS, pass the number of peaks and troughs that you want your squiggle to have
# into OFF, pass a percetange (between 0 and 1) that represents where along the branch you want the squiggle's middle to be. the default is 50%, or the exact center of the branch.
# into COL, pass the color of the edge and its squiggle
# into BG, pass the color of the background of your plot, i erase the edge in a hacky way so you have to pass it the background color
# you may get a warning about 'items to replace is not a multiple of replacement length' if the number of peaks you passed is odd, this is fine
# setting peaks to 0 redraws a straight line in the defined area

squiggle_phy <- function(tree, edge, width, peaks, off = 0.5, col = "black", bg = "white") {
  require(phytools)
  require(ape)
  
  edges <- cbind(tree$edge, tree$edge.length, nodeHeights(tree))
  index <- which(edges[,1] == edge[1] & edges[,2] == edge[2])
  target <- edges[index,]
  
  obj <- get("last_plot.phylo", envir=.PlotPhyloEnv)
  y = obj$yy[edge[2]]
  x1 = obj$xx[edge[1]]
  x2 = obj$xx[edge[2]]
  off = target[3] * off
  mid = target[4] + off
  width = target[3] * width  # the width will be a percentage of the total branch
  x_start = mid - (width/2) # x of the very start on the left
  
  segments(x0 = x_start, y0 = y, x1 = mid + (width/2), y1 = y, col = bg) # hacky erase the existing branch
  
  slice <- rep(width/(peaks+1), peaks+1) |> cumsum()
  x_vals <- c(x_start, x_start + slice)
  y_vals <- rep(y, peaks+2)
  if (length(y_vals) > 2) {
    y_vals[2:(length(y_vals) - 1)] <- c(y + 0.5, y - 0.5) # this can give a warning if peak number is odd, using the vector recycling property
  }
  
  for (i in 1:(peaks+1)) {
    segments(x0 = x_vals[i], y0 = y_vals[i], x1 = x_vals[i+1], y1 = y_vals[i+1], col = col)
  }
}


plot.internal(tree, cex = 0.4)
# fun tests

set.seed(1)
tree <- chronopl(rtree(20), lambda = 1)
plot(tree)
nodelabels(frame = "none", cex = 0.6)
squiggle_phy(tree, edge = c(21,22), width = 0.25, peaks = 2)
squiggle_phy(tree, edge = c(28,29), width = 0.5, peaks = 6)
squiggle_phy(tree, c(32,34), width = 0.15, peaks = 7, off = 0.8)

rect_phy(tree, edge = c(23,24), nrec = 6, height = 1.1, width = 0.005, cols = rainbow(6), border = "black")
rect_phy(tree, edge = c(25,26), nrec = 3, height = 0.5, width = 0.01, cols = c("lightpink2", "darkorchid3", "tan"), border = "white")
         
         
# conversions between polar and cartesian coords
toPolar <- function(x, y) {
  r = sqrt(x^2 + y^2)
  th = atan2(y, x)
  return(list(r = r, th = th))
}

toCart <- function(r, th, deg = T) {
  if (deg) {
    th <- th * (pi/180)
  }
  
  x = r * cos(th)
  y = r * sin(th)
  return(list(x = x, y = y))
}
         
         