# library(data.tree)
# library(rstack)
# library(pryr)

str_to_char <- function(string){
  chars <- strsplit(string, "")[[1]]
  remove <- which(chars == " ")
  if (length(remove) > 0){
    chars <- chars[-remove]
  }
  return(chars)
}

char_to_str <- function(chars){
  string <- paste(chars,sep = '', collapse = '')
  return(string)
}

extract_logical_from_base <- function(chars){
  and <- which(chars == "&")
  or <- which(chars == "|")

  if ((length(and) > 0) & (length(or) > 0)){
    print('not base')
    return()
  }

  if ((length(and) > 0)){
    return('&')
  }

  if ((length(or) > 0)){
    return('|')
  }
}

find_gpr_paths <- function(gprRule){

  gpr_char <- str_to_char(gprRule)
  len_logical <- which(gpr_char %in% c('&', '|'))
  # print(gpr_char)
  # print(len_logical)
  # print(extract_logical_from_base(gpr_char))

  if (length(len_logical) == 0){
    return(gsub("[()]", "", gprRule))
  }

  ast <- capture.output(call_tree(parse(text = gprRule)))

  logic_flag <- FALSE
  i <- 0 #starting idx
  while (!logic_flag){
    i <- i + 1
    if (('&' %in% str_to_char(ast[i])) | ('|' %in% str_to_char(ast[i]))){
      logic_flag <- TRUE
    }
  }

  # print(i)
  #
  # print(matrix(ast, nrow = length(ast), ncol = 1))
  gpr_tree <- tree_building(ast, i)
  # print(gpr_tree)
  paths <- get_paths_from_gpr(gpr_tree)

  return(paths)
}

tree_building <- function(list, index){

  chars <- strsplit(list[index], "")[[1]]
  og_depth <- (which(strsplit(list[index], "")[[1]] == "\\")-1)/2
  og_id <- chars[length(chars)]
  id1 <- "_" #og_id
  id2 <- "_"

  new_node <- Node$new(og_id)

  start_1 <- index + 1
  chars <- strsplit(list[start_1], "")[[1]]
  depth1 <- (which(chars == "\\")-1)/2
  id1 <- chars[length(chars)]

  start_2 <- start_1 + 1
  chars <- strsplit(list[start_2], "")[[1]]
  depth2 <- (which(chars == "\\")-1)/2
  id2 <- chars[length(chars)]
  while (depth2 != depth1 | id2 != id1){
    start_2 <- start_2 + 1
    chars <- strsplit(list[start_2], "")[[1]]
    depth2 <- (which(chars == "\\")-1)/2
    id2 <- chars[length(chars)]
  }


  # depth1 <- og_depth
  # child_ct <- 0
  index1 <- start_1
  # chars <- strsplit(list[index1], "")[[1]]
  # depth1 <- (which(chars == "\\")-1)/2
  # id1 <- chars[length(chars)]

  while (!(id1 %in% c( '|', '&', 'x'))){
    index1 <- index1 + 1
    chars <- strsplit(list[index1], "")[[1]]
    depth1 <- (which(chars == "\\")-1)/2
    id1 <- chars[length(chars)]
  }

  index2 <- start_2
  while (!(id2 %in% c( '|', '&', 'x'))){
    index2 <- index2 + 1
    chars <- strsplit(list[index2], "")[[1]]
    depth2 <- (which(chars == "\\")-1)/2
    id2 <- chars[length(chars)]
  }

  left_node <- Node$new("_")
  right_node <- Node$new("_")
  if (id1 == 'x'){
    chars <- strsplit(list[index1+1], "")[[1]]

    len <- length(chars)
    id1 <- chars[len]
    while (id1 == " "){
      len <- len - 1
      id1 <- chars[len]
    }
    if (chars[len-1] != " "){
      id1 <- paste(chars[len-1], chars[len], sep = '')
    }

    left_node <- Node$new(id1)
  }
  else{
    left_node <- tree_building(list, index1)
    # new_node$AddChildNode(next_node)
  }
  #
  if (id2 == 'x'){
    chars <- strsplit(list[index2+1], "")[[1]]

    len <- length(chars)
    id2 <- chars[len]
    while (id2 == " "){
      len <- len - 1
      id2 <- chars[len]
    }
    if (chars[len-1] != " "){
      id2 <- paste(chars[len-1], chars[len], sep = '')
    }

    right_node <- Node$new(id2)
    # new_node$AddChild(id2)
  }
  else {
    right_node <- tree_building(list, index2)
    # new_node$AddChildNode(next_node)
  }

  if (right_node$name == left_node$name){
    right_node$name <- paste(right_node$name, 2, sep = '')
  }

  new_node$AddChildNode(left_node)
  new_node$AddChildNode(right_node)
  # left_node$AddSiblingNode(right_node)

  return(new_node)
}

get_paths_from_gpr <- function(node){
  paths <- c()
  name <- node$name
  # print(name)
  name <- strsplit(name, '')[[1]][1]
  # if (length(name) == 0){
  #   return(c())
  # }

  if (name != '&' & name != '|'){
    paths <- c(paste('x', '[', node$name, ']', sep = ''))
    return(paths)
  }

  left_paths <- get_paths_from_gpr(node$children[[1]])
  right_paths <- get_paths_from_gpr(node$children[[2]])

  if (name == '&'){
    idx <- 0
    for (i in 1:length(left_paths)){
      for (j in 1:length(right_paths)){
        idx <- idx + 1
        paths[idx] <- list(c(unlist(left_paths[i]), unlist(right_paths[j])))
      }
    }
  }

  if (name == '|'){
    for (i in 1:length(left_paths)){
      paths[i] <- list(unlist(left_paths[i]))
    }

    for (i in 1:length(right_paths)){
      paths[length(left_paths) + i] <- list(unlist(right_paths[i]))
    }
  }

  return(paths)
}

get_all_gpr_paths <- function(gprRules){
  paths <- c()

  for (i in gprRules){
    paths <- c(paths, list(find_gpr_paths(i)))
  }

  return(paths)
}

additional_reaction_count <- function(gpr, gpr_paths){
  binary <- matrix(0, nrow = length(gpr), ncol = 1)
  for (i in 1:length(gpr)){
    if (nchar(gpr[i]) > 0){
      binary[i] <- 1
    }
  }

  paths <- matrix(0, nrow = length(gpr), ncol = 1)
  for (i in 1:length(paths)){
    if (nchar(gpr_paths[i][1]) > 0){
      paths[i] <- length(gpr_paths[[i]])
    }
  }

  return(paths - binary)
}
