# Functions for NGL based viewing of bio3d objects
# 2025-02-14   (11:44:38 PST on Fri, Feb 14)

#' Convert a bio3d PDB object for NGLVieweR
#'
#' Function to take a bio3d structure and use in the NGLVieweR package.
#' 
#' @param pdb a bio3d pdb object as obtained from `read.pdb()`.
#'
#' @return a single element character vector with return characters as required by `NGLVieweR::NGLVieweR()` function with `format=pdb` option.   
#'
#' @author Barry Grant, \email{bjgrant@@ucsd.edu}
#' 
#' @seealso \code{view.pdb()} which uses this function internally, \code{NGLVieweR::NGLVieweR()}, \code{bio3d::read.pdb()}.
#' 
#' @export
#'
#' @examples
#' pdb <- read.pdb("5p21")
#' NGLVieweR::NGLVieweR(pdb2ngl(pdb), format="pdb") |>
#'   NGLVieweR::addRepresentation("cartoon")
#'   
pdb2ngl <- function(pdb) {
  temp <- tempfile()
  on.exit(unlink(temp))
  bio3d::write.pdb(pdb = pdb, file = temp)
  return( paste(readLines(temp), collapse = "\n") )
}



#' Quick interactive PDB viewing using NGLVieweR
#'
#' Generate a quick webGL structure overview of `pdb` objects 
#'  with a number of simple defaults. The returned NGLVieweR object can 
#'  be further added to for custom interactive visualizations.
#'  
#'  The purpose of this function is to quickly view a given PDB
#'  structure object without having to write multiple lines of 
#'  NGLVieweR code. 
#'  
#'  The extra argument `highlight` takes an 
#'   `atom.select()` object to highlight atoms as a given `highlight.style` 
#'   with options including 'ball+stick', 'spacefill', 'licorice', 'line',
#'   'surface' or 'ribbon'.
#'   
#'  To-Do: 
#'  Currently the function does not check for bio3d or NGLVieweR 
#'  availability or work with multi-model structures, trajectory 
#'  or `pdbs` objects.
#' 
#'
#' @param pdb PDB structure object as obtained from `read.pdb()`.
#' @param ligand logical, if TRUE ligands will be rendered as atom colored licorice.
#' @param chain.colors optional color vector for chain based protein cartoon colors. By default `vmd.colors()` are used.
#' @param backgroundColor set the display area background color.
#' @param highlight an optional `atom.select()` object for highlighting.
#' @param highlight.style the representation style to use for selected `highlight` atoms 
#'
#' @return an **NGLVieweR** display object that can be displayed or further added to using `NGLVieweR::addRepresentation()` and friends.
#' 
#' @author Barry Grant, \email{bjgrant@@ucsd.edu}
#' 
#' @seealso \code{pdb2ngl()}, \code{view.pdbs()}, \code{NGLVieweR::NGLVieweR()}, \code{bio3d::read.pdb()}
#' 
#' @export
#'
#' @examples
#'  pdb <- read.pdb("1hsg")
#'  view.pdb(pdb)
#'  view.pdb(pdb, ligand=FALSE, chain.colors=c("pink","aquamarine"))
#'  
#'  #ras <- read.pdb("5p21")
#'  #view.pdb(ras) 
#'  
#'  sele <- atom.select(pdb, resno=c(25, 50))
#'  view.pdb(pdb, highlight = sele, 
#'         chain.colors = c("navy","orange"),
#'         backgroundColor = "pink",
#'         highlight.style = "spacefill")

view.pdb <- function(pdb, ligand=TRUE,
                     chain.colors=vmd_colors(),
                     backgroundColor = "white", # model = NULL
                     highlight=NULL,
                     highlight.style="ball+stick") {
  
  # Find ligand resid/resn
  lig <- atom.select(pdb, "ligand", value=TRUE)
  lig.resid <- unique(lig$atom$resid)
  
  # Do we have multiple protein chains
  #chains <- unique(pdb$atom$chain)
  chains <- paste(":", unique(atom.select(pdb, "protein", value=TRUE)$atom$chain),sep="")
  multi.chain <- length( chains ) > 1
  
  model <- NGLVieweR::NGLVieweR( pdb2ngl(pdb), format="pdb") |>
    NGLVieweR::stageParameters(backgroundColor = backgroundColor)
  
  if(!multi.chain) {
    # Color N-C term spectrum
    model <- model |> NGLVieweR::addRepresentation("cartoon",
                       param = list(colorScheme = "residueindex") )
  } else {
    ## Color by chain
    ## Not clear how to set custom colors when using NGLs colorScheme = "chainid"
    ## so we instead loop through each chain with a for loop here: 
    #model <- model |> NGLVieweR::addRepresentation("cartoon",
    #    param = list(colorScheme = "chainid") )
    
    # Ensure chain.colors has no names set
    names(chain.colors) <- NULL
    
    for(j in 1:length(chains)) {
      model <- model |> 
        addRepresentation("cartoon", 
         param = list(color = chain.colors[j], sele=chains[j])) 
    }
    
  }
  
  # Add ligand as licorice/"spacefill"
  if(length(lig.resid) > 0) {
    if(ligand) {  
      model <- model |> NGLVieweR::addRepresentation("licorice", 
                       param=list(sele=paste(lig.resid, collapse=" "),
                        radius=0.3) )
    } else {
      if(!is.null(lig.resid)) {
        message("Potential ligands found but not displayed, use ligand=T to view")
      }
    }
  }
  # Highlight atoms based on atom.select() object
  if(!is.null(highlight)) {
    #highlight <- atom.select(pdb, resno=90)
    highlight.eleno <- pdb$atom[highlight$atom,]$eleno
    eleno <- paste(paste("@", highlight.eleno, sep=""),collapse = " ")
    model <- model |> 
      addRepresentation(highlight.style, param=list(sele=eleno))
  }
  
  return(model)
}



#' Quick interactive PDBS object viewing using NGLVieweR
#'
#' @param pdbs a multi-structure `pdbs` object as obtained from `pdbaln()`, `read.fasta.pdb()`, etc.
#' @param colors a vector of colors for each structure. If NULL then the output of `vmd.colors()` is used.
#' @param representation the representation style, usefull values are lines, tube and cartoon. 
#' @param backgroundColor set the display area background color.
#'
#' @returns an **NGLVieweR** display object that can be displayed or further added to using `NGLVieweR::addRepresentation()` and friends.
#' 
#' @author Barry Grant, \email{bjgrant@@ucsd.edu}
#' 
#' @seealso \code{view.pdb()}, \code{pdb2ngl()}, \code{NGLVieweR::NGLVieweR()}, \code{bio3d::read.pdb()}
#' 
#' @export
#'
#' @examples
#'  # pth <- "~/Desktop/courses/BIMM143/class10/pdbs/split_chain/"
#'  # files <- list.files(path=pth, full.names = T)
#'  # pdbs <- pdbaln(files, fit=T, exefile="msa")
#'  
#'  view.pdbs(pdbs, representation = "cartoon")
#'  view.pdbs(pdbs, colors = c("red","blue") )
#'  
view.pdbs <- function(pdbs, 
                      colors=NULL, 
                      representation="line",
                      backgroundColor = "white"){
  
  ## Convert to individual PDB objects
  all.pdbs <- pdbs2pdb(pdbs)
  n.pdbs <- length(pdbs$id)
  
  # Setup default function params
  #colors=NULL #vmd_colors()
  #backgroundColor = "white"
  #highlight=NULL
  #highlight.style="ball+stick"
  #representation="line" # "cartoon"
  
  if(is.null(colors)) {
    colors <- vmd_colors()
  }
  if(length(colors) < n.pdbs) {
    warning("Not ennough distinct colors for each structure, recycling")
    colors <- rep(colors, length.out=n.pdbs)
  }
  # Names cause JSON/NGL problems
  names(colors) <- NULL
  
  # Setup stage and viewer with first pdb
  model <- NGLVieweR::NGLVieweR( pdb2ngl( all.pdbs[[1]] ), format="pdb") |>
    NGLVieweR::stageParameters(backgroundColor = backgroundColor) |>
    addRepresentation(representation, param = list(color = colors[1]))
  
  # Work through each remaining structure/pdb
  for(k in 2:n.pdbs) {
    model <- model |> 
      addStructure(pdb2ngl(all.pdbs[[k]]), format="pdb") |>
      addRepresentation(representation, param = list(color = colors[k]))
  }
  return(model)
}



## -- OLDER attempt without intermediate files ---##


#' Alternate convert a bio3d PDB object for NGLVieweR
#'
#' Function to take a bio3d structure and use in the NGLVieweR package.
#' 
#' @param pdb a bio3d pdb object as obtained from `read.pdb()`.
#'
#' @return a single element character vector with return characters as required by `NGLVieweR::NGLVieweR()` function with `format=pdb` option.   
#'
#' @author Barry Grant, \email{bjgrant@@ucsd.edu}
#' 
#' @seealso \code{view.pdb()} which uses this function internally, \code{NGLVieweR::NGLVieweR()}, \code{bio3d::read.pdb()}.
#' 
#' @export
#'
#' @examples
#' pdb <- read.pdb("5p21")
#' NGLVieweR::NGLVieweR(b2n(pdb), format="pdb") |>
#'   NGLVieweR::addRepresentation("cartoon")
#'   
b2n <- function(pdb) {
  # Format each row of pdb$atom to PDB 'ATOM' record lines...
  pdb_lines <- apply(pdb$atom, 1, function(row) {
    sprintf("ATOM  %5d %-4s %3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f          %-2s",
            as.integer(row["eleno"]), row["elety"], 
            row["resid"], row["chain"], as.integer(row["resno"]), 
            as.numeric(row["x"]), as.numeric(row["y"]), as.numeric(row["z"]),
            as.numeric(row["o"]), as.numeric(row["b"]), row["elesy"])
  })
  
  ## To-Do: what about multi-model files?
  
  # Combine all lines into a single element 
  #  character string with rtn as used by NGLVieweR()
  return(paste(pdb_lines, collapse = "\n"))
  
}

# Multi model attempt

bm2n <- function(pdb) {
  # Format each row of pdb$atom to PDB 'ATOM' record lines...
  pdb_lines <- apply(pdb$atom, 1, function(row) {
    sprintf("ATOM  %5d %-4s %3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f          %-2s",
            as.integer(row["eleno"]), row["elety"], 
            row["resid"], row["chain"], as.integer(row["resno"]), 
            as.numeric(row["x"]), as.numeric(row["y"]), as.numeric(row["z"]),
            as.numeric(row["o"]), as.numeric(row["b"]), row["elesy"])
  })
  
  
  ## To-Do: what about multi-model files?
  if(nrow(pdb$xyz) > 1) {
    cat("*** Multi-model file ***\n")
    pdb_lines <- c("MODEL        1", pdb_lines, "ENDMDL")
    xinds <- seq(from=1, to=ncol(pdb$xyz), by=3)
    
    # Do model 2 to last model
    for(i in 2:nrow(pdb$xyz)) {
      cat(i,"\n")
      pdb_lines <- c(pdb_lines, 
                     paste0("MODEL      ",i))
      # \nENDMDL\nMODEL       34\n
      pdb$atom$x <- pdb$xyz[i, xinds]
      pdb$atom$y <- pdb$xyz[i, (xinds+1)]
      pdb$atom$z <- pdb$xyz[i, (xinds+2)]
      
      pdb_lines <- c(pdb_lines, apply(pdb$atom, 1, function(row) {
        sprintf("ATOM  %5d %-4s %3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f          %-2s",
                as.integer(row["eleno"]), row["elety"], 
                row["resid"], row["chain"], as.integer(row["resno"]), 
                as.numeric(row["x"]), as.numeric(row["y"]), as.numeric(row["z"]),
                as.numeric(row["o"]), as.numeric(row["b"]), row["elesy"])
      }) )
      pdb_lines <- c(pdb_lines, "ENDMDL")
    }
  }
  
  # Combine all lines into a single element 
  #  character string with rtn as used by NGLVieweR()
  return(paste(pdb_lines, collapse = "\n"))
  
}
