##' Read acoustic HAC raw data
##'
##' The HAC data format is a binary format containing so-called \code{tuples}.
##' A tuple can hold various sorts of information depending on the tuple type.
##' For instance tuples exist to specify positions, echosounder information and
##' acoustic signal data etc.
##' This R package can read, write and subset the HAC data format.
##'
##' See the description of the ICES HAC standard data exchange format,
##' version 1.60.
##' 
##' @name readHAC-package
##' @aliases readHAC-package
##' @docType package
##' @references
##'     McQuinn, Ian H., et al. Description of the ICES HAC standard
##'     data exchange format, version 1.60. Conseil international pour
##'     l'exploration de la mer, 2005.
##'     \url{http://biblio.uqar.ca/archives/30005500.pdf}
##' @example demo/intro.R
NULL

## Get all needed HAC definition data:
## * tableList:          List of definition tables.
## * tuple2text:         Lookup code -> description text
## * tuple2identifier:   Lookup code -> identifier offset
## * tuple2parent:       Lookup code -> parent offset
## * known.tuples:       Vector of known tuples

## Read tuple definitions
file <- system.file("hacdef.dcf",package="readHAC")
zz <- read.dcf(file,all=TRUE)
zz$offset[zz$offset=="..."] <- NA
zz$offset <- as.numeric(zz$offset)
## FIXME: Some tuples e.g. 10040 have a field "Optional field: Space".
## It is unclear how to tell whether it is in use or not...
## For now discard this field (but beware that if space-field is in use
## it will be read as an additional sample)
zz <- subset(zz, field != "Optional field: Space")
## FIXME: Private tuple (65397) has "length: multiple of 4 bytes"
##        For now un-supported :
zz$length[zz$length=="multiple of 4 bytes"] <- NA
zz$length <- as.numeric(zz$length)
head <- zz$header[!is.na(zz$header)]
regexpr <- ".*\\. (.*)\\((.*)\\).*"
tup <- as.integer(gsub(regexpr,"\\2",head))
txt <- gsub(regexpr,"\\1",head)
tuple2text <- txt; names(tuple2text) <- tup
df <- data.frame(txt=txt,tup=tup,stringsAsFactors=FALSE)
tableList <- base::split(zz,cumsum(!is.na(zz$header))) ## All 30 tables
names(tableList) <- df$tup                       ## Name = tuple code
rm("df","file","head","regexpr","tup","txt","zz")

## Interpretations
## e.g.: Repetition patterns, unique identifier, parent identifier etc.
interpretations <- function(df){
  df$rep <- FALSE
  i <- grep("\\#2",df$field)
  if(length(i)>0){
    df$content[i] <- "idem"
    df$field <- gsub(" \\(.*\\#[12]\\)","",df$field)
  }
  ## Tables with "idem"
  idem <- grep("[iI]dem",df$content)
  if(length(idem)>0){
    stopifnot(all(diff(idem))==1)
    df$rep[idem-length(idem)] <- TRUE ## Repeat previous entries
    df <- df[-idem,]                  ## Remove "idem" entries
  }
  df
}

tableList <- lapply(tableList,interpretations)
## Fast lookup tables
## tuple2identifier
## tuple2parent
## === List of known tuple codes:
## > cat(deparse(as.integer(df$tup)))
known.tuples <- c(20L, 41L, 42L, 100L, 200L, 210L, 901L, 1000L, 1001L,
2000L, 2001L, 2002L, 2100L, 4000L, 9001L, 10000L, 10001L, 10010L,
10011L, 10030L, 10031L, 10040L, 10090L, 10100L, 10140L, 10142L,
11000L, 65397L, 65534L, 65535L)

## Interpret unit field
## x <- unique(zz$unit)
parseUnit <- function(def,datatype){
  x <- def$unit
  regexpr <- "([0-9|\\.]*)[ ]*([a-z|A-Z|/]*).*"
  ## Unit field can be a semicolon separated string of the form "0th
  ## unit;1st unit;2nd unit;3rd unit" etc "datatype" determines the
  ## number.
  if(is.character(datatype)){
      y <- strsplit(x, ";")
      x <- sapply(y, function(y)findBestMatch(datatype, y))
  }
  ## Cleanup string:
  ## Example:
  ##   x <- c("0.000001 deg precision: ~ 2 cm",
  ##          "description: 10 dB")
  ## Should become:
  ##        c("0.000001 deg",
  ##          "10 dB")
  x <- sub("[ ]*", "", x) ## Remove initial whitespace
  x <- sub("^[a-z|A-Z].*:[ ]*", "", x) ## E.g: "description: 10 dB"
  ans <- data.frame(mult=as.numeric(sub(regexpr,"\\1",x)),
                    unit=sub(regexpr,"\\2",x),stringsAsFactors=FALSE)
  ans$mult[ans$mult==0] <- NA
  ans$unit[tolower(ans$unit)=="unitless"] <- ""
  ans
}
## x <- unique(zz$content)
parseContent <- function(def){
  x <- def$content
  spl <- strsplit(x," =[ |\n]")
  myf <- function(x){
    if(length(x)==1)return(NA)
    re <- "(.*) ([0-9]*)$"
    ans <- gsub(re,"\\1",x)[-1]
    ans <- gsub("\n"," ",ans)
    nm <- gsub(re,"\\2",x)
    names(ans) <- nm[-length(nm)]
    ans
  }
  lapply(spl,myf)
}
## Add units to tuple
addUnits <- function(x){
  type <- unique(x[["Tuple type"]])
  stopifnot(length(type)==1)
  ## Hack to get datatype
  i <- grep("software channel identifier",names(x),ignore.case=TRUE)
  if(length(i)>0){
    channel <- unique(x[[i]])
    stopifnot(length(channel)<=1)
    datatype <- channel2unitname(x)[as.character(channel)]
  } else datatype <- 0
  def <- tableList[[as.character(type)]]
  df <- parseUnit(def,datatype)
  co <- parseContent(def)
  df$haslevels <- sapply(co,length)>1
  convert <- function(i){
    df <- df[i,]
    x <- x[[i]]
    if(!is.na(df$mult)){
      x[] <- x[]*df$mult
    } else if(df$haslevels){
      x <- (co[[i]][as.character(x)])
    }
    x
  }
  ans <- lapply(seq(x),convert)
  names(ans) <- paste(names(x)," [",df$unit,"]",sep="")
  class(ans) <- "tuple"
  ans
}
## Utility to help finding the right unit.
## Given a one-length character x find the best match in
## a character vector y.
## Example:
## x <- "Sv [Volume backscattering strength in dB]"
## y <- c("0.001 volts", "Sv or TS: 0.01 dB")
findBestMatch <- function(x, y){
    stopifnot(length(x)==1 &&
              is.character(x) &&
              is.character(y))
    doSplit <- function(x){
        strsplit(gsub("\\[*\\]*\\(*\\)*","",x),"[ ]+")[[1]]
    }
    splx <- doSplit(x)
    score <- function(y){
        sply <- doSplit(y)
        sum(splx %in% sply) +
            sum(tolower(splx) %in% tolower(sply))
    }
    s <- sapply(y, score)
    y[which.max(s)]
}

## Quotes from ICES HAC manual:
## ============================
##
## All tuple types contain five required fields in the following order:
##   (1) the tuple data size (totallength of the data fields),
##   (2) the tuple type code, or file tag,
##   (3) the data fields,
##   (4) the tuple attribute field (e.g. original tuple, edited
##       tuple), and
##   (5) the tuple backlink, which gives the tuple size (i.e. 10 bytes
##       longer than the tuple data size).
##
## A software application is defined as HAC compatible if it can read
## and/or write, and use a minimum number of commonly used basic
## tuples following the HAC syntax rules as outlined in this
## document. These tuple numbers are:
##
## * 20, 100, 200, 901, 1000, 2000, 2001, 2002, 9001, 10000, 10001,
## 10010, 10011, 10030, 10031, 10040, 10100, 65534, 65535
##
## Therefore a file
##
##   (1) must start with the code 172 (=hexadecimal 0xAC, for
##       ACoustics) stored in a ULONG2 word in the HAC signature
##       tuple, to identify the byte encoding mode of the computer
##       platform, and
##
##   (2) must contain at least the seven minimum tuple types. The
##       first tuple in an HAC file must be the HAC signature
##       tuple. The last tuple must be the End of file tuple.
##
## Position tuple
## This tuple type stores the geographic position of the transducer
## deployment platform. The data collection rate of position data can
## be independent of the ping rate. [Reserved tuple type codes: 20 -
## 29].


## =============================== PART I: PARSE ICES DOCUMENTATION
## get Tuple definitions
## Each Table defines a tuple.
## ============== At page 3 of the manual we have: 
## ...... 
## To solve these ambiguities, it was decided to fix the Remark field
## lengths of tuples 901 and 9001 at 40 and 100 bytes, respectively,
## and at 20 bytes for tuple 2002.
## ......
## ==============
## We have to patch the tables: 901 has 100 bytes and 9001 has 40.
## Hm, not solving the problem. In practice the character field can
## have variable length.
if(TRUE){ ## DO_PATCH
  DO_PATCH <- TRUE
  patchTable <- function(x,n){
    i <- grep("Remarks",x$field)
    x$length[i] <- n
    cs <- cumsum(c(0,x$length))
    cs <- cs[-length(cs)]
    x$offset <- as.integer(cs)
    x
  }
  ## tableList[["901"]] <- patchTable(tableList[["901"]],40)
  ## tableList[["9001"]] <- patchTable(tableList[["9001"]],100)
}

## TODO: we can do a lot more:
## * Data values: grep for "Idem". Then we know the data pattern.


## =============================== PART II: READ EXAMPLE FILE

## TODO: Handle NA values correctly: Manual says extreme values are
## used to code NAs.

## Extractor functions (vectorized) from binary format.
USHORT <- function(raw){
  dim(raw) <- c(2L, length(raw) / 2L)
  storage.mode(raw) <- "integer"
  as.integer( t( c(1, 256) ) %*% raw )
}
ULONG <- function(raw){
  ## Largest ulong is 2^32-1. Greater than largest R integer (2^31-1)!
  ## ==> MUST represent ULONG as double !
  dim(raw) <- c(4L, length(raw) / 4L)
  storage.mode(raw) <- "integer"
  as.numeric( t( c(1, 256, 65536, 16777216) ) %*% raw )
}
CHAR <- function(raw){
  nulmatch <- match(as.raw(0),raw)
  if(is.na(nulmatch)) return(rawToChar(raw))
  i <- seq.int(length.out = nulmatch) ## Let \0 interrupt character
  rawToChar(raw[i])
}
SHORT <- function(raw){
  ## Implementation of "twos complement".
  tmp <- USHORT(raw)
  tmp - (tmp >= 2^15) * 2^16
}
LONG <- function(raw){
  tmp <- ULONG(raw)
  tmp - (tmp >= 2^31) * 2^32
}
RAW <- function(raw)raw

## Interpretation of n-bit RLE compression (note: we count bits from 1
## to n as opposed to C-style 0 to n-1).
## Highest bit (n) codes the type of sample:
##   * If 1 the remaining n-1 bits code an unsigned integer giving the
##   number (minus-1) of zero-samples.
##   * If 0 the remaining n-1 bits code a signed integer giving *one*
##   sample.
SHORTrle <- function(raw){
  tmp <- USHORT(raw)
  ## Get highest bit (16).
  b16 <- tmp >= 2^15
  ## Get second highest bit (15).
  b15 <- (tmp - b16 * 2^15) >= 2^14
  ## Case b16=1: Get *unsigned* 15 bit integer
  u <- tmp - b16 * 2^15
  ## Case b16=0: Get *signed* 15 bit integer
  s <- u - b15 * 2^15
  ## Do RLE decompression
  value <- ifelse(b16,   0, s)
  repl  <- ifelse(b16, u+1, 1)
  rep(value, repl)
}
LONGrle <- function(raw){
  tmp <- ULONG(raw)
  ## Get highest bit (32).
  b32 <- tmp >= 2^31
  ## Get second highest bit (31).
  b31 <- (tmp - b32 * 2^31) >= 2^30
  ## Case b32=1: Get *unsigned* 31 bit integer
  u <- tmp - b32 * 2^31
  ## Case b32=0: Get *signed* 31 bit integer
  s <- u - b31 * 2^31
  ## Do RLE decompression
  value <- ifelse(b32,   0, s)
  repl  <- ifelse(b32, u+1, 1)
  rep(value, repl)
}

## Parse tuple
## TODO:
## * Data values. (Idem gives the repetition pattern).
## * We know all tuples end with 8 bytes (attribute and backlink fields).
## * We probably want parsing in two steps:
##   1. (Machine readable) A rough parse where value types are USHORT,ULONG,SHORT,LONG,CHAR as now and
##   2. (Human readable)   A further step where units are added and e.g. time is calculated (Time+Time fraction etc.)
## * Performance consideration: Tuple parsing is probably a bottleneck. Maybe R will do okay since
##   data values can be processed vectorized. However, the number of tuples is high (>30000). Alternative is
##   to "merge" many compatible tuples to a "vectorized tuple". ("Compatible Tuples" here means same raw-length and
##   tuple-code).
## NOTE: A merge of compatible tuples is probably a good idea!!:
parseBinaryTuple <- function(tup){
  code <- USHORT(tup[5:6])
  def <- tableList[[as.character(code)]]
  if(is.null(def)){
    ans <- list("Error:"="This tuple-code does not exist in HAC definition manual!",
                "Tuple type"=code)
    class(ans) <- "tuple"
    return(ans)
  }

  ## Variable character length bug
  if(DO_PATCH){
    i <- which(def$format=="CHAR")
    if(length(i)==1){
      ## Only one tuple has length=4 (2100) and current fix wont work
      ## (not a problem if '2100' uses fixed length CHARs - which
      ## seems to be the case in practice)
      newlen <- ULONG(tup[1:4])+10-sum(def$length[-i]) ## This must be CHAR length !
      def <- patchTable(def,newlen)
    }
  }
  
  ## Allow many compatible tuples to be parsed simultaneously
  size <- ULONG(tup[1:4])+10 ## Size of this tuple type (bytes)
  dim(tup) <- c(size,length(tup)/size)

  ## Fields "Tuple attribute" and "Backlink" always occupy the final 8 bytes:
  def$offset[c(nrow(def)-1,nrow(def))] <- c(size-8,size-4)

  ## If "def$length" still have empty fields then plug in something...
  if(any(is.na(def$length))){
    na <- is.na(def$length)
    stopifnot(sum(na)==1)
    def$length[na] <- size-sum(def$length[!na])
    def$offset[na] <- def$offset[which(na)-1]+def$length[which(na)-1]
    def$format[na] <- "RAW"
  }
  
  ## DATA-blocks
  ## * Block size = sum(def$length[def$rep])
  ## * Begin: Offset        (bytes)
  ## * End:   length(tup)-8 (bytes)
  ## How many blocks are there space for?:
  blocksize <- sum(def$length[def$rep])
  if(blocksize>0){ ## Have data
    offset <- def$offset[def$rep][1]+1
    nblocks <- floor((length(tup[,1])-8-offset)/blocksize)
  }
  sequence <- function(offset,length,rep=FALSE){
    ans <- seq.int(offset,length.out=length)
    if(rep){
      shift <- seq.int(from=0,by=blocksize,length.out=nblocks)
      ans <- outer(ans,shift,"+")
      dim(ans) <- NULL
    }
    ans
  }
  args <- Map(function(x,y,z)tup[sequence(x,y,z),],def$offset+1,def$length,def$rep)
  funs <- lapply(def$format,get, envir = asNamespace("readHAC"), inherits = FALSE)
  ans <- Map(function(x,f)f(x),args,funs)
  setdim <- function(x){
    if(ncol(tup)>1 & length(x)>ncol(tup))dim(x) <- c(length(x)/ncol(tup),ncol(tup))
    x
  }
  ans <- lapply(ans,setdim)
  trim <- function(x)gsub("[ ]*$","",gsub("^[ ]*","",x))
  names(ans) <- trim(def$field)
  class(ans) <- "tuple"
  ans
}

print.tuple <- function(x,show=2,...){
  m <- max(nchar(names(x)))+1
  addspace <- function(qw,n=m-nchar(qw))paste(qw,paste(rep(" ",n),collapse=""))
  myformat <- function(x){
    if(is.vector(x)){
      if(length(x)>15){
        txt <- paste("[length=",length(x),"]",sep="")
        x <- c(txt,head(x,show),"...",tail(x,show))
      }        
    }
    if(is.matrix(x)){
      txt <- paste("[dim=",paste(dim(x),collapse="x"),"]",sep="")
      space <- paste(rep(" ",m+nchar(txt)+1),collapse="")
      x <- c(txt,
             head(x[1,],show),"...",tail(x[1,],show),"\n",
             space,"...\n",
             space,head(x[nrow(x),],show),"...",tail(x[nrow(x),],show))
    }
    x
  }
  for(i in 1:length(x)){
    cat(addspace(names(x)[i]))
    cat(" ")
    cat(myformat(x[[i]]))
    cat("\n")
  }
}

getIdentifierOffset <- function(regexpr="software channel"){
  myf <- function(x,name){
    i <- grep(regexpr,x$field,ignore.case=TRUE)
    if(length(i)>1)stop("More than one match!")
    if(length(i)==0)return(NA)
    x[[name]][i]
  }
  offset <- sapply(tableList,myf,"offset")
  format <- sapply(tableList,myf,"format")
  length <- sapply(tableList,myf,"length")
  field <- sapply(tableList,myf,"field")
  shave <- function(x)unique(x[!is.na(x)])
  attr(offset,"format") <- shave(format)
  attr(offset,"length") <- shave(length)
  attr(offset,"field") <- field
  offset
}

## Now, we can create a "map" to aid extracting relevant tuples.
## A map here is a data.frame with a row for each tuple.
## A HAC subset should be performed in two steps:
## 1. Perform the subsetting on the "map" first.
## 2. Then on the actual HAC byte stream.
## NOTE: 1. is much smaller than 2!
identifierRegexpr <- list(
                          softwarechannel="software channel identifier",
                          echosounder="echo.*sounder document identifier",
                          subchannel="sub-channel identifier",
                          depthsensor="depth sensor identifier",
                          distancesensor="distance sensor identifier",
                          attitudesensor="attitude sensor identifier",
                          time="time cpu ansi c standard",
                          timefraction="time fraction",
                          typeofdata="type of data|Data type"
                          )
identifierOffsetList <- lapply(identifierRegexpr,getIdentifierOffset)

## HAC class
## =========
## * Simply a data.frame with a row for each tuple.
## * With attribute "binary" defining the binary raw data.
##   (Binary data is wrapped into an environment to avoid
##    deep copy when e.g. subsetting a HAC class. )
binary <- function(x){
  m <- attr(x,"binary")$data
  j <- unlist(Map("seq.int",from=x$pointer,length.out=x$length))
  m[,j,drop=FALSE]
}
"binary<-" <- function(x,value){
  e <- new.env()
  e$data <- value
  attr(x,"binary") <- e
  x
}
channel2datatype <- function(x)attr(x,"binary")$channel2datatype
channel2unitname <- function(x)attr(x,"binary")$channel2unitname

## ---------------------------------------------------------------------------
##' Read raw HAC data file
##'
##' This function reads the binary HAC format and locates the tuples.
##' @title Read HAC data into R.
##' @param file File to read.
##' @return HAC object.
## ---------------------------------------------------------------------------
readHAC <- function(file){
  size.bytes <- file.info(file)$size
  ## Read both as raw and integer:
  m <- readBin(file,raw(),size.bytes,size=1,endian="little")
  dim(m) <- c(4,length(m)/4)
  ## Integers (ULONG) indexed by columns 
  c <- readBin(file,integer(),size.bytes/4,size=4,endian="little")
  stopifnot(ncol(m)==length(c))
  stopifnot(c[1]==172L)  ## Must start like this
  ## Find tuple pointers
  ## NOTE: * 32-bit (integer) pointers
  ##       * Indentifies columns of m
  ##       * Elementwise 8-bit (raw) pointers: 4*p-3 
  pointers <- function(){ ## of c
    p <- integer(length(c))
    p[1] <- 2L ## First pointer
    m[1:2,p[1]+1] ## Tuple-type
    for(k in 1:length(p)){
      if(p[k]>=length(c))break;
      p[k+1] <- p[k]+(c[p[k]]+10L)/4L
    }
    p <- p[1:k]
    p
  }
  p <- pointers()
  types <- function(){ ## Of m and p
    mat <- m[1:2,p[-length(p)]+1]
    USHORT(mat)
  }
  typ <- types()
  identifier <- function(idoff){ ## Return identifier vector (by tuple)
    offset <- idoff[as.character(typ)]
    p0 <- 4*p[-length(p)]-3+offset         ## RAW-pointers of identifier
    length <- attr(idoff,"length")         ## Sizeof pointer
    mat <- outer(seq(0,length-1),p0,"+")
    dim(mat) <- NULL
    fun <- get(attr(idoff,"format"), mode="function")
    ans <- fun(m[mat])
    ans[is.na(p0)] <- NA
    ans
  }
  identifiers <- data.frame(lapply(identifierOffsetList,identifier))
  map <- cbind( data.frame(type=typ), identifiers )
  map$pointer <- head(p,-1)
  map$length <- diff(p)
  binary(map) <- m
  class(map) <- c("HAC",class(map))

  ## Workaround: The unit of the key data tuples are stored in other tuples !
  ## E.g. volts versus dB must be looked up in the parent channel of a ping tuple.
  ## This is a problem when subsetting like "subset(x,type==10000)"...
  ## Lets store the channel to datatype information as part of the HAC object:
  env <- attr(map,"binary")
  env$channel2datatype <- tapply(map$typeofdata,map$softwarechannel,function(x)x[1])
  ## Also store the 'unitname' corresponding to the datatype (number)
  channel2type <- tapply(map$type,map$softwarechannel,function(x)x[1])
  myf <- function(datatype, type){
      def <- tableList[[as.character(type)]]
      i <- grep("Type of data|Data type", def$field)
      tab <- parseContent(def[i,])[[1]]
      structure(tab[as.character(datatype)], names=NULL)
  }
  env$channel2unitname <- unlist(Map(myf, env$channel2datatype, channel2type))

  map
}

## ---------------------------------------------------------------------------
##' Write raw HAC data file
##'
##' This function writes the binary HAC format. The output file begins with
##' "ac 00 00 00" followed by the binary tuples defined by the HAC object \code{x}.
##' Note that the function does not perform a check for mandatory tuples.
##' @title Write HAC binary data.
##' @param x HAC object
##' @param file File to write to.
##' @return NULL
## ---------------------------------------------------------------------------
writeHAC <- function(x,file){
  stopifnot(inherits(x,"HAC"))
  x <- c( as.raw(c(172,0,0,0)) , binary(x) )
  writeBin(x,file,size=1,endian="little")
}

##' Extract tuples of HAC object.
##'
##' Extract subset of tuples.
##' For instance x[1:2] extracts the first two tuples.
##' Alternatively the method can be indirectly invoked by
##' the \code{subset} function.
##' @title Extract tuples.
##' @param x HAC object
##' @param i Integer vector
##' @param ... Currently not used
##' @return HAC object
##' @rdname indexSubset
##' @examples
##' \dontshow{
##' hacfile <- system.file("hac","Hac-test_000001.hac",package="readHAC")
##' x <- readHAC(hacfile)
##' }
##' x[1:2]
##' subset(x, type == 10000)
##' split(x, x$type)
"[.HAC" <- function(x, i, ...){
  ans <- as.data.frame(x)[i,,drop=FALSE]
  class(ans) <- class(x)
  attr(ans,"binary") <- attr(x,"binary")
  ans
}

##' Parse binary HAC to a list of data values.
##'
##' HAC parsing can be performed for one or multiple tuples of
##' the same type and length. The binary tuples are translated
##' to data values according to the definition document.
##' @title Parse binary HAC.
##' @param hac Object of class \code{HAC} to be parsed.
##' @param split Force parsing of incompatiple tuples by first splitting the raw data?
##' @param split.by If \code{split=TRUE} then split by this factor.
##' @param units Convert to human readable units?
##' @return Object of class \code{tuple}.
parseHAC <- function(hac,split=FALSE,split.by=paste(hac$type,hac$length),units=TRUE){
  if(split){
    split.by <- substitute(split.by)
    fac <- eval(split.by,hac)
    spl <- split(hac,fac)
    ans <- lapply(spl,parseHAC,split=FALSE)
    names(ans) <- NULL
    ans <- unlist(ans,FALSE)
    class(ans) <- "tuple"
    return(ans)
  }
  if(nrow(hac)>1){
    all.equal <- function(x)diff(range(x))==0
    ok <- all.equal(hac$length) & all.equal(hac$type)
    if(!ok)stop("Vectorized parsing requires compatible tuple tupes (equal type and length)")
  }
  tup <- binary(hac)
  ans <- parseBinaryTuple(tup)
  if(units){
    attr(ans,"binary") <- attr(hac,"binary") ## Unit conversion needs extra information
    ans <- addUnits(ans)
    attr(ans,"binary") <- NULL
  }
  ans
}
print.HAC <- function(x,show.max=10,...){
  txt <- paste("Object of class 'HAC' containing",nrow(x),"tuples")
  txt2 <- gsub(".","=",txt)
  cat(txt,"\n")
  cat(txt2,"\n")
  cat("\n")
  cat("Tuple types:\n")
  df <- as.data.frame(table(x$type))
  names(df) <- c("type","count")
  tuple <- tuple2text[as.character(df$type)]
  names(tuple) <- NULL
  df <- cbind(data.frame(tuple=tuple),df)
  print(df)
  cat("\n")
  cat("Identifiers (to select tuples):\n")
  abbrev <- function(x){
    x <- x[!is.na(x)]
    ux <- unique(x)
    if(length(ux)>show.max){
      ux <- paste(min(ux),max(ux),sep=" - ")
    } else ux <- paste(sort(ux),collapse=", ")
    paste("[",ux,"]")
  }
  m <- max(nchar(names(x)))+1
  addspace <- function(qw,n=m-nchar(qw))paste(qw,paste(rep(" ",n),collapse=""))
  nm <- sapply(names(x),addspace)
  y <- sapply(x,abbrev)
  y <- paste("$",nm,"\t",y,"\n",sep="")
  for(z in y)cat(z)
}
