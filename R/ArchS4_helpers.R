
#Package environments
excludedTerms_env <- new.env(parent = emptyenv())
excludedTerms_env$co <- c("cultured cell","zygote","centrocyte","precursor cell","primary cultured cell","oocyte","stem cell","germ cell",
                          "cap cell","pp cell","somatic stem cell","tr1 cell","be cell","somatic cell","h minus","h plus","native cell",
                          "part of","description")
excludedTerms_env$to <- c("awn","amnion","berry","culture condition","cob","crop","fetus","flower","fruit","excretion","ink","larva","leaf",
                          "myotome","infundibulum","mantle","zygote","oviduct","adult","pith","pollen","pod","pons","node","rind","roe",
                          "larynx","seed","serum","tail","venom","vena cava","plant","mitotic cell","morula","somite","palate","wart",
                          "nucleus accumbens","ampulla","juvenile","pr cell","aro cell","leaflet","s2 cell","coculture","rl cell","groin",
                          "m cell","keloid","na cell","aril","g cell","paw","s9 cell","r1 cell","fin","orbit","null cell","scale","kle cell",
                          "vocal fold","rd cell","wood","awn","column","vertebra","polyp","cell culture","stoma","part of","produced by",
                          "related to","infected cell","culture medium","ink","kb cell","embryonic cell line","leaf","carcinoid cell","adult",
                          "primary culture","roe","seed","primary cell","rko cell","peer cell","t2 cell","na cell","cell lysate","scale",
                          "pos cell","tex cell","t1 cell","animal","stoma","amnion","ascites","bacteriod","berry","infected cell","calyx","cell culture",
                          "culture condition","clove","cob","cone","crop","crypt","cyst","egg","embryo",
                          "fetus","fruit","excretion","germ","germ cell","germ layer","ink","larva","leaf",
                          "zygote","nape","ags cell","adult","pod","spore","rind","root","roe","larynx",
                          "seed","femur","stem","style","tail","tuber","ureter","primary cell","venom",
                          "vena cava","plant","loin","morula","somite","utricle","wart","nucleus accumbens",
                          "oc cell","ampulla","juvenile","pr cell","siha cell","aro cell","pg cell",
                          "leaflet","groin","rch-acv cell","m cell","germinal center","hos cell","g cell",
                          "nccit cell","cell lysate","paw","hn cell","hbe cell","sebocyte","fin",
                          "null cell","scale","rd cell","mnng/hos cell","wood","awn","kgn cell","column",
                          "cerebral organ","2ftgh cell","part of","produced by","related to")
excludedTerms_env$age <- c("-","conception","wpc","gestation","wk","NA","embryo","adult","pcw","NES","weeks",
                           "<","post-infection","child","Differentiated", "\\b(or)\\b","n/a","--","onset",
                           "infant","months","[0-9]+\\s*d", "days","[0-9]+\\s*w", "old", "young", "about", "mean")

#' Title
#'
#' @param version The version of ArchS4 to retrieve. Currently, "v2.1.2" and "v11" are supported.
#' @param org The organism of which to retrieve to the data. Currently, "human" and "mouse" are supported.
#' @param path The path to store the h5 file
#'
#' @return NULL
#' @export
#'
getArchS4Data <- function(version = c("v2.1.2","v11"),
                          org = c("human","mouse"),
                          path){
  version <- match.arg(version,c("v2.1.2","v11"))
  org <- match.arg(org,c("human","mouse"))
  if(version == "v2.1.2"){
    url <- paste0("https://s3.dev.maayanlab.cloud/archs4/archs4_gene_",org,"_",
                  version,
                  ".h5")
  }else if(version == "v11"){
    url <- paste0("https://s3.amazonaws.com/mssm-seq-matrix/",org,"_matrix_v11.h5")
  }

  destination_file <- paste0(path,
                             ifelse(endsWith(path,"/"),"","/"),
                             "archs4_gene_",
                             org,
                             ".h5")
  overwrite <- TRUE
  if(file.exists(destination_file)){
    message(paste0("File already exits: ",destination_file))
    overwrite <- ""
    while(!(base::tolower(overwrite) %in% c("yes","no"))){
      overwrite <- readline(prompt = paste0("Do you want to overwrite it? (Yes/No) "))
    }
    overwrite <- base::tolower(overwrite) == "yes"
  }
  if(overwrite){
    utils::download.file(url = url,
                         destfile = destination_file,
                         quiet = FALSE,
                         mode = 'wb',
                         method = "curl")
  }
  return(destination_file)
}

#' Title
#'
#' @param archs4_file The path to the ArchS4 data file as obtained from getArchS4Data.
#' @param scCutoff Cutoff for probability that the sample is from single cell data. (Default: 0.1)
#' @param numAligned Cutoff for the number of aligned reads per sample. Only samples having at least this many reads will be selected.
#' @param org The organism of which to select the samples. Currently, "human" and "mouse" are supported. (Default: human)
#'
#' @return NULL
#' @export
#'
selectArchS4Samples <- function(archs4_file = NULL,
                                scCutoff = 0.1,
                                numAligned = 10^7,
                                org = c("human","mouse")
                                ){
  org_orig <- match.arg(org,c("human","mouse"))
  if(org_orig == "human"){
    org <- "Homo sapiens"
    taxid <- "9606"
  }else if(org_orig == "mouse"){
    org <- "Mus musculus"
    taxid <- "10090"
  }
  samples <- rhdf5::h5read(archs4_file, "meta/samples/geo_accession")
  scProb <- rhdf5::h5read(archs4_file, "meta/samples/singlecellprobability")
  readsTotal <- rhdf5::h5read(archs4_file, "meta/samples/readstotal")
  readsAligned <- rhdf5::h5read(archs4_file, "meta/samples/readsaligned")
  libSelection <- rhdf5::h5read(archs4_file, "meta/samples/library_selection")
  libStrategy <- rhdf5::h5read(archs4_file, "meta/samples/library_strategy")
  libSource <- rhdf5::h5read(archs4_file, "meta/samples/library_source")
  organism <- rhdf5::h5read(archs4_file, "meta/samples/organism_ch1")
  molecule <- rhdf5::h5read(archs4_file, "meta/samples/molecule_ch1")
  taxid <- rhdf5::h5read(archs4_file, "meta/samples/taxid_ch1")
  instrument <- rhdf5::h5read(archs4_file, "meta/samples/instrument_model")
  source <- rhdf5::h5read(archs4_file, "meta/samples/source_name_ch1")
  char <- rhdf5::h5read(archs4_file, "meta/samples/characteristics_ch1")
  title <- rhdf5::h5read(archs4_file, "meta/samples/title")
  series <- rhdf5::h5read(archs4_file, "meta/samples/series_id")
  geoid <- rhdf5::h5read(archs4_file,"meta/samples/geo_accession")

  idx_polyA <- base::which(scProb <= scCutoff &
                     readsAligned > numAligned &
                     libSelection == "cDNA" &
                     libStrategy == "RNA-Seq" &
                     libSource == "transcriptomic" &
                     organism == org &
                     molecule == "polyA RNA" &
                     taxid == taxid
                    )

  idx_total <- base::which(scProb <= scCutoff &
                     readsAligned >= numAligned &
                     libSelection == "cDNA" &
                     libStrategy == "RNA-Seq" &
                     libSource == "transcriptomic" &
                     organism == org &
                     molecule == "total RNA" &
                     taxid == taxid
                    )

  co_terms <- MultiTIMER::co$name
  co_terms <- co_terms[-1]
  co_terms <- unname(co_terms)

  to_terms <- MultiTIMER::to$name
  to_terms <- to_terms[-1]
  to_terms <- unname(to_terms)

  source_lower <- base::tolower(source)
  char_lower <- base::tolower(char)
  co_to_samples <- sapply(base::tolower(co_terms),function(x){base::which(grepl(x, source_lower, fixed = TRUE))})
  idx <- base::which(sapply(co_to_samples,base::length) > 0)
  co_to_samples <- co_to_samples[idx]

  co_to_samples_char <- sapply(base::tolower(co_terms),function(x){base::which(grepl(x, char_lower, fixed = TRUE))})
  idx <- base::which(sapply(co_to_samples_char,base::length) > 0)
  co_to_samples_char <- co_to_samples_char[idx]

  to_to_samples <- sapply(base::tolower(to_terms),function(x){base::which(grepl(x, source_lower, fixed = TRUE))})
  idx <- base::which(sapply(to_to_samples,base::length) > 0)
  to_to_samples <- to_to_samples[idx]

  to_to_samples_char <- sapply(base::tolower(to_terms),function(x){base::which(grepl(x, char_lower, fixed = TRUE))})
  idx <- base::which(sapply(to_to_samples_char,base::length) > 0)
  to_to_samples_char <- to_to_samples_char[idx]

  idx_co <- union(setdiff(unique(do.call("c",co_to_samples)),unique(do.call("c",co_to_samples[which((names(co_to_samples) %in% MultiTIMER::excludedTerms_env$co))]))),
                  setdiff(unique(do.call("c",co_to_samples_char)),unique(do.call("c",co_to_samples_char[which((names(co_to_samples_char) %in% MultiTIMER::excludedTerms_env$co))])))
                 )
  idx_to <- union(setdiff(unique(do.call("c",to_to_samples)),unique(do.call("c",to_to_samples[which((names(to_to_samples) %in% MultiTIMER::excludedTerms_env$to))]))),
                  setdiff(unique(do.call("c",to_to_samples_char)),unique(do.call("c",to_to_samples_char[which((names(to_to_samples_char) %in% MultiTIMER::excludedTerms_env$to))])))
                 )

  samples_to_consider <- intersect(idx_total,union(idx_to,idx_co))

  batchSize <- as.data.frame(table(series[samples_to_consider]))
  idx_batchGr2 <- which(series %in% batchSize$Var1[batchSize$Freq > 2])

  samples_to_consider <- intersect(samples_to_consider,idx_batchGr2)
  samples_to_consider_geoid <- geoid[samples_to_consider]

  #Select wild type samples
  controlsamples <- Reduce(union,list(which(grepl("[Cc]ontrol",char)),
                                      which(grepl("[Hh]ealthy",char)),
                                      which(grepl("[Cc]ontrol",source)),
                                      which(grepl("[Hh]ealthy",source)),
                                      which(grepl("WT",char)),
                                      which(grepl("WT",source)),
                                      which(grepl("[Ww]ild[- ]?[Tt]ype",char)),
                                      which(grepl("[Ww]ild[- ]?[Tt]ype",source)),
                                      which(grepl("[Uu]ntreated",char)),
                                      which(grepl("[Uu]ntreated",source)),
                                      which(grepl("[Nn]ormal",char)),
                                      which(grepl("[Nn]ormal",source)),
                                      which(grepl("[Tt]reatment[ ]?:[ ]?[Nn]one",char)),
                                      which(grepl("[Tt]reatment[ ]?:[ ]?[Nn]one",source))
                                     )
                           )
  controlsamples_geoid <- geoid[controlsamples]

  withAge <- which(grepl("\\b([Aa]ge)\\b",char))
  age_info <- stringr::str_extract(char[withAge],"(?<=\t|^)\\s*\\b([Aa]ge)\\b(.)*?(?=\t|$)")
  names(age_info) <- geoid[withAge]

  idx_excludedPhrases <- sapply(MultiTIMER::excludedTerms_env$age,function(x){which(grepl(x,age_info,ignore.case = T))})
  idx_excludedPhrases <- unique(do.call("c",idx_excludedPhrases))
  idx_excludedPhrases <- unique(c(idx_excludedPhrases,which(is.na(age_info))))

  idx_retained <- setdiff(seq(1,length(age_info)),idx_excludedPhrases)
  age_info[idx_retained]

  withAge_retained <- withAge[idx_retained]
  age_info_retained <- age_info[idx_retained]
  age_info_num_retained <- age_info_retained
  age_info_num_retained <- (sapply(age_info_num_retained,function(x){
    as.numeric(unlist(base::regmatches(x,
                                 gregexpr("[[:digit:]]+\\.*[[:digit:]]*",x))
    )      )[1]
  }))

  age_info_num_retained <- age_info_num_retained[which(!is.na(age_info_num_retained))]

  seriesIDs <- series[samples_to_consider]

  return(list(samplesToConsider = geoid[samples_to_consider],
              controlSamples = controlsamples_geoid,
              samplesWithAge = age_info_num_retained,
              seriesIDs = seriesIDs,
              parameters = list(archs4_file = archs4_file,
                                scCutoff = scCutoff,
                                numAligned = numAligned,
                                org = org_orig)
              ))
}

#' Title
#'
#' @param archs4File The path to the ArchS4 data file as obtained from getArchS4Data.
#' @param backgroundSamples All samples obtained from ArchS4 as returned in the \emph{samplesToConsider} field by \link{selectArchS4Samples}.
#' @param samplesWithAge The samples with corresponding age information as returned by \link{selectArchS4Samples}.
#' @param controlSamples The control samples as returned by \link{selectArchS4Samples}.
#'
#' @return List of training samples and associated age information
#' @export
#'
getTrainingSamples <- function(archs4File = NULL,
                               backgroundSamples = c(),
                               samplesWithAge = c(),
                               controlSamples = c()
                               ){
  trainSamples <- intersect(backgroundSamples,intersect(names(samplesWithAge),controlSamples))
  trainAge <- samplesWithAge[trainSamples]

  return(list(trainingSamples = trainSamples,
              trainingAge = trainAge
              ))
}

#' Title
#'
#' @param archs4File The path to the ArchS4 data file as obtained from getArchS4Data.
#' @param SamplesToConsider All samples obtained from ArchS4 as returned in the \emph{samplesToConsider} field by \link{selectArchS4Samples}.
#' @param seriesIDs The GEO series IDs as returned by \link{selectArchS4Samples}.
#' @param legacy Should the legacy version be run that filters out genes based on a fixed threshold of mean expression? (default: FALSE)
#'
#' @return List of training samples and associated age information
#' @export
#'
getSelectedSamplesArchS4 <- function(archs4File = NULL,
                                     SamplesToConsider = c(),
                                     seriesIDs = c(),
                                     legacy = FALSE){

  genes = rhdf5::h5read(archs4File, "meta/genes/genes")
  geoid <- rhdf5::h5read(archs4File,"meta/samples/geo_accession")
  idx <- match(SamplesToConsider,geoid)
  expression = t(rhdf5::h5read(archs4File, "data/expression", index=list(idx, 1:length(genes))))
  rhdf5::H5close()

  rownames(expression) = genes
  colnames(expression) = geoid[idx]

  #TMM normalization
  expression_dgelist <- edgeR::DGEList(expression)
  expression_dgelist <- edgeR::calcNormFactors(expression_dgelist,
                                               method = "TMMwsp")

  rm(expression)

  expression_tmm_cpm <- edgeR::cpm(expression_dgelist,
                                   normalized.lib.sizes = T,
                                   log = T)

  rm(expression_dgelist)

  if(legacy){
    geExpMeans <- rowMeans(expression_tmm_cpm)
    idx_1 <- names(which(geExpMeans > -3.8))
    expression_tmm_cpm <- expression_tmm_cpm[idx_1,]
  }

  series <- seriesIDs[idx]

  batchid = match(series, unique(series))
  correctedExpression <- sva::ComBat(dat=expression_tmm_cpm, batch=batchid, par.prior=TRUE, prior.plots=FALSE)

  rm(expression_tmm_cpm)

  return(list(correctedExpression = correctedExpression,
              batchIDs = batchid,
              seriesIDs = series,
              parameters = list(archs4File = archs4File,
                                SamplesToConsider = SamplesToConsider,
                                seriesIDs = seriesIDs,
                                legacy = legacy
                               )
              )
         )
}





