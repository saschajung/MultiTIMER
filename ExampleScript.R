#assumes repository cloned in the home of a mac/linux filesytem
#change ~ to your actual path if not, but you need MultiTIMER folder

#we get all scripts in R, which must be sourced
r_files <- list.files('~/MultiTIMER/R/', pattern = '\\.R$', full.names = TRUE)
#and all objects in the data folder, which must be loaded
rda_files <- list.files('~/MultiTIMER/data/', pattern = '\\.rda$', full.names = TRUE)
#objects in the data file must be loaded in the *global* environment, not general
#this way they're there for both R and the functions

sapply(rda_files, load, envir = .GlobalEnv)

# Source each R file
sapply(r_files, source)

#warning: if not downloaded, big, may take 90 minutes in fast connection
#if the userhome/Downloads folder doesn't exist, create it or replace path
f <- getArchS4Data(version = "v11",
                   org = "human",
                   path = "~/Downloads/")

samps_list <- selectArchS4Samples(archs4_file = f,
                                  scCutoff = 0.1,
                                  numAligned = 10^7,
                                  org = "human")

train_samps <- getTrainingSamples(archs4File = f,
                                  backgroundSamples = samps_list$samplesToConsider,
                                  samplesWithAge = samps_list$samplesWithAge,
                                  controlSamples = samps_list$controlSamples)

expr <- getSelectedSamplesArchS4(archs4File = f,
                                 SamplesToConsider = train_samps$trainingSamples,
                                 seriesIDs = samps_list$seriesIDs,
                                 legacy = T)

#example data
inp <- expr$correctedExpression
age <- unname(train_samps$trainingAge[colnames(inp)])
age <- transformAge(age,method = "identity")

model <- trainModel(trainExpr = inp,
                    trainAge = age)

#up to here must be ran as in the example script. in this line, import your
#input data as inp. or new variable name to not get confused
preds <- predictAge(model = model,
                    predExpr = inp)
