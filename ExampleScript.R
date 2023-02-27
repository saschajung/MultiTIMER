library(MultiTIMER)
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

inp <- expr$correctedExpression
age <- unname(train_samps$trainingAge[colnames(inp)])
age <- transformAge(age,method = "identity")

model <- trainModel(trainExpr = inp,
                    trainAge = age)

preds <- predictAge(model = model,
                    predExpr = inp)





