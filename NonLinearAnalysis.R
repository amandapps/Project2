source("scaling_region_estimation.R")
library(RHRV)
library(dplyr)

records = c(16265,16273,16420,16483,16539,16773,16786,16795,17052,17453,18177,18184,19090,19093,19140)
methodEstimation = c("first.zero", "first.e.decay","first.minimum", "first.value")
porcentajes=c(0,50,100,200) 

dfmaestro = data.frame()
trace(estimateEmbeddingDim, edit=TRUE)
for (k in 1:length(records)) {
  hrv.data = CreateHRVData()
  hrv.data = SetVerbose(hrv.data, TRUE)
  filename = records[k]
  hrv.data = LoadBeatWFDB(
    hrv.data,
    RecordName = filename,
    RecordPath = "mit-bih-normal-sinus-rhythm-database-1.0.0/",
    annotator = "atr"
  )
 # hrv.data$Beat = hrv.data$Beat[1:500, , drop=FALSE ] #mas o menos 30 mins
  hrv.data=BuildNIHR(hrv.data)
  hrv.data = FilterNIHR(hrv.data)
  hrv.data = InterpolateNIHR (hrv.data, freqhr = 4)
  hrv.data = SetVerbose(hrv.data, TRUE)
  

  hrv.data=CreateNonLinearAnalysis(hrv.data)
  TimeLagEstimation = CalculateTimeLag(hrv.data, technique = "ami")
  EmbeddingDimEstimation = CalculateEmbeddingDim(hrv.data, numberPoints = 10000,
                                                 timeLag = TimeLagEstimation,
                                                 maxEmbeddingDim = 18)
  
  for (i in 1:length(porcentajes)) { 
  analysis_index = length(hrv.data$NonLinearAnalysis) 
  porcentajeVariacion=porcentajes[i]*TimeLagEstimation/100
  TimeLagEstimation2=as.integer(TimeLagEstimation+porcentajeVariacion) 
  
  #CorrDim
  hrv.datanoise=NonLinearNoiseReduction(hrv.data,
                                    embeddingDim = EmbeddingDimEstimation, 
                                    radius = NULL)
  hrv.data1=CalculateCorrDim(hrv.datanoise,indexNonLinearAnalysis = length(hrv.data$NonLinearAnalysis),
                            minEmbeddingDim = EmbeddingDimEstimation - 1,
                            maxEmbeddingDim = EmbeddingDimEstimation + 2,
                            timeLag = TimeLagEstimation2,minRadius = 1,
                            maxRadius = 100, pointsRadius = 100,
                            theilerWindow = 20, doPlot = TRUE)
  

  PlotCorrDim(hrv.data1)
  
  cd=hrv.data1$NonLinearAnalysis[[analysis_index]]$correlation$computations
  filteredCd= nltsFilter(cd, threshold=0.99)
  cdScalingRegion = estimate_scaling_region(filteredCd, numberOfLinearRegions = 4, doPlot = TRUE)
  hrv.data1=EstimateCorrDim(hrv.data1,regressionRange=cdScalingRegion) 
  
  
  #SampleEntropy
  se = sampleEntropy(cd, do.plot = F)
  plot(se, type = "o", add.legend = FALSE, log = "x")
  abline(v = cdScalingRegion, lwd = 3, lty = 2)
  hrv.data1=CalculateSampleEntropy(hrv.data1,indexNonLinearAnalysis = length(hrv.data$NonLinearAnalysis),doPlot=FALSE)
  hrv.data1=EstimateSampleEntropy(hrv.data1,regressionRange = cdScalingRegion) 
  
  #Lyapunov
  hrv.data2=CalculateMaxLyapunov(hrv.data,indexNonLinearAnalysis = length(hrv.data$NonLinearAnalysis), 
                                 minEmbeddingDim = EmbeddingDimEstimation,
                                 maxEmbeddingDim = EmbeddingDimEstimation+2,
                                 timeLag = TimeLagEstimation2, 
                                 radius = 40, theilerWindow = 50, doPlot = TRUE)
  # regression ----------------------------------------------------------------
  sampling.period = diff(hrv.data$Beat$Time)[1]
  indx = which.min(abs(colMeans(cd$corr.matrix) - 1e-3))
  useRadius = cd$radius[indx]
  ml=hrv.data2$NonLinearAnalysis[[analysis_index]]$lyapunov$computations
  plot(ml, type = "l")
  hrv.data2=EstimateMaxLyapunov(hrv.data2, indexNonLinearAnalysis = length(hrv.data$NonLinearAnalysis), 
                                useEmbeddings = (EmbeddingDimEstimation):(EmbeddingDimEstimation+2),
                                doPlot=TRUE) 

  #RQA (no regression!!)
  indx = which.min(abs(colMeans(cd$corr.matrix) - 1e-3))
  useRadius = cd$radius[indx]
  hrv.data3=RQA(hrv.data, indexNonLinearAnalysis = length(hrv.data$NonLinearAnalysis), embeddingDim = EmbeddingDimEstimation,
                timeLag = TimeLagEstimation2, radius=useRadius, doPlot = TRUE)
  
  #analysis_index = length(hrv.data$NonLinearAnalysis)
  resultados=hrv.data$NonLinearAnalysis[[analysis_index]]
  resultados$file=filename
  resultados$TimeLag=TimeLagEstimation
  resultados$PorcentajeVariacion=porcentajes[i]
  resultados$TimeLagUsado=TimeLagEstimation2
  resultados$EDim=EmbeddingDimEstimation
  resultados$CorrelationDimension=hrv.data1$NonLinearAnalysis[[analysis_index]]$correlation$statistic
  resultados$SampleEntropy=mean(hrv.data1$NonLinearAnalysis[[analysis_index]]$sampleEntropy$statistic) #mean pq sino devuelve 3 valores
  resultados$MaxExponentLyapunov= hrv.data2$NonLinearAnalysis[[analysis_index]]$lyapunov$statistic
  resultados$RQA_LAM=hrv.data3$NonLinearAnalysis[[analysis_index]]$rqa$LAM
  resultados$RQA_DET=hrv.data3$NonLinearAnalysis[[analysis_index]]$rqa$DET
  resultados$RQA_LMAX=hrv.data3$NonLinearAnalysis[[analysis_index]]$rqa$Lmax
  #resultados$RQA_LMEAN=hrv.data3$NonLinearAnalysis[[analysis_index]]$rqa$Lmean
  resultados$RQA_ENTROPY=hrv.data3$NonLinearAnalysis[[analysis_index]]$rqa$ENTR
  resultados$RQA_VMAX=hrv.data3$NonLinearAnalysis[[analysis_index]]$rqa$Vmax
  resultados$RQA_VMEAN=hrv.data3$NonLinearAnalysis[[analysis_index]]$rqa$Vmean
  
  
  df = as.data.frame(resultados)
  dfmaestro = rbind(dfmaestro, df)
  }
  
}


