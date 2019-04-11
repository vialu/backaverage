library(scales)

TS00033_posture <- read.table("~/Documents/NIH/Proyects/time_frequency/TS00033_posture.txt", quote="\"", comment.char="")
data   <-TS00033_posture[,4]

wind  <- 100
onset <- 50

time  <-seq(onset*-1,onset)

plot(data,type = "l")
emg     <- rescale(data)
plot(emg,type = "l")

emg1    <- emg > 0.5
emgtrl  <- diff(emg1)
emgon   <- which(emgtrl %in%  1)  # beggening of the burst
emgoff  <- which(emgtrl %in% -1) #finish of the burst
ons     <- which(diff(emgon)>20) # for real data has to be >40

on_after  <-0;
on_before <-0;

for (li in 1:length(ons)){
    if ((ons[li]-20)<1) next
  on_after[li]  <- mean(emg[ons[li]:(ons[li]+20)]) # should adjust according to sampling rate?
  on_before[li] <- mean(emg[(ons[li]-20):ons[li]])
}

ons2 <- ons[which(on_after>0.15 & on_before <0.1 )] # original numbers are 0.04 and 0.03
beg    <- ons2 -onset
ending <- ons2 + (wind-onset)
trial <- data.frame( beg[which(beg>0)], ending [which(beg>0)])


# Segmenting the data

segDat <-matrix(nrow = dim(trial)[1], ncol = wind+1)
for (ti in 1:dim(trial)[1]){
    segDat[ti,]<- data[trial[ti,1]:trial[ti,2]]
  }

aver <-apply(segDat,2,mean)
plot(time , aver,type = "l")
