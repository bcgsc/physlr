#install.packages("ggplot2")
#install.packages("optparse")
#install.packages("dplyr")
library("optparse")
library("ggplot2")
library("dplyr")
options <- list(make_option(c("-t", "--timefile"), help = "time file to plot"),
                make_option(c("-m", "--memoryfile"), help = "memory file to plot"),
                make_option(c("-n", "--minN"), type = "integer", help = " "),
                make_option(c("-N", "--maxN"), type = "integer", help = " "),
                make_option(c("-C", "--maxC"), type = "integer", help = " ") )

arguments <- parse_args(OptionParser(usage = "%prog [options] file file", 
                                     option_list=options,
                                    prog = "plot time-memory usages of Physlr"))
#time plot:
file_path<-sprintf("%s/%s", getwd(), arguments$timefile)
key<-strsplit(arguments$timefile, "[.]")[[1]][1]
pdf(paste0(key,".time.pdf"))
time_file<-read.csv(file_path, sep = "\t", header=T)
totaltime<-time_file %>% summarize(totaltime = sum(Elapsedtime)) %>% as.double() %>% round(digits=2)
time_file$Physlrstep<-factor(time_file$Physlrstep, levels=time_file$Physlrstep)
ggplot(time_file, aes(time_file$Elapsedtime, time_file$Physlrstep)) + 
  geom_segment(size=2, aes(x = time_file$Begin, y = time_file$Physlrstep, xend = time_file$End, yend = time_file$Physlrstep))+ 
  xlab("Time (min)") + ylab("Physlr step") + theme_bw() + theme(text=element_text(size=8))+ 
  ggtitle(paste("n", arguments$minN, "-", arguments$maxN, ", C", arguments$maxC, ": Total= ", totaltime, " min")) + 
  geom_text(size=2.5, aes(x=(time_file$Begin + time_file$End)/2, y=time_file$Physlrstep), hjust=0.5, vjust=-1, label=paste(round(time_file$Elapsedtime, digits = 2), " min", sep=""))
dev.off()

#memory plot:
gb=1000000
file_path<-sprintf("%s/%s", getwd(), arguments$memoryfile)
key<-strsplit(arguments$memoryfile, "[.]")[[1]][1]
pdf(paste0(key,".mem.pdf"))
memory_file<-read.csv(file_path, sep = "\t", header=T)
totalmem<-max((memory_file$Memory)/gb) %>% 
  as.double() %>% round(digits=2)
memory_file$Physlrstep<-factor(memory_file$Physlrstep, levels=memory_file$Physlrstep)
ggplot(memory_file, aes(memory_file$Memory/gb, memory_file$Physlrstep)) + 
  geom_segment(size=2, aes(x=0, y=memory_file$Physlrstep, xend = memory_file$Memory/gb, yend = memory_file$Physlrstep)) + 
  xlab("Memory (GB)") +ylab("Physlr step") + theme_bw()+theme(text=element_text(size=10))+
  ggtitle(paste("n", arguments$minN, "-", arguments$maxN, ", C", arguments$maxC, ": Peak= ", totalmem, " GB")) + 
  geom_text(size=2.5, aes(x=memory_file$Memory/gb, y=memory_file$Physlrstep), hjust=0.5, vjust=-1, label=paste(round(memory_file$Memory/gb, digits = 2), " GB", sep=""))
dev.off()

