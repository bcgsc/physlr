install.packages("ggplot2")
install.packages("optparse")
install.packages("dplyr")
install.packages("tidyverse")
library("optparse")
library("ggplot2")
library("dplyr")
library("tidyverse")
library("lubridate")
library("stringr")

#read time and memory files, convert time into start-interval-end format for plotting, and plot time and memory benchmarks!

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
time_orig <- read_tsv(file_path, col_types = "cc") #read the step name and time as characters
time_df <- time_orig %>%
  mutate(
    Minutes = as.duration(
      ifelse(str_count(time, ":") == 2, as.numeric(hms(time, quiet=TRUE))/60, as.numeric(ms(time, quiet=TRUE))/60)),
    stop = cumsum(Minutes),
    start = stop - Minutes)
time_df$Minutes<-as.numeric(time_df$Minutes)
time_df$stop<-as.numeric(time_df$stop)
time_df$start<-as.numeric(time_df$start)
pdf(paste0(key,".time.pdf"))
totaltime<-time_df %>% summarize(totaltime = sum(Minutes)) %>% as.double() %>% round(digits=2)
time_df$step<-factor(time_df$step, levels=time_df$step)
p<-ggplot(time_df, aes(time_df$Minutes, time_df$step)) +
  geom_segment(size=2, aes(x = time_df$start, y = time_df$step, xend = time_df$stop, yend = time_df$step))+
  xlab("Time (min)") + ylab("Physlr step") + theme_bw() + theme(text=element_text(size=10))+
  ggtitle(paste("n", arguments$minN, "-", arguments$maxN, ", C", arguments$maxC, ": Total= ", totaltime, " min")) +
  geom_text(size=2.5, aes(x=(time_df$start + time_df$stop)/2, y=time_df$step), hjust=0.5, vjust=-1, label=paste(round(time_df$Minutes, digits = 2), " min", sep=""))
dev.off()

#memory plot:
gb=1000000
file_path<-sprintf("%s/%s", getwd(), arguments$memoryfile)
key<-strsplit(arguments$memoryfile, "[.]")[[1]][1]
memory_file <- read_tsv(file_path)
pdf(paste0(key,".mem.pdf"))
maxmem<-max((memory_file$memory)/gb) %>%
  as.double() %>% round(digits=2)
memory_file$step<-factor(memory_file$step, levels=memory_file$step)
ggplot(memory_file, aes(memory_file$memory/gb, memory_file$step)) +
  geom_segment(size=2, aes(x=0, y=memory_file$step, xend = memory_file$memory/gb, yend = memory_file$step)) +
  xlab("Memory (GB)") +ylab("Physlr step") + theme_bw()+theme(text=element_text(size=10))+
  ggtitle(paste("n", arguments$minN, "-", arguments$maxN, ", C", arguments$maxC, ": Peak= ", maxmem, " GB")) +
  geom_text(size=2.5, aes(x=memory_file$memory/gb, y=memory_file$step), hjust=0.5, vjust=-1, label=paste(round(memory_file$memory/gb, digits = 2), " GB", sep=""))
dev.off()
