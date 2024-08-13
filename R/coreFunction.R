
library(sangerseqR)
library(Biostrings)
library(tidyverse)

#' Title
#'
#' @param ab1 ab1 file
#' @param signal_noise_ratio defalt[noise/signal = 0.33]
#' @param ab1_trim5 defalt[50]
#' @param ab1_trim3 defalt[100]
#' @param ref_seq class[strings vectior]
#'
#' @return pairwiseAlignment
#' @export 
#'
#' @examples ......
core_sanggerseq_split <- function(ab1,
                                  signal_noise_ratio = 0.33,
                                  ab1_trim5 = 50,
                                  ab1_trim3 = 20,
                                  ref_seq
                                  
){
  
  
  #message('-------------------------------------------------')
  #message(ab1)
  

  # ab1='../TSS20240409-023-00846全部结果文档/成功/AC-1_AF_TSS20240409-023-00846_F01.ab1'
  # signal_noise_ratio = 0.33
  # ab1_trim5 = 50
  # ab1_trim3 = 20
  # ref_in = readDNAStringSet('../ap.txt') %>% toString()


hetsangerseq <- readsangerseq(ab1)
hetcalls <- makeBaseCalls(hetsangerseq, signal_noise_ratio) 
f <- primarySeq(hetcalls)
s <- secondarySeq(hetcalls)


#####确定方向
plus <- pwalign::pairwiseAlignment(ref_seq,f,gapExtension = 6.66,type = 'local')
munius <- pwalign::pairwiseAlignment(ref_seq %>% DNAString() %>% reverseComplement(),f,gapExtension = 6.66,type = 'local')
direction <- if_else(plus@score > munius@score,'+','-')
shinytoastr::toastr_info(message = paste0('direction is >',direction) ,position = "bottom-right",showDuration =5000)

#####转换测序结果，以ref方向为准

if(direction != '+') {f <- reverseComplement(f);s <- reverseComplement(s)}


#############
#单峰

if(f[ab1_trim5:(length(f)-ab1_trim3)]  != s[ab1_trim5:(length(f)-ab1_trim3)]) 
{
#############
seed <- map_dfr(ab1_trim5:(length(f)-ab1_trim3),function(i){
  
  seed_raw <- subseq(f,start = i,width = 10) %>% toString()
  
  loc <- str_locate(ref_seq ,seed_raw)
  
  if(!nrow(loc) <= 1) loc[1,'start'] <- NA 
  
  tibble(seed = seed_raw ,fstart = i,ref = loc[1,'start'])
  
})

info <- seed %>% filter(!is.na(ref)) %>% slice_min(order_by = fstart)


anyq <- map_chr(info$fstart:length(f),function(i){
  
  #i=50
 
 
  
 fq <- f[i]
 sq <- s[i]
 rq <- DNAString(ref_seq)[i + info$ref - info$fstart]
 
 if(rq == fq ) return(sq %>% toString()) else return(fq %>% toString())
 
}) %>% paste0(collapse = '')

fanyqr <- paste0(DNAString(ref_seq)[ (info$ref - info$fstart -ab1_trim5 + 1):(info$ref - info$fstart)],
                 anyq,
                 DNAString(ref_seq)[ (length(f) + info$ref - info$fstart+1) : (length(f) + info$ref - info$fstart + ab1_trim3) ]
)

} else 
{ 
  fanyqr <-  f %>% toString()
  
 }

return(fanyqr)

#pwalign::pairwiseAlignment(ref_seq,fanyqr,gapExtension = 6.66) %>% writePairwiseAlignments(block.width = 100) 

}
