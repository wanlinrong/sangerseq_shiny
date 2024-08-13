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
core_sanggerseq_split_new <- function(ab1,
                                  signal_noise_ratio = 0.33,
                                  ab1_trim5 = 50,
                                  ab1_trim3 = 100,
                                  ref_seq
                                  
){
  
  
  # ab1='../TSS20240409-023-00846全部结果文档/成功/AC-1_AF_TSS20240409-023-00846_F01.ab1'
  # signal_noise_ratio = 0.33
  # ab1_trim5 = 50
  # ab1_trim3 = 100
  # ref_seq = readDNAStringSet('../ap.txt') %>% toString()
  
    seq <- readsangerseq(ab1) %>%  makeBaseCalls(ratio = signal_noise_ratio ) %>% 
            setAllelePhase(ref_seq, trim5 = ab1_trim3, trim3 = ab1_trim3)
           
  
  # pwalign::pairwiseAlignment(ref_seq,secondarySeq(seq),gapExtension = 10) %>% writePairwiseAlignments(block.width = 100)
  # pwalign::pairwiseAlignment(primarySeq(seq),secondarySeq(seq),gapExtension = 10) %>% writePairwiseAlignments(block.width = 100)
  # pwalign::pairwiseAlignment(ref_seq,primarySeq(seq),gapExtension = 20) %>% writePairwiseAlignments(block.width = 100)
}
