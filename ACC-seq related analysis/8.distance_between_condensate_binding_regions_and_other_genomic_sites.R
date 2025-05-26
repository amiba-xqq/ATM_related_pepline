# The distance between condensate binding peaks or non-DSB pNBS1 peaks and their nearest ATM-activated DSB sites was analyzed.
bedtools closest -a pOHT_Ushape.bed -b DSB.bed -d > pOHT_Ushape_distance_to_total_DSB.bed
bedtools closest -a DMSO_Ushape.bed -b DSB.bed -d > DMSO_Ushape_distance_to_total_DSB.bed
bedtools closest -a pOHTIRAK1i_Ushape.bed -b DSB.bed -d > pOHTIRAK1i_Ushape_distance_to_total_DSB.bed

# The relative distance between condensate binding peaks and the nearest non-DSB pNBS1 peaks was assessed. 
bedtools reldist -a pOHT_Ushape.bed -b pNBS1_noDSB_peaks.bed > pOHT_Ushape_distance_to_pNBS1_noDSB_peak_reldist.txt
bedtools reldist -a DMSO_Ushape.bed -b pNBS1_noDSB_peaks.bed > DMSO_Ushape_distance_to_pNBS1_noDSB_peak_reldist.txt
bedtools reldist -a pOHTIRAK1i_Ushape.bed -b pNBS1_noDSB_peaks.bed > pOHTIRAK1i_Ushape_distance_to_pNBS1_noDSB_peak_reldist.txt
