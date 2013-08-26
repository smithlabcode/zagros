vpath %.pwm ../../output

SORTBED=/home/rcf-47/bahramis/xpipe/bin/sortbed
EXTRACTSEQ=/home/rcf-47/bahramis/rmap/bin/extractseq
GENERATE_LOGO=generate_logo.R

ILLIDAN=/home/rcf-47/bahramis/cmb/illidan/bin/em

CHROM=/home/rcf-47/bahramis/panfs/cache/human/chromFaMasked/

EXTRACT_RADIUS=150

%_sorted.bed : %.bed
	$(SORTBED) -o $@ $^

%_sig_1000.bed : %.bed
	for i in {1,0.7,0.1,0.01,0.001,0.0001,0.00001,0.000001,0.0000001,0.00000001,0.000000001}; do export i; count=`awk -v val=$$i '{if ($$7 < val) print $$0;}' $^ | awk '{if (NR==1) {print $$0;c=$$1;s=$$2} if ($$1!=c || $$2-s>120) {print $$0;c=$$1;s=$$2;}}' | wc -l`; export count; if [ $$count -lt 1000 ];then echo $$count; awk -v val=$$i '{if ($$7 < val) print $$0;}' $^ | awk '{if (NR==1) {print $$0;c=$$1;s=$$2} if ($$1!=c || $$2-s>120) {print $$0;c=$$1;s=$$2;}}' > $@; break; fi; done;

%_expanded.bed : %_sig_1000.bed
	awk '{r = $(EXTRACT_RADIUS); mid = $$2+($$3-$$2)/2; if (mid-r > 0) print $$1"\t"int((mid-r))"\t"int((mid+r))"\t"$$4"\t"$$5"\t"$$6;}' $^ > $@

%_seqs.fa : %_expanded.bed
	$(EXTRACTSEQ) -o $@ -v -c $(CHROM) $^

%_filtered.fa : %_seqs.fa
	awk '{if (NR%2 == 1) a=$$0; else { found = 0;for (i=1;i<=length($$0);i+=1) if (substr($$0,i,1) == "N") found = 1; if ( found == 0) {print a;print $$0;}}}' $^ > $@

%.ild : %_filtered.fa
	$(ILLIDAN) -i $^ -l 7 -n 6 -v &> $@

%.pwm : %.ild
	awk '{if ($$0=="iteration 6...") i=NR; if (i!=0 && NR>=i+5 && NR<=i+13) {split($$0,a,","); for (j=1;j<length(a);j++) {if (NR<=i+8) pwm[NR-i-4,j]+=a[j]; if (NR>=i+10) pwm[NR-i-9,j]+=a[j];}}} END{ for (l=1; l< length(a);l++) {sum=0; for (k=1;k<5;k++) sum +=pwm[k,l]; for (k=1;k<5;k++) pwm[k,l]+=((1-sum)/4);} for (k=1;k<5;k++) {s="" pwm[k,1]; for (l=2; l< length(a);l++) {s=s "," pwm[k,l];} print s;}}' $^ > $@

%.pdf : %.pwm
	cat $(GENERATE_LOGO) | R --vanilla $^ ft

FILES=test_C6_em test_C6_gs

MOTIF=$(addsuffix .ild, $(FILES))
SIG=$(addsuffix _sig_1000.bed, $(FILES))
FILTERED=$(addsuffix _filtered.fa, $(FILES))
PWMS=$(addsuffix .pwm, $(FILES))
LOGOS=$(addsuffix .pdf, $(FILES))


significant : $(SIG)
filtered : $(FILTERED)
motif : $(MOTIF)
pwms : $(PWMS)
logos : $(LOGOS)
