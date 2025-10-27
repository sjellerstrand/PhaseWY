### Outputs the read distribution in heterozygotes
# - plot with the depth of the genotype against the count of the minor allele
# - histogram of the proportion of the minor allele reads in the total reads

## Set parameters
r=$RANDOM
f=$VCF_OUT.vcf.gz

## Creates individual file
if [ $NONDIPLOIDS == Yes ]; then
	cat $METADATA/ploidy.txt | awk -F'\t' '{if($2==2) print $1}' > $OUTDIR/$r".tmp.indv"
else
	bcftools view -h $f | tail -n1 | cut -f 10- | sed 's/\t/\n/g' > $OUTDIR/$r".tmp.indv"
fi;
samples=`seq 1 1 $(cat $OUTDIR/$r".tmp.indv" | wc -l)`;

## Writes the AD fields from all heterozygous genotypes with minor allele count i into a file hets.i
echo "Looping through individuals:"
n=$(wc -l $OUTDIR/$r".tmp.indv" | tr -s " " | cut -f 1 -d " ")
c=1
for j in $(seq 1 1 $(cat $OUTDIR/$r".tmp.indv" | wc -l)); do
i=$(cat $OUTDIR/$r".tmp.indv" | head -n $j | tail -n 1)
k=$(cat $OUTDIR/$r".tmp.indv" | head -n $j | tail -n 1)
echo -ne "$i | $c of $n"\\r
vcftools --gzvcf $f --mac 1 --max-mac 1 --indv $i --recode --stdout |\
grep -v "^#" | grep "0[/|]1" |\
awk '{split($9,a,":");for(i=1;i<=10;i++){if(match(a[i],"AD")){adidx=i}};for(i=10;i<=NF;i++){if(match($i,"0[/|]1")==1){split($i,b,":");print b[adidx]}}}' | \
awk -F',' '{if(NF == 2) print}' > $OUTDIR/$r".hets."$k
((c=c+1))
done

## Runs r-script to plot plots outlined above
echo -e "\nRscript outputting your plots"
Rscript -e 'args=commandArgs(TRUE);options(warn=-1);'\
'base<-args[1];s<-args[2];inds<-as.character(read.table(paste(s,"tmp.indv",sep="."))$V1);'\
'pdf(paste(base,"hetIndStats.pdf",sep="."),width=10,height=5);'\
'rcbpal<-c("#ffeda0","#fed976","#feb24c","#fd8d3c","#fc4e2a","#e31a1c","#bd0026","#800026");'\
'par(mfrow=c(1,2));rf<-colorRampPalette((rcbpal));r<-rf(12);'\
'for(i in inds){'\
'if(file.exists(paste(s,"hets",i,sep="."))){if(file.info(paste(s,"hets",i,sep="."))$size>0){'\
'd<-read.table(paste(s,"hets",i,sep="."),sep=",");'\
'd<-cbind(d,max=apply(d,1,max),min=apply(d,1,min));'\
'd<-cbind(d,minreadprop=d$min/(d$min+d$max),depth=d$min+d$max);'\
'x.bin<-seq(floor(min(d$depth,na.rm=T)),ceiling(max(d$depth,na.rm=T)),length=ceiling(max(d$depth,na.rm=T))-floor(min(d$depth,na.rm=T))+1);'\
'y.bin<-seq(floor(min(d$min,na.rm=T)),ceiling(max(d$min,na.rm=T)),length=ceiling(max(d$min,na.rm=T))-floor(min(d$min,na.rm=T))+1);'\
'freq<-as.data.frame(table(findInterval(d$depth,x.bin),findInterval(d$min,y.bin)));'\
'freq[,1] <- as.numeric(as.character(freq[,1]));'\
'freq[,2] <- as.numeric(as.character(freq[,2]));'\
'freq2D<-matrix(0,nrow=length(x.bin),ncol=length(y.bin));'\
'freq2D[cbind(freq[,1], freq[,2])] <- freq[,3];'\
'image(x.bin,y.bin,log10(freq2D),col=r,breaks = seq(0,max(log10(freq2D),na.rm=T),length.out=length(r)+1),'\
'main=paste(i,sep=""),xlab="DP",ylab="Reads minor alleles",xlim=c(0,max(x.bin,na.rm=T)),ylim=c(0,max(x.bin,na.rm=T)/2));'\
'barval<-round(10^(seq(0,max(log10(freq2D),na.rm=T),length.out = length(r)+1)));'\
'legend("topleft",legend=rev(barval[round(seq(2,length(barval),length.out = length(r)))]),'\
'fill=rev(r[round(seq(1,length(barval[-1]),length.out = length(r)))]),bty="n",ncol=3,title="heterozygous genotypes count");'\
'abline(0,1/2);abline(0,1/5,col="#00FF0088",lty=1);abline(0,3/10,col="#0000FF88",lty=1);'\
'legend("left",c("50:50 both alleles","70:30 common:rare","80:20"),'\
'lwd=1,lty=c(1,1,1),col=c(1,"#0000FF88","#00FF0088"),bty="n",cex=0.8);'\
'hist(d$minreadprop,breaks=seq(0,0.5,by=0.02),col=1,xlim=c(0,0.5),'\
'xlab="Proportion minor allele reads",main=paste(i,sep=""));'\
'abline(v=0.2,col="#00FF0088",lty=1);abline(v=0.3,col="#0000FF88",lty=1)'\
'}}'\
'};dev.off()' ${f%.vcf.gz} $OUTDIR/$r >/dev/null

## Removes temporary files
rm $OUTDIR/$r.hets.*;
rm $OUTDIR/$r.tmp.indv*;
