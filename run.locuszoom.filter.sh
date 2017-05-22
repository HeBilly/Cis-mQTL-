#!/bin/bash -l

echo "=========================================================="
echo "Starting on : $(date)"
echo "Running on node : $(hostname)"
echo "Current directory : $(pwd)"
echo "Current job ID : $JOB_ID"
echo "Current job name : $JOB_NAME"
echo "Task index number : $SGE_TASK_ID"
echo "=========================================================="
start_date="$(date)"
echo
dir_gaw20="/restricted/projectnb/gaw20/user/jiayiwu/assoc.geno"
echo "sh $dir_gaw20/run.locuszoom.filter.sh $*"

fn_assoc0=$1
chr=$2
sbp=$3
ebp=$4
ind_p=$5
ind_freq=$6
ind_info=$7
gname=$8
ethnic=$9

if [ ! -f $fn_assoc0 ]; then ls $fn_assoc0; exit; fi

temp_fn=`basename $fn_assoc0`
temp_nf=`echo $temp_fn | sed 's/\./ /g' | awk '{print NF}'`; temp_nf=$((temp_nf-1))
fn_out=`echo $temp_fn | cut -d"." -f-${temp_nf}`
fn_out="$fn_out.$gname.$chr.$sbp.$ebp.$ethnic.txt"
prefix_out="`echo $fn_out | sed 's/.txt$//g'`"

fn_assoc=`readlink -e $fn_assoc0 | awk '{print $NF}'`
#########################################################
####	1. make a temporary directory and move there
#########################################################
dir_pwd=`pwd`
dir_temp="$dir_pwd/temp_${prefix_out}"
rm -rf $dir_temp; mkdir $dir_temp; cd $dir_temp


#########################################################
####	2. extract SNPs within the region given by the user
#########################################################
echo
echo "_____________________________________________________________________"
echo "CHR $chr : $sbp - $ebp"
echo "GENE $gname"
echo "PREFIX $prefix_out"
echo "_____________________________________________________________________"
echo
echo "... extracting [chr bp p-value]"
	cmd="tabix -h $fn_assoc ${chr}:${sbp}-${ebp} | sed 's/#//g' | awk -v ind_freq=$ind_freq -v ind_info=$ind_info 'NR==1 || (\$ind_freq>=0.01 && \$ind_freq<=0.99 && \$ind_info>=0.4) {print \$0}' | grep -v NA > $fn_out.0"
#	cmd="tabix -h $fn_assoc ${chr}:${sbp}-${ebp} | sed 's/#//g' | awk -v ind_freq=$ind_freq -v ind_info=$ind_info 'NR==1 || (\$ind_freq>=0.05 && \$ind_freq<=0.95 && \$ind_info>=0.4) {print \$0}' | grep -v NA > $fn_out.0"
	echo $cmd
	eval $cmd

	cat $fn_out.0 | awk -v ind_p=$ind_p '{print $1,$6,$ind_p}' | grep -v NA | sort -gk3 > $fn_out
	rm -rf $fn_out.0
        ls $fn_out

#########################################################
####	3. annotatie the SNPs
#########################################################
echo "... annotating: `cat $fn_out | wc -l` snps"
echo

fn_in="${fn_out}"
header="T"

fn_temp_header="${fn_out}.header"
markername="MarkerName"
pvalue="P.value"

index_chr=0
index_bp=1
ind_p=4

echo "MarkerName CHR BP P.value WITHIN UP DOWN" > $fn_temp_header

cmd="sh /restricted/projectnb/adgc/_compute.gwas.nb/script.to.compute/annot/run.do.annotate.chr.bp.sh $fn_in $index_chr $index_bp $header"
echo $cmd
eval $cmd

cat $fn_in.annot | sed 's/ /\t/g' | cut -f-3,6- > $fn_in.annot.0
cat $fn_temp_header $fn_in.annot.0 | sed 's/ /\t/g' > $fn_in.annot

rm -rf  $fn_in.annot.0
rm -rf $fn_temp_header
rm -rf ${fn_out}


#########################################################
####	4. locuszoom
#########################################################
echo "... locuszoom"
echo $fn_in.annot

#cmd="/usr3/graduate/jychung/bin/locuszoom --meta $fn_in.annot --chr $chr --start $sbp --end $ebp --pop $ethnic --build hg19 --source 1000G_March2012 --prefix $prefix_out  --plotonly --pvalcol P.value --markercol MarkerName --verbose --delim tab geneFontSize=.8 smallDot=.3 largeDot=.9 format=pdf ymax=10 legend=none metalRug='Plotted SNPs' rfrows=4 weightCol=Weight showRecomb=TRUE warnMissingGenes=T showAnnot=FALSE --snpset=NULL"
cmd="/restricted/projectnb/adgc/_compute.gwas.nb/library/locuszoom/locuszoom/bin/locuszoom --meta $fn_in.annot --chr $chr --start $sbp --end $ebp --pop $ethnic --build hg19 --source 1000G_March2012 --prefix $prefix_out  --plotonly --pvalcol P.value --markercol MarkerName --verbose --delim tab --snpset=NULL --cache None geneFontSize=.8 smallDot=.3 largeDot=.9 format=pdf ymax=10 showRecomb=TRUE warnMissingGenes=T showAnnot=FALSE rugColor=white "

echo $cmd
eval $cmd || exit


gzip -9f $fn_in.annot

echo "... converting pdf to jpg"
fn_pdf_date=`date +"%y%m%d"`
fn_pdf="${prefix_out}_${fn_pdf_date}_chr${chr}_${sbp}-${ebp}"
echo
ls $fn_pdf.pdf || exit
cmd="convert -density 300 $fn_pdf.pdf[0] $fn_pdf.jpg"
echo $cmd
eval $cmd || exit

ls $fn_pdf.jpg || exit
rm -rf $fn_pdf.pdf
#rm -rf $fn_in.annot.gz

mv $fn_pdf.jpg $dir_pwd || exit
cd $dir_pwd
rm -rf  $dir_temp
echo
echo "...Successful!"
echo "=========================================================="
echo "Started  on : $start_date"
echo "Finished on : $(date)"
echo "=========================================================="
