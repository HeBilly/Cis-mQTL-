rm(list=ls() )
# option: Allow the user to set and examine a variety of global options which affect the way in which R computes and displays its results.
# width controls the maximum number of columns on a line used in printing vectors, matrices and arrays, and when filling by cat.
options (width = 200)
library(kinship2)
library(coxme)
#Provides access to a copy of the command line arguments supplied when this R session was invoked.
#This is especially useful with the --args command-line flag to R, as all of the command line after that flag is skipped.
#If trailingOnly = TRUE, a character vector of those arguments (if any) supplied after --args.
t.args<-  commandArgs(trailingOnly = T)
print(t.args)
#"cg05313009" "22"
cgid.input<-  t.args[ 1 ]
chr<-         as.integer( t.args[ 2 ] )

####
#cgid.input<-  "cg00000165"

#-----------------------------------------------------------------
#---	optimization by J Chung
#---	Find out the temporary directory whcih SCC automatically allocates for a job
#---	This temporary directory is located in each core (or node), which is not in the main hard disk.
#-----------------------------------------------------------------

# try is a wrapper to run an expression that might fail and allow the user's code to handle error-recovery.
# Sys.getenv obtains the values of the environment variables.
tmp.dir<-  try( Sys.getenv("TMPDIR"), silent=TRUE )
if ("try-error" %in% class(tmp.dir) | tmp.dir == "")     {
	print("...working not in queue")
#getwd returns an absolute filepath representing the current working directory of the R process; setwd(dir) is used to set the working directory to dir.
	tmp.dir<-  getwd()
}	else {
#################$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ why use <5
	if (nchar(tmp.dir) < 5) stop("TMPDIR variable is not set")
}
#---	This is the temporary directory for your job.
print( tmp.dir )
# "/scratch/5293182.1.kulisgpu-pub"

######read PED, cov and PC
fid_id<-  read.table("/restricted/projectnb/gaw20/user/jiayiwu/plink/fid_id.txt",as.is=T,header=T,sep=" ")
#head(fid_id[,1:4])
  #GAWSUBJ GPEDID GAWDAD GAWMOM
  #1       8    477      0      0
  #2      30    368   6611   7973
  #3      37    352   3360   3263
  #4      47    475    878   4467
  #5      56     71   5437   5640
  #6      66    276   4580    145
  
data.cov<-  read.csv("/restricted/projectnb/gaw20/user/jiayiwu/plink/covar.csv.gz", as.is=T)
data.cov<-  data.cov[,1:4]
 #head(data.cov[,1:4])*/
   #/*GAWSUBJ age center      X_smoking*/
   #/*1       6  40      0 a Never smoker*/
   #/*2       8  62      1 a Never smoker*/
   #/*3      30  51      1 a Never smoker*/
   #/*4      37  44      1 a Never smoker*/
   #/*5      47  43      1 a Never smoker*/
   #/*6      56  31      1 a Never smoker*/
   
fn.pcs<-  read.table("/restricted/projectnb/gaw20/PCA/GAW_filt_evecs.txt",as.is=T, header=T)
#head(fn.pcs[,1:4])
            #PC1          PC2          PC3          PC4
            #8   0.015071412 -0.030225789 -0.027096736  0.016521928
            #30  0.005488187  0.019368955  0.007138088  0.007779159
            #37  0.004884281  0.011810126  0.013812618  0.013220119
            #47  0.003949411  0.008420765  0.001646705 -0.019787052
            #56  0.001487723  0.023058138  0.001228932 -0.016201773
            #66 -0.020463283  0.059196821 -0.009762994  0.052187625

##-----------------------------------------------------------------
# read in sample id and merge the dosage with sample id
#-----------------------------------------------------------------
sample<-read.table("/restricted/projectnb/gaw20/user/jiayiwu/plink.geno/gaw20.geno.chr1.sample", as.is=T, sep=" ", header=T) 
colnames(sample)[1]<-"GAWSUBJ"
print(head(sample))
  #GAWSUBJ ID_2 missing sex phenotype
  #1       0    0       0   D         B
  #2       8    8       0   0      <NA>
  #3      30   30       0   0      <NA>
  #4      37   37       0   0      <NA>
  #5      47   47       0   0      <NA>
  #6      56   56       0   0      <NA>
sample<-sample[-1,]
# combine IDs, COV, PCs
fn.pcs$GAWSUBJ<-rownames(fn.pcs)
head(fn.pcs)
#            PC1          PC2          PC3          PC4           PC5          PC6          PC7          PC8         PC9         PC10         PC11         PC12         PC13         PC14         PC15
#8   0.015071412 -0.030225789 -0.027096736  0.016521928  0.0644492215 -0.002772859 -0.020180721  0.019085818  0.02242878 -0.028046120  0.018944507 -0.009270432 -0.047213932  0.022153637 -0.072016176
#30  0.005488187  0.019368955  0.007138088  0.007779159 -0.0517168807  0.049316665 -0.009530737 -0.039743794 -0.02242190  0.002723822 -0.017027176  0.024431899  0.028830545 -0.020449751 -0.015092092
#37  0.004884281  0.011810126  0.013812618  0.013220119  0.0022721005  0.020417614  0.009144295 -0.022299136 -0.02216175 -0.021436685  0.003591056  0.018373844 -0.009351711 -0.006717418 -0.008638328
#47  0.003949411  0.008420765  0.001646705 -0.019787052  0.0224542170 -0.044512239 -0.043897834  0.015456549 -0.11554012 -0.064968498 -0.037954119  0.004773336  0.016235295  0.097127182  0.005090790
#56  0.001487723  0.023058138  0.001228932 -0.016201773  0.0005112628  0.005689348  0.015768398 -0.009818731 -0.02567243 -0.009028544  0.020568159  0.016184612  0.015721692 -0.011135050 -0.001741268
#66 -0.020463283  0.059196821 -0.009762994  0.052187625 -0.0079905022 -0.102179577  0.206788720  0.170400380  0.20717735  0.157797097  0.008983123 -0.116320500 -0.069624987 -0.009318259 -0.093336984
#           PC16         PC17         PC18         PC19         PC20 GAWSUBJ
#8   0.036734379 -0.004411523 -0.039702516  0.017685992  0.030548340       8
#30 -0.005687367  0.049480289 -0.022138664  0.005546721 -0.003455965      30
#37 -0.015510976 -0.007776686  0.002442928 -0.007496588 -0.004507958      37
#47  0.006624046 -0.090393135 -0.023148785  0.154357336  0.028254224      47
#56  0.002902404 -0.010224235  0.004257404  0.029587731 -0.009055841      56
#66 -0.131948773 -0.018294647 -0.029141343  0.011581517 -0.068048985      66

id_cov<-merge(fid_id, data.cov, by="GAWSUBJ")

id_cov_pc<-merge(id_cov, fn.pcs, by="GAWSUBJ")
print(head(id_cov_pc)[,1:10])
#  GAWSUBJ GPEDID GAWDAD GAWMOM sex age center        X_smoking          PC1          PC2
#1       8    477      0      0   2  62      1   a Never smoker  0.015071412 -0.030225789
#2      30    368   6611   7973   1  51      1   a Never smoker  0.005488187  0.019368955
#3      37    352   3360   3263   2  44      1   a Never smoker  0.004884281  0.011810126
#4      47    475    878   4467   2  43      1   a Never smoker  0.003949411  0.008420765
#5      56     71   5437   5640   2  31      1   a Never smoker  0.001487723  0.023058138
#6      66    276   4580    145   2  38      0 c Current smoker -0.020463283  0.059196821



# load CpG probe info
load("/restricted/projectnb/gaw20/user/jiayiwu/rawdata/chr_bp_phenotyp.RData")
## extract sbp, ebp based on cgxxxx
methyl_rate_in_aChr<-list_methyl_rate_by_chr[[chr]]
 #head(list_methyl_rate_by_chr[[22]][,1:6])
          #snp            Genes       BP chrom GAWSUBJ_722 GAWSUBJ_8081
          #1 cg00014104            MAPK1 22220367    22  0.59206546 -0.164691418
          #2 cg00017461              OSM 30663316    22  0.01902528  0.018991190
          #3 cg00021762             <NA> 27832684    22  0.90863829 -0.056790003
          #4 cg00034755             <NA> 50918361    22  0.83235815  0.087117663
          #5 cg00047287 CPT1B;CHKB-CPT1B 51016899    22  0.97047732 -0.025508537
          #6 cg00047815           PHF21B 45311501    22  0.96932144 -0.001386367

# get the 1st line 
methyl_rate_in_aChr<- subset( methyl_rate_in_aChr, snp==cgid.input )

head(methyl_rate_in_aChr$Genes)
#1246 Levels: A4GNT AADAC AADACL2 ABCC5 ABCF3 ABHD10 ABHD14A ABHD14B ABHD14B;ABHD14A ABHD14B;PCBP4 ABHD5 ABHD6 ABI3BP ABTB1 ABTB1;PODXL2 ACAA1 ACAA1;MYD88 ACAD11 ACAD11;CCRL1 ACAD9 ACAP2 ACOX2 ... ZXDC
methyl_rate_in_aChr$Genes<- as.character(  methyl_rate_in_aChr$Genes  )
# head(methyl_rate_in_aChr$Genes)
#[1] "KREMEN1"

###get basepair column
head(methyl_rate_in_aChr[,1:6])
#           snp   Genes       BP chrom GAWSUBJ_722 GAWSUBJ_8081
#956 cg05313009 KREMEN1 29469033    22  0.05768789   0.02616524
methyl_rate_in_aChr_bp<-    as.integer( methyl_rate_in_aChr[,c("BP")] )
methyl_rate_in_aChr_gene<-  methyl_rate_in_aChr[,c("Genes")]
#print(length(methyl_rate_in_aChr_bp)) #[1] 1

# index each one and calculate 1Mb around it
sbp<-  methyl_rate_in_aChr_bp -1000000
ebp<-  methyl_rate_in_aChr_bp +1000000

print(sbp);print(ebp)
#[1] 90194674
#[1] 92194674

# extract sbp, ebp from dosage file
fn.dosage1<-  "/restricted/projectnb/gaw20/user/jiayiwu/plink.geno/gaw20.gen.chrXXX.impute2.dos1.gz"
fn.dosage1<-  gsub("XXX",chr,fn.dosage1)
#-----------------------------------------------------------------
#---	Define the output file location and names.
#-----------------------------------------------------------------
#dir.out<- 			"/restricted/projectnb/gaw20/user/jiayiwu/assoc.geno/chrXXXX"
#dir.out<-				gsub( "XXXX", chr, dir.out)
	dir.out<-				"/restricted/projectnb/gaw20/user/jiayiwu/assoc.geno/test"
	fn.out.in.dir<-  		paste( dir.out, paste("gaw20",chr,methyl_rate_in_aChr_bp,cgid.input,"assoc.out.gz", sep="."), sep="/" )	#--- It is a gzip
	fn.out.in.tmp.dir<-  	paste( tmp.dir, paste("gaw20",chr,methyl_rate_in_aChr_bp,cgid.input,"assoc.out", sep="."), sep="/" )	#--- It is not a gzip
	print(paste( "The final output file name will be ", fn.out.in.dir ) )
	print(paste( "The temporary output file name will be ", fn.out.in.tmp.dir ) )


#open an output file, and start the writing mode
	file.out<-  file(fn.out.in.tmp.dir,"w")

#---  write a header line into a file
	header.line<-   c("#Cg_chr","Cg_bp","Cg_id","Cg_gene","CHR","BP","LOC_ID","A1","A2","N","freq","info","beta","SE","p")
# diff b/t sep and collapse: https://gist.github.com/briandk/d9231ba1e2603eed0df1
	header.line<-   paste( header.line, collapse="\t" )
	cat( header.line, file=file.out, fill=256 )


#fill:a logical or (positive) numeric controlling how the output is broken into successive lines.
# If FALSE (default), only newlines created explicitly by "\n" are printed. Otherwise, the output is broken into 
#lines with print width equal to the option width if fill is TRUE, or the value of fill if this is numeric.
# Non-positive fill values are ignored, with a warning.

#------------------------------read dosage 100 lines at a time----------------------------------
buffer.nlines<- 100	     ### the number of lines to read dosage lines at a time
cmd.tabix<-	paste( "tabix ", fn.dosage1, paste(chr,":",sbp,"-",ebp,sep="") )
print( cmd.tabix )

#-----------------------------running analysis 100snps at a time--------------------------------
##### pipe: open connection and in reading mode
pipe.tabix<-    pipe( cmd.tabix, "r" )
print("...computing regression")
k<-	0
while( length( stream<- readLines( pipe.tabix, n=buffer.nlines) ) > 0 )    {

#'do.call' constructs and executes a function call from a name or a
#function and a list of arguments to be passed to it.
#take the lines and put it into dataframe
	stream.df<-     data.frame(do.call(rbind,strsplit(stream,"\t")),stringsAsFactors=F)

	k<-     k+1;	cat('.'); cat("--------------------------------------------------------",k, fill=T)
	print(dim(stream.df))
	print(head(stream.df[,1:10]))
#  X1       X2           X3 X4 X5 X6 X7 X8 X9 X10
#1 22 29977094 rs4820043XXX  A  G  0  0  1  1   1
#2 22 29982338 rs11089493XX  A  C  0  1  1  0   0
#3 22 29984823 rs7287402XXX  G  A  1  1  2  1   1
#4 22 29986265 rs5753527XXX  A  G  0  0  0  0   0
#5 22 29989101 rs739427XXXX  G  C  1  1  2  1   1
#6 22 29994496 rs16989459XX  C  G  0  1  1  0   0



#-----------------------------------------------------------------
#---	Extract the SNPs info from the stream.df  dosage data
#-----------------------------------------------------------------
# divide dosage file into 2 parts, transpose the 2nd part
	i.data.info<-   stream.df[, 1:5]; colnames( i.data.info )<- c('CHR','BP','SNP','A1','A2')
# add dummy vairable SNP names to i.data.info
	i.data.info$SNP_ID<-paste0("SNP",1:nrow(i.data.info))
	print("a")
	print(head(i.data.info))
#HR       BP          SNP A1 A2 SNP_ID
#1  22 28473206 rs140123XXXX  G  A   SNP1
#2  22 28477038 rs140126XXXX  A  G   SNP2
#3  22 28480549 rs140129XXXX  C  T   SNP3
#4  22 28485113 rs131278XXXX  A  G   SNP4
#5  22 28493146 rs16988015XX  T  C   SNP5
#6  22 28508305 rs140145XXXX  T  C   SNP6

#-----------------------------------------------------------------
#---	Extract the SNP dosages from the stream.df  dosage data
#-----------------------------------------------------------------
#---	What is i.data.tdos1
#------	what is row?
#------	what is column?
	i.data.tdos1<-  as.data.frame(t( stream.df[, 6:ncol(stream.df)] ),stringAsFactors=F)
	print(head(i.data.tdos1[,1:6]))
#.1    V1 V2 V3 V4 V5 V6
#X6   1  1  1  1  0  0
#X7   1  1  1  0  0  1
#X8   2  2  2  1  0  0
#X9   1  1  1 NA  0  1
#X10  1  1  1  0  0  0
#X11  1  1  1  1  0  1
#.2    V1 V2 V3 V4 V5 V6
#X6   1  0  1  1  0  0
#X7   2  0  1  1  0  0
#X8   0  0  0  0  2  0
#X9   1  0  1  1  0  0
#X10  0  1  0  0  0  1
#X11  1  0  1  1  1  0
#.3    V1 V2 V3 V4 V5 V6
#X6   1  0  0  0  0  1
#X7   1  1  0  0  0  1
#X8   0  0  0  0  0  0
#X9   1  0  0  0  0  0
#X10  1  0  0  0  0  1
#X11  0  1  0  0  0  0
#.4    V1 V2 V3 V4 V5 V6
#X6   0  0  1  0  1  0
#X7   0  1  1  0  1  1
#X8   1  1  2  0  2  1
#X9   1  0  1  0  1  0
#X10  1  0  1  0  1  0
#X11  0  2  2  0  2  1


#-----------------------------------------------------------------
# add dummy variable column "SNP_ID" to dosage file
#-----------------------------------------------------------------
	colnames(i.data.tdos1)<-paste0("SNP", 1:ncol(i.data.tdos1))
	print(head(i.data.tdos1[,1:6]))
#  SNP1 SNP2 SNP3 SNP4 SNP5 SNP6
#1    1    1    1    1    0    0
#2    1    1    1    0    0    1
#3    2    2    2    1    0    0
#4    1    1    1   NA    0    1
#5    1    1    1    0    0    0
#6    1    1    1    1    0    1


## $$$$$$$$$$$$$$$$$$$$$$the value naturally is integer, why need to change to numeric?
	i.data.tdos1<-  apply( i.data.tdos1, 2, as.numeric )
	i.data.tdos1<-  as.data.frame( i.data.tdos1, as.is=T, stringsAsFactors=F )

	i.data.tdos1$GAWSUBJ<-sample$GAWSUBJ
### add in the end the colnames with GAWSUBJ
	print(head(i.data.tdos1[,1:6]))

#index 1st cpg probe
	cgid_mehtyl0<-methyl_rate_in_aChr
	head(methyl_rate_in_aChr[,1:10])
#           snp   Genes       BP chrom GAWSUBJ_722 GAWSUBJ_8081 GAWSUBJ_6430 GAWSUBJ_7674 GAWSUBJ_5646 GAWSUBJ_5553
#956 cg05313009 KREMEN1 29469033    22  0.05768789   0.02616524   0.00357622    0.0463826    0.0131715    0.0450444

	cgid_methyl<-cgid_mehtyl0[,5:ncol(cgid_mehtyl0)]
	t_cgid_methyl<-as.data.frame(t(cgid_methyl), stringsAsFactors=F)
	colnames(t_cgid_methyl)[1]<-cgid.input
	t_cgid_methyl$GAWSUBJ<-gsub("GAWSUBJ_","",rownames(t_cgid_methyl))
	head(t_cgid_methyl)
#                cg05313009 GAWSUBJ
#GAWSUBJ_722  0.05768789     722
#GAWSUBJ_8081 0.02616524    8081
#GAWSUBJ_6430 0.00357622    6430
#GAWSUBJ_7674 0.04638260    7674
#GAWSUBJ_5646 0.01317150    5646
#GAWSUBJ_5553 0.04504440    5553

#-----------------------------------------------------------------
# merge the 1st cpg methylation value to id_cov_pc
#-----------------------------------------------------------------
	id_cov_pc_meth<-merge(id_cov_pc, t_cgid_methyl,by ="GAWSUBJ")
	head(id_cov_pc_meth)
#  GAWSUBJ GPEDID GAWDAD GAWMOM sex age center      X_smoking         PC1          PC2          PC3          PC4           PC5          PC6          PC7           PC8          PC9         PC10
#1       8    477      0      0   2  62      1 a Never smoker 0.015071412 -0.030225789 -0.027096736  0.016521928  0.0644492215 -0.002772859 -0.020180721  0.0190858180  0.022428778 -0.028046120
#2      30    368   6611   7973   1  51      1 a Never smoker 0.005488187  0.019368955  0.007138088  0.007779159 -0.0517168807  0.049316665 -0.009530737 -0.0397437942 -0.022421904  0.002723822
#3      47    475    878   4467   2  43      1 a Never smoker 0.003949411  0.008420765  0.001646705 -0.019787052  0.0224542170 -0.044512239 -0.043897834  0.0154565488 -0.115540123 -0.064968498
#4      56     71   5437   5640   2  31      1 a Never smoker 0.001487723  0.023058138  0.001228932 -0.016201773  0.0005112628  0.005689348  0.015768398 -0.0098187309 -0.025672431 -0.009028544
#5      82    478   4705   1604   1  43      1 a Never smoker 0.006240577 -0.005591322  0.012293364 -0.009777538  0.0103198916 -0.013509362  0.034505827 -0.0166344742 -0.004625279  0.003975996
#6      99     79   7122   5308   2  75      0 a Never smoker 0.031766515 -0.022118373  0.021640755  0.041103477 -0.0124621312  0.001407743  0.008804465  0.0009062255 -0.022088197 -0.024463233
#         PC11         PC12         PC13         PC14         PC15         PC16         PC17         PC18          PC19         PC20  cg05313009
#1  0.01894451 -0.009270432 -0.047213932  0.022153637 -0.072016176  0.036734379 -0.004411523 -0.039702516  0.0176859924  0.030548340  0.01976041
#2 -0.01702718  0.024431899  0.028830545 -0.020449751 -0.015092092 -0.005687367  0.049480289 -0.022138664  0.0055467211 -0.003455965  0.01462155
#3 -0.03795412  0.004773336  0.016235295  0.097127182  0.005090790  0.006624046 -0.090393135 -0.023148785  0.1543573360  0.028254224 -0.00759868
#4  0.02056816  0.016184612  0.015721692 -0.011135050 -0.001741268  0.002902404 -0.010224235  0.004257404  0.0295877312 -0.009055841  0.01235163
#5  0.00225225 -0.022992690  0.003723964 -0.041557795  0.003807874 -0.025029029 -0.008476069 -0.051755341 -0.0002234098  0.009461303  0.04735358
#6  0.03899414  0.024576372 -0.046200217  0.003689658  0.063043691 -0.055576011 -0.008124302 -0.149199337 -0.0700030007  0.062690524  0.01165548

# add 1 snp dosage to ped_cov_pheno file
#for (i in 1:col(i.data.tdos1))

#compute sample size
	sample.size<-nrow(na.omit(id_cov_pc_meth))

#calculate kinship matrix
	kmat<-makekinship(id_cov_pc_meth$GPEDID, id_cov_pc_meth$GAWSUBJ, id_cov_pc_meth$GAWMOM, id_cov_pc_meth$GAWDAD)

#open an output file, and start the writing mode
#	file.out<-  file(fn.out.in.tmp.dir,"w")

#for (i in 1:10){
	for (i in 1:nrow(i.data.info)){
  		i.marker<-			paste0("SNP",i)


  		snp_table<-			as.data.frame(i.data.tdos1[, ncol(i.data.tdos1)],stringsAsFactors=F)
  		snp_table$SNP_dummy<-	i.data.tdos1[,i]
  		colnames(snp_table)<-	c("GAWSUBJ",i.marker)
		head(snp_table)
#  GAWSUBJ SNP1
#1      6436    1
#2       722    1
#3      8081    2
#4      6430    1
#5      7674    1
#6      2933    1
#7      5646    2
#8      5553    0


  # merge calculated snp with ped, cov, cg probe
 		 id_cov_pc_meth_snp<-	merge(id_cov_pc_meth,snp_table, by ="GAWSUBJ")
 		head(id_cov_pc_meth_snp)
#  GAWSUBJ GPEDID GAWDAD GAWMOM sex age center      X_smoking         PC1          PC2          PC3          PC4           PC5          PC6          PC7           PC8          PC9         PC10
#1       8    477      0      0   2  62      1 a Never smoker 0.015071412 -0.030225789 -0.027096736  0.016521928  0.0644492215 -0.002772859 -0.020180721  0.0190858180  0.022428778 -0.028046120
#2      30    368   6611   7973   1  51      1 a Never smoker 0.005488187  0.019368955  0.007138088  0.007779159 -0.0517168807  0.049316665 -0.009530737 -0.0397437942 -0.022421904  0.002723822
#3      47    475    878   4467   2  43      1 a Never smoker 0.003949411  0.008420765  0.001646705 -0.019787052  0.0224542170 -0.044512239 -0.043897834  0.0154565488 -0.115540123 -0.064968498
#4      56     71   5437   5640   2  31      1 a Never smoker 0.001487723  0.023058138  0.001228932 -0.016201773  0.0005112628  0.005689348  0.015768398 -0.0098187309 -0.025672431 -0.009028544
#5      82    478   4705   1604   1  43      1 a Never smoker 0.006240577 -0.005591322  0.012293364 -0.009777538  0.0103198916 -0.013509362  0.034505827 -0.0166344742 -0.004625279  0.003975996
#6      99     79   7122   5308   2  75      0 a Never smoker 0.031766515 -0.022118373  0.021640755  0.041103477 -0.0124621312  0.001407743  0.008804465  0.0009062255 -0.022088197 -0.024463233
#         PC11         PC12         PC13         PC14         PC15         PC16         PC17         PC18          PC19         PC20  cg05313009 SNP1
#1  0.01894451 -0.009270432 -0.047213932  0.022153637 -0.072016176  0.036734379 -0.004411523 -0.039702516  0.0176859924  0.030548340  0.01976041    0
#2 -0.01702718  0.024431899  0.028830545 -0.020449751 -0.015092092 -0.005687367  0.049480289 -0.022138664  0.0055467211 -0.003455965  0.01462155    1
#3 -0.03795412  0.004773336  0.016235295  0.097127182  0.005090790  0.006624046 -0.090393135 -0.023148785  0.1543573360  0.028254224 -0.00759868    2
#4  0.02056816  0.016184612  0.015721692 -0.011135050 -0.001741268  0.002902404 -0.010224235  0.004257404  0.0295877312 -0.009055841  0.01235163    2
#5  0.00225225 -0.022992690  0.003723964 -0.041557795  0.003807874 -0.025029029 -0.008476069 -0.051755341 -0.0002234098  0.009461303  0.04735358    1
#6  0.03899414  0.024576372 -0.046200217  0.003689658  0.063043691 -0.055576011 -0.008124302 -0.149199337 -0.0700030007  0.062690524  0.01165548    1
#
#compute snp frequency
 		 i.freq<-	mean(id_cov_pc_meth_snp[, i.marker],na.rm=TRUE)/2


# compute imputation score of the snp
# 	if (is.na(i.freq) | i.freq==0 | i.freq==1){ 
 		if (is.na(i.freq) | i.freq<=0.05 | i.freq>=0.95){ 
			i.var <-"NA"
			}  else {
			i.var<- var ( i.data.tdos1[, c(i.marker)], na.rm=T ) / (2*i.freq*(1-i.freq))
		}

#print(head(id_cov_pc_meth_snp))
#for(i in sbp:ebp){
 		if (is.na(i.freq) | i.freq<=0.05 | i.freq>=0.95 | i.var<0.3){ 
			mod<-"NA"
			} else {
			mod<-lmekin(id_cov_pc_meth_snp[,cgid.input]~ id_cov_pc_meth_snp[,i.marker]+age + center +X_smoking + sex + PC1 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 
                  + (1|GAWSUBJ),data=id_cov_pc_meth_snp,varlist=list(kmat),na.action=na.omit)
		}
   # print(mod)
   # mod.out<-capture.output(mod)
   # r<-mod.out[11:11] 
#        if (is.na(i.freq) | i.freq==0 | i.freq==1){
 		 if (is.na(i.freq) | i.freq<=0.05 | i.freq>=0.95 | i.var<0.3){ 
      		r<-as.vector(c("NA", "NA", "NA", "NA"))
 			 } else {
     	 	r <- c(beta=mod$coef$fixed[2],se=sqrt(mod$var[2,2])) #regression coefficient and SE
      		rp<-pchisq((r[1]/r[2])^2,1,lower.tail=FALSE) #Wald test p-value
   #Given a list structure x, unlist simplifies it to produce a vector which contains all the atomic components which occur in x
	 		r<-as.vector(unlist(c(r,rp)))
  		}   ##returns regression coefficient, se, and Wald test p-value}

   #r <- c(beta=mod$coef$fixed[2],se=sqrt(mod$var[2,2])) #regression coefficient and SE
 		 i.outline <-c(c(chr,methyl_rate_in_aChr_bp,cgid.input,methyl_rate_in_aChr_gene), as.vector(unlist(i.data.info[i, 1:5])), sample.size,i.freq, i.var, r )
    #colnames(a)
    #result<-as.data.frame(c(result, a))
 		 i.outline<- paste( i.outline, collapse="\t" )
  		cat(i.outline, file=file.out, fill=500) #want to see what it look like
#rm(r);rm(mod);rm(rp);rm(a)
	}
#write.table(a, file="a.txt")

}
####fill=T or a number means after cat, it will go to a new line to print
#cat( cmd.tabix, fill=T )
#system( cmd.tabix )

#-----------------------------------------------------------------
#---	closs the access of the temporary file
#-----------------------------------------------------------------
close(pipe.tabix)
close(file.out)

#-----------------------------------------------------------------
#---	gzip the temporary output file
#-----------------------------------------------------------------
system( paste( "gzip -f ", fn.out.in.tmp.dir, sep="" ) )

   #r <- c(beta=mod$coef$fixed[2],se=sqrt(mod$var[2,2])) #regression coefficient and SE
    #rp<-pchisq((r[1]/r[2])^2,1,lower.tail=FALSE) #Wald test p-value
    #r<-c(r,rp)   ##returns regression coefficient, se, and Wald test p-value
    #c(length(r),r)
#}

#-----------------------------------------------------------------
#---	Move the temporary output to the final place
#-----------------------------------------------------------------
cmd.move<-	paste0( "mv ", fn.out.in.tmp.dir, ".gz", " ", fn.out.in.dir )
cat( cmd.move, fill=1000 )# want to see what it look like
system( cmd.move )



#-----------------------------------------------------------------
#---	Let's make a notice the R worked and ended without errors
#-----------------------------------------------------------------
cat("__END__in_R_without_errors\n")


