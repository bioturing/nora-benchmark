#running in paried end mode, $1 $2 output path, $3 mean fragment length, $4 fragment length standard deviation
THREAD=32
mkdir $2

mkdir $2/star_rsem
\time -v -o $2/star_rsem/time.txt rsem-calculate-expression --fragment-length-mean $3 --fragment-length-sd $4 --no-bam-output --star -p $THREAD $1 ./index/star_rsem/grch38 $2/star_rsem/result 2> $2/star_rsem/screen.log 

mkdir $2/bowtie2_rsem
\time -v -o $2/bowtie2_rsem/time.txt rsem-calculate-expression --fragment-length-mean $3 --fragment-length-sd $4 --no-bam-output  --bowtie2 -p $THREAD $1 ./index/bowtie2_rsem/grch38 $2/bowtie2_rsem/result 2> $2/bowtie2_rsem/screen.log 
 
mkdir $2/kallisto
\time -v -o  $2/kallisto/time.txt kallisto quant -i ./index/kallisto -o $2/kallisto -t $THREAD -l $3 -s $4 --single $1 2> $2/kallisto/screen.log

mkdir $2/salmon
\time -v -o $2/salmon/time.txt salmon quant --index ./index/salmon --libType A --fldMean $3 --fldSD $4 -r $1 -p $THREAD --output $2/salmon 2> $2/salmon/screen.log

mkdir $2/hera
\time -v -o $2/hera/time.txt hera quant -i ./index/hera/grch38 -1 $1 -t $THREAD -o $2/hera 2> $2/nora/screen.log

mkdir $2/nora
\time -v -o $2/nora/time.txt Nora quant -x ./index/nora/grch38 -1 $1 --fmu $3 --fsigma $4  -t $THREAD -o $2/nora 2> $2/nora/screen.log
