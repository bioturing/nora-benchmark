#running in paried end mode, $1 $2 are read files path, $3 output path
THREAD=32
mkdir $3

mkdir $3/star_rsem
\time -v -o $3/rsem/time.txt rsem-calculate-expression --no-bam-output --star --star-gzipped-read-file --paired-end -p $THREAD $1 $2 ./index/star/grch38 $3/star_rsem/result 2> $3/star_rsem/screen.log 

mkdir $3/bowtie2_rsem
\time -v -o $3/rsem/time.txt rsem-calculate-expression --no-bam-output  --bowtie2 --paired-end -p $THREAD $1 $2 ./index/rsem/grch38 $3/bowtie2_rsem/result 2> $3/bowtie2_rsem/screen.log 
 
mkdir $3/kallisto
\time -v -o  $3/kallisto/time.txt kallisto quant -i ./index/kallisto -o $3/kallisto -t $THREAD $1 $2 2> $3/kallisto/screen.log

mkdir $3/salmon
\time -v -o $3/salmon/time.txt salmon quant --index ./index/salmon --libType A -1 $1 -2 $2 -p $THREAD --output $3/salmon 2> $3/salmon/screen.log

mkdir $3/hera
\time -v -o $3/hera/time.txt hera quant -i ./index/hera/grch38 -1 $1 -2 $2 -t $THREAD -o $3/hera 2> $3/nora/screen.log

mkdir $3/nora
\time -v -o $3/nora/time.txt Nora quant -x ./index/nora/grch38 -1 $1 -2 $2 -t $THREAD -o $3/nora 2> $3/nora/screen.log
