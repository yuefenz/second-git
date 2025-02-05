####################################################################
#              基因学苑VIP课程（第4季）:基因组分析                   #
#                         王通 2024/11                             #
####################################################################

#####################################################
#              八十八、三代高质量测序数据拼接        #
#####################################################
#=========================
#       pacbio hifi      #
#=========================
#https://downloads.pacbcloud.com/public/dataset/2021-11-Microbial-96plex/

#1 使用smrtlink转换数据
#激活环境
mamba create -n smrtlink -y
mamba activate smrtlink
mamba install -c hcc -y smrtlink-tools

#建立索引
pbindex bc2080.bam
bam2fastq -o pacbio bc2080.bam 

#使用samtools转换数据
samtools fastq bc2008.bam > bc2008.fastq
samtools fastq bc2032.bam > bc2032.fastq
samtools fastq bc2080.bam > bc2080.fastq
samtools fastq bc2056.bam > bc2056.fastq

#压缩
pigz -p 20 *.fastq

#fastqc质控
mamba activate genome
mkdir qc
fastqc -f fastq -o qc -t 12 *.fastq.gz
multiqc -d qc -o multiqc

#拼接基因组
tmux new-session -s flye
time ~/mambaforge/bin/flye --pacbio-hifi bc2008.fastq.gz --out-dir bc2008 --genome-size 4.6m --threads 16


#=========================
#       nanopore q20+   #
#=========================
mamba activate genome
fasterq-dump ERR8958564.sra 
#fastqc质控
mkdir qc
fastqc -f fastq -o qc -t 12 ERR8958564.fastq 

#提取长度大于5K序列
seqkit seq -m 15000 ERR8958564.fastq > nanopore.fastq

#nanopore Q20+序列拼接
time ~/mambaforge/bin/flye --nano-hq nanopore.fastq --out-dir flye --genome-size 4.6m --threads 16

#####################################################
#                     五、基因预测              #
#####################################################
#=========================
#       安装分析软件     #
#=========================
mamba install -y blast
mamba install -y blat
mamba install -y lastz
mamba install -y muscle
mamba install -y clustalw
mamba install -y bwa
mamba install -y bwa-mem2
mamba install -y samtools
mamba install -y bowtie2
mamba install -y prodigal
mamba install -y glimmer3
mamba install -y augustus
mamba install -y trnascan-se
mamba install -y trf
mamba install -y repeatmasker

#=========================
#     原核生物基因预测    #
#=========================
mkdir gene

#1 prodigal基因预测
prodigal -a MGH78578.pep -d MGH78578.cds -f gff -g 11 -o MGH78578.gff -p single -s MGH78578.stat -i MGH78578.fasta  

#2 glimmer3基因预测
#运行软件
sed -e '/>/d' MGH78578.fasta |tr -d '\n' |awk 'BEGIN {print ">wholefile"}{print $0}' >wholefile
long-orfs -n -t 1.15 wholefile tagname.longorfs  1>/dev/null 2>/dev/null
extract -t wholefile tagname.longorfs > tagname.train  2>/dev/null
build-icm  -r tagname.icm < tagname.train 1>/dev/null 2>/dev/null
glimmer3 -o50 -g110 -t30 MGH78578.fasta tagname.icm ref

#3 同源比对预测
#下载参考序列基因集
https://www.ncbi.nlm.nih.gov/genome/?term=NC_009648
wget -c https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/240/185/GCF_000240185.1_ASM24018v2/GCF_000240185.1_ASM24018v2_protein.faa.gz
#建立索引
makeblastdb -in GCF_000240185.1_ASM24018v2_protein.faa -dbtype prot -parse_seqids -out GCF_000240185.1_ASM24018v2_protein.faa
#blastx比对
blastx -query MGH78578.fasta -out blast.out -db GCF_000240185.1_ASM24018v2_protein.faa -outfmt 6 -evalue 1e-5
#提取比对区域，生成bed文件
awk '{if ($7 < $8) print $1"\t"$7-1"\t"$8;else print $1"\t"$8-1"\t"$7}' blast.out >gene.bed
#根据比对位点，提取序列
seqkit subseq --bed gene.bed MGH78578.fasta  >MGH78578_gene.ffn

#=========================
#     真核生物基因预测    #
#=========================
#安装augustus软件
mamba create -n augustus -y augustus=3.4.0
#激活环境
conda activate augustus
#查看软件自带模型
augustus --species=help
augustus --strand=both --genemodel=partial --singlestrand=false --protein=on --introns=on --start=on --stop=on --cds=on --codingseq=on --alternatives-from-evidence=true --gff3=on --UTR=on --outfile=out.gff --species=arabidopsis ninanjie.fa

#####################################################
#                     六、基因功能注释             #
#####################################################
#安装eggnog-emapper
mamba install -c bioconda -y eggnog-mapper python=2.7

#下载数据库，比较慢
#下载地址：http://eggnog5.embl.de/download/emapperdb-5.0.2/
#http://eggnog5.embl.de/download/eggnog_5.0/

mkdir eggnog_database;
download_eggnog_data.py -y --data_dir eggnog_database

#镜像数据下载
ftp://download.nmdc.cn/tools/eggnog/eggnog.db.gz
ftp://download.nmdc.cn/tools/eggnog/eggnog_proteins.dmnd.gz

#基因功能注释
emapper.py -i MGH78578.pep --output annotation -m diamond --data_dir eggnog_database


#####################################################
#                   七、ncRNA分析            #
#####################################################
#1 核糖体RNA（rrna）预测
#安装软件 rnammer
#首先需要下载安装hmmer
mamba install -y hmmer=2.3.2 
mamba install -y perl-xml-simple
mamba install -y perl-getopt-long

#rnammer 需使用教育edu邮箱单独申请
https://services.healthtech.dtu.dk/cgi-bin/sw_request

#下载之后解压缩
mkidr rnammer-1.2
tar zxvf rnammer-1.2.src.tar.gz -C rnammer-1.2

#修改rnammer程序路径
$INSTALL_PATH
$HMMSEARCH_BINARY

#运行程序
mkdir ncrna;cd ncrna;
rnammer -S bac -m tsu,lsu,ssu  -gff MGH78578.gff -f MGH78578_rrna.frn MGH78578.fa 


#2 转运RNA(trna)分析
#安装软件trnascan-se
mamba install -y trnascan-se

#检查默认perl版本
perl ~/miniconda3/bin/tRNAscan-SE
perl ~/miniconda3/bin/tRNAscan-SE -B  -o tRNAScan.out -f tRNAScan.out.structure -m stat.list MGH78578.fasta

#提取序列
perl get_tRNA.pl tRNAScan.out MGH78578.fasta MGH78578_trna.ffn

#查看二级结构
http://rna.tbi.univie.ac.at/forna/

#在线分析
>molecule_name
CGCUUCAUAUAAUCCUAAUGAUAUGGUUUGGGAGUUUCUACCAAGAGCCUUAAACUCUUGAUUAUGAAGUG
...(((((((..((((((.........))))))......).((((((.......))))))..))))))...


#####################################################
#                   八、重复序列分析            #
#####################################################

#1 串联重复序列分析 
#串联重复序列预测
trf MGH78578.fasta  2 7 7 80 10 50 500 -f -d -m

#2 与数据库比对预测重复序列
#安装软件RepeatMasker
#repeatmasker
mamba install -y repeatmasker

#数据库下载，需注册
#https://www.girinst.org/server/RepBase/index.php
RepBaseRepeatMaskerEdition-20181026.tar

#配置数据库
#数据库路径：/ifs1/User/meta/miniconda3/share/RepeatMasker/Libraries/
tar -zxvf RepBaseRepeatMaskerEdition-20181026.tar
mv Libraries/* /ifs1/User/meta/miniconda3/share/RepeatMasker/Libraries/
#RepeatModel
http://www.repeatmasker.org/RepeatModeler/
mamba install -y repeatmodeler
mkdir repeatmasker 
RepeatMasker -pa 2 -species bacteria -q  -html -gff -dir repeatmasker  MGH78578.fasta

#3 从头预测RepeatModeler
#安装软件
mamba install -y repeatmodeler

#软件使用
#下方作为测试，只使用了主要的参数
#建立基因组索引
BuildDatabase -name MGH78578 -engine ncbi MGH78578.fasta 
#从头预测
RepeatModeler -pa 8 -engine ncbi -database MGH78578

#####################################################
#                   九、blast比对            #
#####################################################
#安装blast
mamba install -y blast

#ncbi ftp
ftp.ncbi.nlm.nih.gov/blast/db/

#下载blast数据库
for i in {00..50};do echo "~/.aspera/connect/bin/ascp -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh --overwrite=diff -QTr -l6000m anonftp@ftp.ncbi.nlm.nih.gov:blast/db/nt.${i}.tar.gz  ./ ";done;

#自己构建nt库
~/.aspera/connect/bin/ascp -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh --overwrite=diff -QTr -l6000m anonftp@ftp.ncbi.nlm.nih.gov:blast/db/FASTA/nt.gz  ./

#自己构建数据库
gunzip nt.gz
makeblastdb -in nt -dbtype nucl -parse_seqids -out nt

#nt库比对
blastn -db /ifs1/MetaDatabase/nt_20200716/nt -query assembly.fasta -out blast.out  -outfmt 6 -evalue 1e-5 -num_threads 12


#质粒数据库
ftp://ftp.ncbi.nih.gov/refseq/release/plasmid/  #质粒

#下载质粒数据库
alias 'aspera=~/.aspera/connect/bin/ascp -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh -QTr -k 1 -l 50m'
#批量下载
for i in {1..10};do echo aspera anonftp@ftp.ncbi.nlm.nih.gov:refseq/release/plasmid/plasmid.${i}.1.genomic.fna.gz  ./;done;
#与质粒库进行blast比对
blastn -query assembly.fasta -db ../data/plasmid/plasmid.fa -out plasmid18.txt -outfmt 6 -html -evalue 1e-5 -num_threads 12
#提取id
cat plasmid.txt | awk '{if ($3 >80 && $4 > 100 ) print $1}' | sort | uniq >id.list
#提取序列
seqkit grep -r -f id.list assembly.fasta -w 0 >plasmid_4.fa

#筛选共有基因
#建立索引
makeblastdb -in ref.faa -dbtype prot -parse_seqids -out ref.faa
#blastp比对
blastp -query query.faa -out blast.out -db ref.faa -outfmt 6 -evalue 1e-5
#筛选id
cat blast.out |awk '{if ($3 >=80 && $4 >=100 ) print $2}' | sort  |uniq >id.list
#提取序列
seqkit grep -r -f id.list ref.faa

#安装diamond
#bioconda安装
mamba install -y diamond

#现在预编译程序 
wget http://github.com/bbuchfink/diamond/releases/download/v2.0.13/diamond-linux64.tar.gz
tar xzf diamond-linux64.tar.gz

#blastp比对
#建立索引
diamond makedb --in ref.faa --db ref
diamond比对
diamond blastp -q query.faa -d ref -o blastp.txt -p 12 -f 6

#物种鉴定
#检查数据库版本
diamond dbinfo -d /ifs1/Database/nr_diamond/nr.dmnd


#下载物种分类信息
~/.aspera/connect/bin/ascp -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh --overwrite=diff -QTr -l6000m anonftp@ftp.ncbi.nlm.nih.gov:pub/taxonomy/accession2taxid/prot.accession2taxid.FULL.gz .
~/.aspera/connect/bin/ascp -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh --overwrite=diff -QTr -l6000m anonftp@ftp.ncbi.nlm.nih.gov:pub/taxonomy/taxdmp.zip ./
diamond makedb --in ref.faa --db ref

#diamond比对
diamond blastx -q P15.fastq.gz --db /ifs1/MetaDatabase/diamond_20210825/nr -o P15  -p 12 -f 100

#可视化
diamond view -a blastx.daa