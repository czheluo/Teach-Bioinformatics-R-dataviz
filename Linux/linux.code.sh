
注：两个##好代表的是注释说明

##part1

##写自己中文名拼音

mkdir name 

mkdir -p test test1 test2

mkdir -v trash 

##part2
cd /
cd 
cd -
cd ..

#part3 
##

cp /home/teacher/1.Linux/mengluo/pop.test.vcf ./

cp /home/teacher/1.Linux/mengluo/bsa.result ./

cp -i /home/teacher/1.Linux/mengluo/pop.test.vcf ./

cp -r /home/teacher/1.Linux/mengluo/data /home/teacher/1.Linux/mengluo/demo /home/teacher/1.Linux/mengluo/other ./


##or

ln /home/teacher/1.Linux/mengluo/data/S120.clean.1.fastq.gz

ln -s /home/teacher/1.Linux/mengluo/data/S120.clean.2.fastq.gz

ln -s /home/teacher/1.Linux/mengluo/data/S120.clean.2.fastq.gz 2.fastq.gz

##

chmod 755 test1

chmod u-w,g-w,o-r pop.test.vcf

chmod u=rwx，g=rx，o=rx pop.test.vcf


##part4 to check the file 

mv test trash

mv -f test1 trash

mv -i test2 trash 

##

rm -i S120.clean.2.fastq.gz

rm -f S120.clean.2.fastq.gz

rm -r test2

#part5 tree

tree
tree -d 
##显示两级目录和文件
tree -L 2 

##只显示两级目录下的目录信息

tree -L 2 -d  

tree -P t


##check file 

more pop.test.vcf

less -S pop.test.vcf 

less -N pop.test.vcf 
less -N -S pop.test.vcf 

##|| two command

head -n 20 pop.test.vcf >head20

tail -n 30 pop.test.vcf >tail30

tail -n 30 pop.test.vcf >>head20

head -n -3 pop.test.vcf |tail -n +11

head -n -3 pop.test.vcf |tail -n +11 >head.tail


#cat

cat head20

cat -n head20

cat head20 tail30 

cat head20 tail30 > cat.file


##如果你想查找特定的字符

grep "##" pop.test.vcf  

grep -c "##" pop.test.vcf  

grep -n "##" pop.test.vcf 

grep -i "chr" pop.test.vcf |less

grep -v "chr" pop.test.vcf |less



##查找文件

find ./ -name "*.sh"
find ./ -type d -name "*d*"


##文件切割

cut -f 1 bsa.result 

cut -f 1-3 bsa.result 

cut -f 2 bsa.result >1.result 

cut -f 3 bsa.result >2.result


##文件合并

paste 1.result 2.result

paste -d ":" 1.result 2.result

paste -s 1.result 2.result 


##排序命令

sort -n ./demo/num.txt

sort -r -n ./demo/num.txt  

sort -u ./demo/num.txt 

##去重命令

uniq ./demo/num.txt
uniq -c2 ./demo/num.txt 
uniq -d ./demo/num.txt 
uniq -u ./demo/num.txt 

cat ./demo/num.txt |uniq
cat ./demo/num.txt |sort |uniq
cat ./demo/num.txt |uniq -d 


sort ./demo/num.txt |uniq

sort -n ./demo/num.txt |uniq


##统计与计算命令

wc -l pop.test.vcf

split -l 50 1.result

split -n 3 1.result

echo "100/100"|bc 

##文件打包与解包命令

解包：tar -xvf FileName.tar
打包：tar -cvf FileName.tar DirName
例子：

tar -cvf mengluo20181112.tar mengluo

tar -xvf mengluo20181112.tar   

tar -zcvf mengluo20181112.tar.gz mengluo

tar -zxvf mengluo20181112.tar.gz 

tar -jcvf mengluo20181112.tar.bz2 mengluo

tar -jxvf mengluo20181112.tar.bz2 


##高阶命令

grep '@chr' region.threshold.vcf.total|sed 's/\t/,/g' >table.csv

less RAxML_info.pop|grep ": Time"|cut -d " " -f3|awk '{s+=$l}END{print s}'



