
ע������##�ô������ע��˵��

##part1

##д�Լ�������ƴ��

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

chmod u=rwx��g=rx��o=rx pop.test.vcf


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
##��ʾ����Ŀ¼���ļ�
tree -L 2 

##ֻ��ʾ����Ŀ¼�µ�Ŀ¼��Ϣ

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


##�����������ض����ַ�

grep "##" pop.test.vcf  

grep -c "##" pop.test.vcf  

grep -n "##" pop.test.vcf 

grep -i "chr" pop.test.vcf |less

grep -v "chr" pop.test.vcf |less



##�����ļ�

find ./ -name "*.sh"
find ./ -type d -name "*d*"


##�ļ��и�

cut -f 1 bsa.result 

cut -f 1-3 bsa.result 

cut -f 2 bsa.result >1.result 

cut -f 3 bsa.result >2.result


##�ļ��ϲ�

paste 1.result 2.result

paste -d ":" 1.result 2.result

paste -s 1.result 2.result 


##��������

sort -n ./demo/num.txt

sort -r -n ./demo/num.txt  

sort -u ./demo/num.txt 

##ȥ������

uniq ./demo/num.txt
uniq -c2 ./demo/num.txt 
uniq -d ./demo/num.txt 
uniq -u ./demo/num.txt 

cat ./demo/num.txt |uniq
cat ./demo/num.txt |sort |uniq
cat ./demo/num.txt |uniq -d 


sort ./demo/num.txt |uniq

sort -n ./demo/num.txt |uniq


##ͳ�����������

wc -l pop.test.vcf

split -l 50 1.result

split -n 3 1.result

echo "100/100"|bc 

##�ļ������������

�����tar -xvf FileName.tar
�����tar -cvf FileName.tar DirName
���ӣ�

tar -cvf mengluo20181112.tar mengluo

tar -xvf mengluo20181112.tar   

tar -zcvf mengluo20181112.tar.gz mengluo

tar -zxvf mengluo20181112.tar.gz 

tar -jcvf mengluo20181112.tar.bz2 mengluo

tar -jxvf mengluo20181112.tar.bz2 


##�߽�����

grep '@chr' region.threshold.vcf.total|sed 's/\t/,/g' >table.csv

less RAxML_info.pop|grep ": Time"|cut -d " " -f3|awk '{s+=$l}END{print s}'



