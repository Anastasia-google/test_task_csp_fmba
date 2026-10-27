# FP_SNPs: восстановление референсных аллелей (REF/ALT)

## 1. Назначение

Скрипт `add_ref_allele.py` преобразует файл со списком FingerPrinting SNPs (FP SNPs) из формата `#CHROM POS ID allele1 allele2` в формат `#CHROM POS ID REF ALT`.

Исходные данные — файл `FP_SNPs.txt` из архива программы GRAF версии 2.4 (NCBI). Файл содержит 10 000 аутосомных FP SNPs и 1 000 FP SNPs с хромосомы X 
(в GRAF она обозначена номером 23). По условию задачи варианты с X-хромосомы удаляются.

Скрипт просматривает входной файл формата `#CHROM POS ID allele1 allele2` и обращается к референсным данным генома GRCh38.d1.vd1, какое основание на данной
позиции является референсным, а какое альтернативным. Например allele 1 = C, allele 2 = T во входном файле, в референсном геноме указано основание T, значит,
в выходном файле будет REF = T, ALT = C. Если указана некорректная позиция аллеля, нет fasta-файла референсного генома или ни один из аллелей исходного файла
не соответствует референсному основанию, скрипт выведет сообщение на экране и лог-файле.

## 2. Требования к окружению

Скрипт работает **внутри Docker-контейнера**, собранного из `Dockerfile` этого репозитория. На хосте не требуется ни Python, ни `pysam`, ни `samtools`.

В образе установлены:
- Python 3 и `pysam==0.24.1`;
- `samtools`, `bcftools`, `vcftools` (для подготовки референса и сопутствующих задач - подробнее в README.md, который был создан при выполнении задания №2);
- сам скрипт `add_ref_allele.py` в `/usr/local/bin/`.

Референсный геном GRCh38.d1.vd1 должен лежать **на хосте** в `/mnt/data/ref/GRCh38.d1.vd1_mainChr/sepChrs/` и пробрасываться в контейнер
в `/ref/GRCh38.d1.vd1_mainChr/sepChrs/`. В папке должны быть файлы `chr[1-22,M,X,Y].fa` и `chr[1-22,M,X,Y].fa.fai`.

## 3. Предподготовка входного файла

Предподготовка **не является частью скрипта**. Она приводит исходный `FP_SNPs.txt` к формату `#CHROM  POS  ID  allele1  allele2`.
Предподготовка выполнялась в локальной директории */home/anastasssikk/test_task_cspfmba* - при необходимости можно заменить в последующем коде
на название необходимой директории.

Переходим в необходимую директорию, скачиваем архив GRAF 2.4 с официального сайта NCBI, распаковываем, переходим в папку
*/home/anastasssikk/test_task_cspfmba/data* (появляется при распаковке архива) и смотрим первые строки файла `FP_SNPs.txt`

```bash
cd test_task_cspfmba/
wget -O GrafPkg.tar.gz "http://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/GetZip.cgi?zip_name=GRAF_files.zip"
tar -xzf GrafPkg.tar.gz
cd data/
head FP_SNPs.txt
```

```
rs#     chromosome      GB37_position   GB38_position   allele1 allele2
2887286 1       1156131 1220751 C       T
6685064 1       1211292 1275912 T       C
2840528 1       2283896 2352457 A       G
3890745 1       2553624 2622185 C       T
1798246 1       3080855 3164291 C       A
1181875 1       3681831 3765267 C       T
6663840 1       3743319 3826755 A       G
4233262 1       4236982 4176922 C       T
693734  1       4364226 4304166 C       T
```

Колонки:
1. `rs#` — rsID без префикса `rs`;
2. `chromosome` — номер хромосомы (1–22, 23 = X);
3. `GB37_position` — координата в GRCh37 (не нужна, удаляется);
4. `GB38_position` — координата в GRCh38;
5. `allele1` — первый аллель;
6. `allele2` — второй аллель.

Преобразовываем файл с помощью awk и записываем выход в файл `FP_SNPs_10k_GB38_twoAllelsFormat.tsv`.

```bash
awk -F'\t' 'BEGIN{OFS="\t"}
    # заголовок заменить на #CHROM POS ID allele1 allele2
    NR==1 {print "#CHROM","POS","ID","allele1","allele2"; next}
    # пропустить строки с хромосомой 23 (X)
    $2 == 23 {next}
    {
        # добавить префикс "chr" к номеру хромосомы
        chr = "chr" $2;
        # добавить префикс "rs" к rsID
        id  = "rs" $1;
        # вывести номер хромосомы, координаты GRCh38 (колонка 4), rsID и аллели 1 и 2 (5 и 6 колонки, соответственно)
        print chr, $4, id, $5, $6;
    }' FP_SNPs.txt > FP_SNPs_10k_GB38_twoAllelsFormat.tsv
```

Проверяем полученный файл.

```bash
head FP_SNPs_10k_GB38_twoAllelsFormat.tsv
```

```
#CHROM  POS     ID      allele1 allele2
chr1    1220751 rs2887286       C       T
chr1    1275912 rs6685064       T       C
chr1    2352457 rs2840528       A       G
chr1    2622185 rs3890745       C       T
chr1    3164291 rs1798246       C       A
chr1    3765267 rs1181875       C       T
chr1    3826755 rs6663840       A       G
chr1    4176922 rs4233262       C       T
chr1    4304166 rs693734        C       T
```

## 4. Описание скрипта `add_ref_allele.py`

### 4.1. Аргументы командной строки
Справку можно получить при выполнении в docker-контейнере команды `python3 /usr/local/bin/add_ref_allele.py -h`.

|Аргумент|Описание|
|---|---|
|--input|Обязателен. Входной TSV-файл в формате `#CHROM POS ID allele1 allele2`|
|--output|Обязателен. Выходной TSV-файл в формате `#CHROM POS ID REF ALT`|
|--ref-dir|Обязателен. Директория с файлами референсного генома, разбитого на хромосомы chr*.fa и chr*.fa.fai|
|--log|Необязателен. Путь к лог-файлу (по умолчанию `add_ref_allele.log`)|

### 4.2. Логика работы

1. Проверяет существование входного файла и каталога референса. При ошибке пишет понятное сообщение в лог и завершается с кодом 1.

2. Проверяет заголовок входного файла: он должен быть ровно `#CHROM POS ID allele1 allele2`.
   Если заголовок другой — сообщает, что ожидалось и что получено, и завершается.

3. Для каждой строки данных:

    - извлекает хромосому, позицию, ID и два аллеля;
    - приводит аллели к верхнему регистру;
    - открывает FASTA соответствующей хромосомы (с кэшированием — один раз на хромосому);
    - извлекает основание в позиции `POS` (с учётом того, что в файле координаты 1-based, а pysam использует 0-based);
    - сравнивает основание с `allele1` и `allele2`;
    - если совпало с `allele1` — `REF = allele1`, `ALT = allele2`;
    - если совпало с `allele2` — `REF = allele2`, `ALT = allele1`;
    - если не совпало ни с одним — пишет WARNING с указанием позиции и пропускает строку;
    - записывает результат в выходной файл.

4. Ведёт лог с временными метками: старт, каждые 1000 записанных строк, итоговая статистика, финиш.

5. Итоговая строка лога содержит разбивку пропусков по причинам: нет FASTA, некорректная позиция, несовпадение с референсом.


### 4.3. Обработка ошибок

* Нет входного файла → `ERROR: Входной файл не найден: ...`, выход с кодом 1.
* Нет каталога референса → `ERROR: Каталог референса не найден: ...`, выход с кодом 1.
* Неверный заголовок → `ERROR: Неверный заголовок. Ожидалось: ...`, выход с кодом 1.
* Строка с неверным числом колонок → `WARNING` и пропуск строки.
* Позиция не целое число → `WARNING` и пропуск строки.
* Нет FASTA для хромосомы → `WARNING` и пропуск строки.
* Ни один аллель не совпал с референсом → `WARNING` с координатой и пропуск строки.

## 5. Запуск скрипта

1. Клонируйте данный репозиторий `git clone https://github.com/Anastasssikk/test_task_cspfmba.git`
2. По умолчанию создаётся папка `test_task_cspfmba/`, перейдите в неё `cd test_task_cspfmba/`
3. Сборка образа `docker build -t bio-tools:latest .`
4. Запуск контейнера. Монтируется корень репозитория в `/data` и референс — в `/ref/GRCh38.d1.vd1_mainChr/sepChrs/`
   ```bash
   docker run -it --rm \
   -v /home/anastasssikk/test_task_cspfmba:/data \
   -v /mnt/data/ref/GRCh38.d1.vd1_mainChr/sepChrs:/ref/GRCh38.d1.vd1_mainChr/sepChrs \
   bio-tools:latest
   ```
5. Запуск скрипта внутри контейнера.
   ```bash
   python3 /usr/local/bin/add_ref_allele.py \
   --input /data/FP_SNPs_10k_GB38_twoAllelsFormat.tsv \
   --output /data/FP_SNPs_10k_GB38_with_REF_ALT.tsv \
   --ref-dir /ref/GRCh38.d1.vd1_mainChr/sepChrs \
   --log /data/add_ref_allele.log
   ```

## 6. Результат работы скрипта

```
Всего строк: 10000
Записано: 9991
Пропущено: нет FASTA=0, некорректная позиция=0, несовпадение с референсом=9
```

То есть 9 991 из 10 000 SNP успешно сконвертированы, 9 пропущены по причине «ни один аллель не совпал с референсом».
Пропусков из-за отсутствия FASTA или некорректных позиций нет.

|Строка|Позиция|Аллели (allele1, allele2)|Референс в GRCh38|
|---|---|---|---|
|420|chr1:145899155|T, C|G|
|5967|chr10:47420743|T, C|G|
|5968|chr10:47320207|C, A|T|
|5969|chr10:47268342|T, C|A|
|5970|chr10:47099482|G, A|T|
|5977|chr10:46031829|A, G|C|
|6508|chr11:54678670|C, T|G|
|8062|chr15:22907410|T, C|G|
|8063|chr15:22832212|T, C|A|

Во всех 9 случаях ни один из двух аллелей GRAF не совпал с основанием в GRCh38.
Скрипт корректно обнаруживает такие случаи и пропускает их с предупреждением, не прерывая работу.

Выходной файл `FP_SNPs_10k_GB38_with_REF_ALT.tsv` имеет следующий формат.

```
#CHROM  POS  ID  REF  ALT
chr1    1220751  rs2887286  T  C
chr1    1275912  rs6685064  C  T
chr1    2352457  rs2840528  A  G
```

* `#CHROM` — хромосома с префиксом chr;
* `POS` — координата в GRCh38 (1-based);
* `ID` — rsID с префиксом rs;
* `REF` — референсный аллель (совпал с основанием в GRCh38);
* `ALT` — альтернативный аллель (другой аллель из входного файла).

Количество строк: 9 992 (9 991 SNP + заголовок).

## 7. Примечание о референсном геноме
Файлы референсного генома не входят в репозиторий. Предполагается, что они лежат на хосте в /mnt/data/ref/GRCh38.d1.vd1_mainChr/sepChrs/ и пробрасываются
в контейнер в /ref/GRCh38.d1.vd1_mainChr/sepChrs/. В папке должны быть файлы chr[1-22,M,X,Y].fa и chr[1-22,M,X,Y].fa.fai.

Референсный геном был скачан не по данной в задании ссылке (https://gdc.cancer.gov/about-data/data-harmonization-and-generation/gdc-reference-files)
в связи с проблемами с доступом. Референсный геном был скачан по альтернативному пути по ссылке
https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG/technical/reference_genome/GRCh38/GRCh38.d1.vd1.fa

Подготовка референсного генома:

```bash
# Скачать референс GRCh38.d1.vd1
wget https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG/technical/reference_genome/GRCh38/GRCh38.d1.vd1.fa

# Создать целевую папку на хосте
sudo mkdir -p /mnt/data/ref/GRCh38.d1.vd1_mainChr/sepChrs
sudo chown -R $USER:$USER /mnt/data

# Скопировать скачанный FASTA в целевую папку
cp ./GRCh38.d1.vd1.fa /mnt/data/ref/GRCh38.d1.vd1_mainChr/sepChrs/

# Запустить контейнер, пробросив папку с референсом
docker run -it --rm \
  -v /mnt/data/ref/GRCh38.d1.vd1_mainChr/sepChrs:/ref/GRCh38.d1.vd1_mainChr/sepChrs \
  bio-tools:latest

# Внутри контейнера: проиндексировать исходный файл
cd /ref/GRCh38.d1.vd1_mainChr/sepChrs
samtools faidx GRCh38.d1.vd1.fa

# Внутри контейнера: разбить на отдельные хромосомы и проиндексировать каждую
for chr in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM; do
    samtools faidx GRCh38.d1.vd1.fa "$chr" > "${chr}.fa"
    samtools faidx "${chr}.fa"
done

# Внутри контейнера: удалить исходный большой файл и его индекс
rm GRCh38.d1.vd1.fa GRCh38.d1.vd1.fa.fai

# Выйти из контейнера
exit
```
