#!/usr/bin/env python3
"""
add_ref_allele.py — восстановление референсного аллеля для FP SNPs.

Вход: TSV с колонками #CHROM, POS, ID, allele1, allele2.
Выход: TSV с колонками #CHROM, POS, ID, REF, ALT.
Референсный геном: GRCh38.d1.vd1, файлы chr*.fa в каталоге --ref-dir.
"""

import argparse
import logging
import os
import sys
from datetime import datetime

try:
    import pysam
except ImportError:
    sys.stderr.write("Ошибка: не установлен pysam. Установите: pip install pysam\n")
    sys.exit(1)

# настройка логгера, вывод текста на и экран и сохранение в лог-файл в заданном формате
def setup_logger(log_path):
    logger = logging.getLogger("add_ref_allele")
    logger.setLevel(logging.DEBUG)
    fmt = logging.Formatter("[%(asctime)s] %(levelname)s: %(message)s",
                            datefmt="%Y-%m-%d %H:%M:%S")
    fh = logging.FileHandler(log_path, mode="w", encoding="utf-8")
    fh.setFormatter(fmt)
    ch = logging.StreamHandler(sys.stdout)
    ch.setFormatter(fmt)
    logger.addHandler(fh)
    logger.addHandler(ch)
    return logger

# парсинг аргументов
def parse_args():
    parser = argparse.ArgumentParser(
        description="Восстановление REF/ALT для файла FP SNPs с использованием GRCh38.d1.vd1.",
        epilog="Пример: python add_ref_allele.py --input in.tsv --output out.tsv "
               "--ref-dir /ref/GRCh38.d1.vd1_mainChr/sepChrs --log run.log"
    )
    parser.add_argument("--input", required=True, help="Входной TSV-файл")
    parser.add_argument("--output", required=True, help="Выходной TSV-файл")
    parser.add_argument("--ref-dir", required=True,
                        help="Каталог с файлами chr*.fa и chr*.fa.fai")
    parser.add_argument("--log", default="add_ref_allele.log",
                        help="Путь к лог-файлу (по умолчанию: add_ref_allele.log)")
    return parser.parse_args()

 # принимает директорию с референсом и номером хромосомы, выдаёт fa-файл и имя последовательности в данном файле
def open_fasta(ref_dir, chrom):
    """Открывает FASTA для хромосомы; возвращает объект или None."""
    # chrom может быть '1', 'chr1', 'X', 'chrX' и т.п.
    candidates = []
    if chrom.startswith("chr"):
        candidates.append(chrom)
        candidates.append(chrom[3:])
    else:
        candidates.append(chrom)
        candidates.append("chr" + chrom)
    seen = set()
    for name in candidates:
        if name in seen:
            continue
        seen.add(name)
        fa_path = os.path.join(ref_dir, name + ".fa")
        if os.path.isfile(fa_path):
            try:
                return pysam.FastaFile(fa_path), name
            except Exception as e:
                raise RuntimeError(f"Не удалось открыть {fa_path}: {e}")
    return None, None

 # непосредственно работа скрипта
def main():
    args = parse_args()
    logger = setup_logger(args.log)
    logger.info("Старт: %s", datetime.now().isoformat(timespec="seconds"))

    if not os.path.isfile(args.input):
        logger.error("Входной файл не найден: %s", args.input)
        sys.exit(1)
    if not os.path.isdir(args.ref_dir):
        logger.error("Каталог референса не найден: %s", args.ref_dir)
        sys.exit(1)

    # Проверка заголовка
    with open(args.input, "r", encoding="utf-8", newline="") as f:
        first = f.readline().rstrip("\r\n")
    expected = ["#CHROM", "POS", "ID", "allele1", "allele2"]
    got = first.split("\t")
    if got != expected:
        logger.error("Неверный заголовок. Ожидалось: %s", "\t".join(expected))
        logger.error("Получено: %s", first)
        sys.exit(1)
    logger.info("Заголовок корректен: %s", first)

    total = 0 # всего строк обработано
    written = 0 # строк записано в выходной файл
    skipped_no_fasta = 0 # строк пропущено из отсутствия fa-файла
    skipped_bad_pos = 0 # строк пропущено из-за некорректной позиции
    skipped_ref_mismatch = 0 # строк пропущено из-за того, что ни один аллель не совпал с референсом
    fasta_cache = {} # словарь-кэш

    with open(args.input, "r", encoding="utf-8", newline="") as fin, \
         open(args.output, "w", encoding="utf-8", newline="") as fout:
        header = fin.readline()
        fout.write("#CHROM\tPOS\tID\tREF\tALT\n")

        for line_no, line in enumerate(fin, start=2):
            line = line.rstrip("\r\n")
            if not line:
                continue
            total += 1
            parts = line.split("\t")
            if len(parts) != 5:
                logger.warning("Строка %d: ожидалось 5 колонок, получено %d — пропуск",
                               line_no, len(parts))
                continue
            chrom, pos_s, snp_id, a1, a2 = parts
            a1 = a1.upper()
            a2 = a2.upper()

            try:
                pos = int(pos_s)
            except ValueError:
                logger.warning("Строка %d: POS не целое число (%s) — пропуск", line_no, pos_s)
                skipped_bad_pos += 1
                continue

            key = chrom
            if key not in fasta_cache:
                try:
                    fasta, resolved = open_fasta(args.ref_dir, chrom)
                except Exception as e:
                    logger.error("Строка %d: %s", line_no, e)
                    fasta, resolved = None, None
                fasta_cache[key] = (fasta, resolved)

            fasta, resolved = fasta_cache[key]
            if fasta is None:
                logger.warning("Строка %d: нет FASTA для хромосомы '%s' — пропуск",
                               line_no, chrom)
                skipped_no_fasta += 1
                continue

            try:
                # pysam.fetch использует 0-based координаты
                ref_base = fasta.fetch(resolved, pos - 1, pos).upper()
            except Exception as e:
                logger.warning("Строка %d: ошибка извлечения (%s:%d): %s — пропуск",
                               line_no, resolved, pos, e)
                skipped_bad_pos += 1
                continue

            if len(ref_base) != 1:
                logger.warning("Строка %d: получено '%s' вместо одного основания — пропуск",
                               line_no, ref_base)
                skipped_bad_pos += 1
                continue

            if ref_base == a1:
                ref, alt = a1, a2
            elif ref_base == a2:
                ref, alt = a2, a1
            else:
                logger.warning("Строка %d: ни один аллель (%s,%s) не совпал с референсом '%s' "
                               "в %s:%d — пропуск", line_no, a1, a2, ref_base, resolved, pos)
                skipped_ref_mismatch += 1
                continue

            fout.write(f"{chrom}\t{pos}\t{snp_id}\t{ref}\t{alt}\n")
            written += 1
            if written % 1000 == 0:
                logger.info("Обработано %d строк, записано %d", total, written)

    logger.info("Готово. Всего строк: %d, записано: %d, пропущено: "
                "нет FASTA=%d, некорректная позиция=%d, несовпадение с референсом=%d",
                total, written, skipped_no_fasta, skipped_bad_pos, skipped_ref_mismatch)
    logger.info("Выходной файл: %s", args.output)
    for f in fasta_cache.values():
        if f[0] is not None:
            f[0].close()
    logger.info("Финиш: %s", datetime.now().isoformat(timespec="seconds"))


if __name__ == "__main__":
    main()