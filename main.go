package main

import (
	"bufio"   // Для построчного чтения файла
	"fmt"     // Для вывода сообщений в консоль
	"os"      // Для работы с файлами и открытия FASTA
	"strings" // Для работы со строками и преобразования в верхний регистр
)

// dnaStats хранит статистику по нуклеотидам ДНК.
type dnaStats struct {
	total int // Общее количество валидных букв
	a     int // Сколько букв A
	t     int // Сколько букв T
	c     int // Сколько букв C
	g     int // Сколько букв G
}

// Таблица кодонов (генетический код)
// Полная таблица генетического кода (мРНК → аминокислота)
var codonTable = map[string]string{
	// Фенилаланин
	"UUU": "F", "UUC": "F",

	// Лейцин
	"UUA": "L", "UUG": "L",
	"CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L",

	// Изолейцин
	"AUU": "I", "AUC": "I", "AUA": "I",

	// Метионин (старт)
	"AUG": "M",

	// Валин
	"GUU": "V", "GUC": "V", "GUA": "V", "GUG": "V",

	// Серин
	"UCU": "S", "UCC": "S", "UCA": "S", "UCG": "S",
	"AGU": "S", "AGC": "S",

	// Пролин
	"CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P",

	// Треонин
	"ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T",

	// Аланин
	"GCU": "A", "GCC": "A", "GCA": "A", "GCG": "A",

	// Тирозин
	"UAU": "Y", "UAC": "Y",

	// Стоп-кодоны
	"UAA": "Stop",
	"UAG": "Stop",
	"UGA": "Stop",

	// Гистидин
	"CAU": "H", "CAC": "H",

	// Глутамин
	"CAA": "Q", "CAG": "Q",

	// Аспарагин
	"AAU": "N", "AAC": "N",

	// Лизин
	"AAA": "K", "AAG": "K",

	// Аспарагиновая кислота
	"GAU": "D", "GAC": "D",

	// Глутаминовая кислота
	"GAA": "E", "GAG": "E",

	// Цистеин
	"UGU": "C", "UGC": "C",

	// Триптофан
	"UGG": "W",

	// Аргинин
	"CGU": "R", "CGC": "R", "CGA": "R", "CGG": "R",
	"AGA": "R", "AGG": "R",

	// Глицин
	"GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

func formatProtein(protein string) string {
	if len(protein) == 0 {
		return ""
	}

	var b strings.Builder
	for i, aa := range protein {
		if i > 0 {
			b.WriteString(" ➔ ")
		}
		b.WriteRune(aa)
	}
	return b.String()
}
func analyzeProtein(protein string) {
	if strings.HasPrefix(protein, "M") {
		fmt.Println("✅ Белок валиден (начинается с Метионина)")
	} else {
		fmt.Println("⚠️ Ошибка: отсутствует старт-кодон Метионин!")
	}

}
func translate(rna string) string {
	protein := ""
	// Шагаем по строке РНК строго по 3 буквы
	for i := 0; i < len(rna)-2; i += 3 {
		codon := rna[i : i+3]

		// Ищем аминокислоту в таблице
		amino, exists := codonTable[codon]
		if !exists {
			protein += "?" // Если кодон неизвестен
			continue
		}

		// Если наткнулись на Стоп-кодон — мгновенно останавливаемся!
		if amino == "Stop" {
			break
		}

		protein += amino
	}
	return protein
}

// findORF ищет первый старт-кодон ATG и затем идёт по кодонам по 3 нуклеотида
// до первого стоп-кодона TAA, TAG или TGA.
// Возвращает: (кодирующая последовательность, начало_в_1-based, конец_в_1-based).
func findORF(dna string) (string, int, int) {
	cleanDNA := strings.ToUpper(dna)

	for i := 0; i+2 < len(cleanDNA); i++ {
		if cleanDNA[i] == 'A' && cleanDNA[i+1] == 'T' && cleanDNA[i+2] == 'G' {
			start := i + 1 // 1-based индекс старт-кодона

			for j := i + 3; j+2 < len(cleanDNA); j += 3 {
				codon := cleanDNA[j : j+3]
				if codon == "TAA" || codon == "TAG" || codon == "TGA" {
					return cleanDNA[i : j+3], start, j + 3 // j+3 — 1-based индекс последнего нуклеотида стоп-кодона
				}
			}

			fmt.Println("[INFO] Старт-кодон найден, но стоп-кодон не обнаружен.")
			return "", 0, 0
		}
	}

	fmt.Println("[INFO] Старт-кодон ATG не найден.")
	return "", 0, 0
}

// getReverseComplement строит обратную комплементарную цепь ДНК.
func getReverseComplement(dna string) string {
	runes := []rune(dna)

	// Переворачиваем строку задом наперёд.
	for i, j := 0, len(runes)-1; i < j; i, j = i+1, j-1 {
		runes[i], runes[j] = runes[j], runes[i]
	}

	// Заменяем каждую букву на комплементарную.
	for i, r := range runes {
		switch r {
		case 'A':
			runes[i] = 'T'
		case 'T':
			runes[i] = 'A'
		case 'C':
			runes[i] = 'G'
		case 'G':
			runes[i] = 'C'
		default:
			runes[i] = r
		}
	}

	return string(runes)
}

// truncateDNA ограничивает строку до нужной длины и добавляет многоточие.
func truncateDNA(seq string, limit int) string {
	if len(seq) <= limit {
		return seq
	}
	return seq[:limit] + "..."
}

// splicing удаляет все строчные буквы (интроны) и оставляет только заглавные (экзоны).
func splicing(dna string) string {
	var b strings.Builder
	for _, r := range dna {
		if r >= 'A' && r <= 'Z' {
			b.WriteRune(r)
		}
	}
	return b.String()
}

func transcribe(dna string) string {
	return strings.ReplaceAll(strings.ReplaceAll(dna, "T", "U"), "t", "u")
}

func removeTransposon(dna string) (string, bool) {
	if !strings.Contains(dna, "GGGGGGCCCCCC") {
		return dna, false
	}
	return strings.ReplaceAll(dna, "GGGGGGCCCCCC", ""), true
}

// countDNAStats считает количество A/T/C/G в строке.
func countDNAStats(dna string) dnaStats {
	stats := dnaStats{}
	for _, r := range strings.ToUpper(dna) {
		switch r {
		case 'A':
			stats.a++
			stats.total++
		case 'T':
			stats.t++
			stats.total++
		case 'C':
			stats.c++
			stats.total++
		case 'G':
			stats.g++
			stats.total++
		}
	}
	return stats
}

func main() {
	// ANSI-коды для цветного вывода
	purple := "\u001b[35m"
	green := "\u001b[32m"
	cyan := "\u001b[36m"
	reset := "\u001b[0m"

	fmt.Printf("%sBIOGO-CORE v1.3%s\n", purple, reset)

	// Открываем файл с последовательностью ДНК
	file, err := os.Open("gfp_bloated.fasta")
	if err != nil {
		fmt.Println("Не удалось открыть файл:", err)
		return
	}
	defer file.Close()

	fmt.Printf("%s[SUCCESS] Файл gfp_bloated.fasta успешно прочитан.%s\n", green, reset)

	// Считываем файл построчно, убирая заголовок и переносы строк
	scanner := bufio.NewScanner(file)
	var dnaSequence strings.Builder
	position := 0

	for scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())
		if len(line) == 0 || line[0] == '>' {
			continue
		}

		cleanLine := strings.ToUpper(line)
		dnaSequence.WriteString(cleanLine)

		for _, r := range cleanLine {
			position++
			switch r {
			case 'A', 'T', 'C', 'G':
				// валидный нуклеотид
			default:
				fmt.Printf("Обнаружен мусорный символ %c на позиции %d\n", r, position)
			}
		}
	}

	if err := scanner.Err(); err != nil {
		fmt.Println("Ошибка при чтении файла:", err)
		return
	}

	cleanDNA := dnaSequence.String()
	cleanDNA, found := removeTransposon(cleanDNA)
	if found {
		fmt.Printf("%s[WARNING] Обнаружен транспозон! Проведена очистка генома.%s\n", "\u001b[33m", reset)
	}
	orfSequence, startIdx, endIdx := findORF(cleanDNA)
	if orfSequence == "" {
		fmt.Println("[INFO] ORF не найден. Завершение программы.")
		return
	}

	fmt.Printf("\n%s=== ОТКРЫТАЯ РАМКА СЧИТЫВАНИЯ (ORF) ===%s\n", purple, reset)
	fmt.Printf("%s• Начало ORF:%s %d\n", green, reset, startIdx)
	fmt.Printf("%s• Конец ORF:%s %d\n", green, reset, endIdx)
	fmt.Printf("%s• Длина ORF:%s %d\n", green, reset, len(orfSequence))

	stats := countDNAStats(orfSequence)
	fmt.Printf("\n%s=== Статистика по найденному гену ===%s\n", purple, reset)
	fmt.Printf("%sA:%s %d\n", green, reset, stats.a)
	fmt.Printf("%sT:%s %d\n", green, reset, stats.t)
	fmt.Printf("%sC:%s %d\n", green, reset, stats.c)
	fmt.Printf("%sG:%s %d\n", green, reset, stats.g)

	if stats.total > 0 {
		gcPercent := (float64(stats.g+stats.c) / float64(stats.total)) * 100
		fmt.Printf("%sGC-состав:%s %.2f%%%s\n", cyan, reset, gcPercent, reset)
	}

	reverseComplement := getReverseComplement(orfSequence)
	first60 := truncateDNA(orfSequence, 60)
	second60 := truncateDNA(reverseComplement, 60)

	fmt.Printf("\n%s=== ДВОЙНАЯ СПИРАЛЬ ORF (Первые 60 bp) ===%s\n", purple, reset)
	fmt.Printf("%s5'- %s -3'%s\n", green, first60, reset)
	fmt.Printf("    %s\n", strings.Repeat("|", len(first60)))
	fmt.Printf("%s3'- %s -5'%s\n", cyan, second60, reset)

	mrna := transcribe(orfSequence)
	fmt.Printf("\n%s=== ТРАНСКРИПЦИЯ ===%s\n", purple, reset)
	fmt.Printf("%sДНК(первые 60):%s %s\n", green, reset, truncateDNA(orfSequence, 60))
	fmt.Printf("%sмРНК(первые 60):%s %s\n", "\u001b[35m", reset, truncateDNA(mrna, 60))
	// Перевод мРНК в белок
	protein := translate(mrna)

	fmt.Printf("\n%s=== ТРАНСЛЯЦИЯ ===%s\n", purple, reset)
	fmt.Printf("%sБелок GFP:%s %s\n", cyan, reset, formatProtein(protein))
	fmt.Printf("%sДлина белка:%s %d аминокислот\n", green, reset, len(protein))
	analyzeProtein(protein)

	testDNA := "ATGgctgtaCAGatcgtaTAA"
	fmt.Printf("\n%s=== СПЛАЙСИНГ (демо) ===%s\n", purple, reset)
	fmt.Printf("%sИсходная ДНК:%s %s\n", "\u001b[31m", reset, testDNA)
	fmt.Printf("%sОчищенная ДНК:%s %s\n", green, reset, splicing(testDNA))
}
