# src/scorer.py

import numpy as np
from dataclasses import dataclass


class BayesianScorer:
    """贝叶斯评分系统"""

    def __init__(self, min_cds_len: int = 50):
        """
        初始化贝叶斯评分器

        Args:
            min_cds_len: 最小CDS长度阈值
        """
        self.min_cds_len = min_cds_len

        # 设置权重
        self.weights = {
            "tis_prob": 0.30,
            "tts_prob": 0.50,
            "kozak": 0.04,
            "cai": 0.04,
            "gc_score": 0,
        }

        # 设置Kozak序列PWM矩阵
        self.kozak_pwm = {
            -6: {"A": 0.22, "C": 0.28, "G": 0.32, "T": 0.18},
            -5: {"A": 0.20, "C": 0.30, "G": 0.30, "T": 0.20},
            -4: {"A": 0.18, "C": 0.32, "G": 0.30, "T": 0.20},
            -3: {"A": 0.25, "C": 0.15, "G": 0.45, "T": 0.15},  # 重要位点
            -2: {"A": 0.20, "C": 0.35, "G": 0.25, "T": 0.20},
            -1: {"A": 0.20, "C": 0.35, "G": 0.25, "T": 0.20},
            0: {"A": 1.00, "C": 0.00, "G": 0.00, "T": 0.00},  # A of ATG
            1: {"A": 0.00, "C": 0.00, "G": 0.00, "T": 1.00},  # T of ATG
            2: {"A": 0.00, "C": 0.00, "G": 1.00, "T": 0.00},  # G of ATG
            3: {"A": 0.20, "C": 0.20, "G": 0.40, "T": 0.20},  # +4位点
        }
        self.max_score = (
            self._calculate_max_score()
        )  # 定义最大可能得分（理想Kozak序列的得分）

        # 设置密码子使用频率
        self.codon_usage = {
            "TTT": 0.45,
            "TTC": 0.55,
            "TTA": 0.07,
            "TTG": 0.13,
            "TCT": 0.18,
            "TCC": 0.22,
            "TCA": 0.15,
            "TCG": 0.06,
            "TAT": 0.43,
            "TAC": 0.57,
            "TAA": 0.28,
            "TAG": 0.20,
            "TGT": 0.45,
            "TGC": 0.55,
            "TGA": 0.52,
            "TGG": 1.00,
            "CTT": 0.13,
            "CTC": 0.20,
            "CTA": 0.07,
            "CTG": 0.41,
            "CCT": 0.28,
            "CCC": 0.33,
            "CCA": 0.27,
            "CCG": 0.11,
            "CAT": 0.41,
            "CAC": 0.59,
            "CAA": 0.25,
            "CAG": 0.75,
            "CGT": 0.08,
            "CGC": 0.19,
            "CGA": 0.11,
            "CGG": 0.21,
            "ATT": 0.36,
            "ATC": 0.48,
            "ATA": 0.16,
            "ATG": 1.00,
            "ACT": 0.24,
            "ACC": 0.36,
            "ACA": 0.28,
            "ACG": 0.12,
            "AAT": 0.46,
            "AAC": 0.54,
            "AAA": 0.42,
            "AAG": 0.58,
            "AGT": 0.15,
            "AGC": 0.24,
            "AGA": 0.20,
            "AGG": 0.20,
            "GTT": 0.18,
            "GTC": 0.24,
            "GTA": 0.11,
            "GTG": 0.47,
            "GCT": 0.26,
            "GCC": 0.40,
            "GCA": 0.23,
            "GCG": 0.11,
            "GAT": 0.46,
            "GAC": 0.54,
            "GAA": 0.42,
            "GAG": 0.58,
            "GGT": 0.16,
            "GGC": 0.34,
            "GGA": 0.25,
            "GGG": 0.25,
        }
        #
        self.genetic_code = {
            "ATA": "I",
            "ATC": "I",
            "ATT": "I",
            "ATG": "M",
            "ACA": "T",
            "ACC": "T",
            "ACG": "T",
            "ACT": "T",
            "AAC": "N",
            "AAT": "N",
            "AAA": "K",
            "AAG": "K",
            "AGC": "S",
            "AGT": "S",
            "AGA": "R",
            "AGG": "R",
            "CTA": "L",
            "CTC": "L",
            "CTG": "L",
            "CTT": "L",
            "CCA": "P",
            "CCC": "P",
            "CCG": "P",
            "CCT": "P",
            "CAC": "H",
            "CAT": "H",
            "CAA": "Q",
            "CAG": "Q",
            "CGA": "R",
            "CGC": "R",
            "CGG": "R",
            "CGT": "R",
            "GTA": "V",
            "GTC": "V",
            "GTG": "V",
            "GTT": "V",
            "GCA": "A",
            "GCC": "A",
            "GCG": "A",
            "GCT": "A",
            "GAC": "D",
            "GAT": "D",
            "GAA": "E",
            "GAG": "E",
            "GGA": "G",
            "GGC": "G",
            "GGG": "G",
            "GGT": "G",
            "TCA": "S",
            "TCC": "S",
            "TCG": "S",
            "TCT": "S",
            "TTC": "F",
            "TTT": "F",
            "TTA": "L",
            "TTG": "L",
            "TAC": "Y",
            "TAT": "Y",
            "TAA": "*",
            "TAG": "*",
            "TGC": "C",
            "TGT": "C",
            "TGA": "*",
            "TGG": "W",
        }

    def _calculate_max_score(self) -> float:
        """计算理想Kozak序列的最大得分（用于归一化）"""
        max_score = 1.0
        for pos in self.kozak_pwm:
            max_base = max(self.kozak_pwm[pos].values())
            max_score *= max_base
        return max_score

    def calculate_kozak_score(self, sequence: str, tis_pos: int) -> float:
        """计算Kozak序列得分，并归一化到0-1范围"""
        if tis_pos < 6 or tis_pos + 4 >= len(sequence):
            return 0.0

        kozak_region = sequence[tis_pos - 6 : tis_pos + 4].upper()
        if len(kozak_region) != 10:
            return 0.0

        # 计算原始得分
        raw_score = 1.0
        for i, base in enumerate(kozak_region):
            pos = i - 6
            if pos in self.kozak_pwm and base in self.kozak_pwm[pos]:
                raw_score *= self.kozak_pwm[pos][base]

        # 归一化到0-1范围
        normalized_score = raw_score / self.max_score
        return normalized_score

    def calculate_kozak_score_(self, sequence: str, tis_pos: int) -> float:
        """计算Kozak序列得分"""
        if tis_pos < 6 or tis_pos + 4 >= len(sequence):
            return 0.0

        kozak_region = sequence[tis_pos - 6 : tis_pos + 4].upper()
        if len(kozak_region) != 10:
            return 0.0

        score = 1.0
        for i, base in enumerate(kozak_region):
            pos = i - 6
            if pos in self.kozak_pwm and base in self.kozak_pwm[pos]:
                score *= self.kozak_pwm[pos][base]

        return score * 10000

    def calculate_cai(self, sequence: str) -> float:
        """计算密码子适应指数(CAI)"""
        if len(sequence) < 3:
            return 0.0

        w_values = []
        for i in range(0, len(sequence) - 2, 3):
            codon = sequence[i : i + 3].upper()
            if codon in self.codon_usage:
                w_values.append(np.log(self.codon_usage[codon]))

        if not w_values:
            return 0.0

        cai = np.exp(np.mean(w_values))
        return min(cai, 1.0)

    def calculate_gc_score(self, sequence: str) -> float:
        """计算GC含量得分"""
        if not sequence:
            return 0.0

        gc_count = sequence.count("G") + sequence.count("C")
        gc_content = gc_count / len(sequence)
        gc_score = 2 * np.exp(-0.5 * ((gc_content - 0.42) / 0.2) ** 2) - 1
        return gc_score

    def calculate_integrated_score(
        self,
        tis_prob: float,
        tts_prob: float,
        kozak_score: float,
        cai_score: float,
        gc_score: float,
    ) -> float:
        """计算整合得分"""
        return (
            self.weights["tis_prob"] * tis_prob
            + self.weights["tts_prob"] * tts_prob
            + self.weights["kozak"] * kozak_score
            + self.weights["cai"] * cai_score
            + self.weights["gc_score"] * gc_score
        )
